"""H/G specialist for ambiguous residues in the eight-state baseline."""

from dataclasses import dataclass
import hashlib
from pathlib import Path

import numpy as np
import torch
from torch import nn

from .baseline import PredictionStore, PredictionWriter, dataset_entries
from .constants import G_IDX, H_IDX
from .dataset import DEFAULT_MODEL, ProteinDataset, ROOT

RADIUS = 2
MAX_GAP = 0.80
FOLDS = {"oof": [0, 1, 2, 3, 4], "val": [-1], "cb513": [-2]}
KINDS = {"oof": "out_of_fold", "val": "validation", "cb513": "cb513"}


class HGSpecialist(nn.Sequential):
    """Estimate P(G | H or G) from the central embedding and five Q8 vectors."""

    def __init__(self, embedding_dim=1280):
        super().__init__(
            nn.Linear(embedding_dim + 8 * (2 * RADIUS + 1), 64),
            nn.ReLU(), nn.Dropout(0.2), nn.Linear(64, 1), nn.Flatten(start_dim=0),
        )


def candidate_positions(probs):
    g, h = probs[:, G_IDX], probs[:, H_IDX]
    top_two = np.minimum(g, h) > probs[:, 2:].max(axis=1)
    return np.flatnonzero(top_two & (g != h) & (np.abs(g - h) <= MAX_GAP))


def features(embedding, probs, positions):
    # Pad each chain separately. All-zero rows indicate missing neighbors.
    padded = np.pad(probs, ((RADIUS, RADIUS), (0, 0)))
    window = padded[positions[:, None] + np.arange(2 * RADIUS + 1)]
    return np.concatenate((embedding[positions], window.reshape(len(positions), 8 * (2 * RADIUS + 1))), axis=1)


def confusion(truth, predicted):
    return np.bincount(truth.astype(np.int64) * 8 + predicted, minlength=64).reshape(8, 8)


@dataclass
class Candidates:
    store: PredictionStore
    entries: list
    positions: np.ndarray  # global row numbers in this store
    inputs: torch.Tensor
    truth: np.ndarray  # -1 where masked
    base_labels: np.ndarray
    gaps: np.ndarray
    baseline_confusion: np.ndarray


def read_split(kind, baseline_dir):
    ds = ProteinDataset(FOLDS[kind])
    store = PredictionStore(Path(baseline_dir) / f"baseline_{kind}.npy")
    digest = hashlib.sha256((ROOT / "splits/splits.json").read_bytes()).hexdigest()
    entries = dataset_entries(ds)
    if (store.metadata.get("kind") != KINDS[kind] or store.metadata.get("split_sha256") != digest
            or store.metadata.get("esm_model") != DEFAULT_MODEL or set(store.uids) != {e["uid"] for e in entries}):
        raise ValueError(f"{kind}: baseline source, split, or protein IDs do not match")
    matrix = np.zeros((8, 8), dtype=np.int64)
    positions, inputs, truth, labels, gaps = [], [], [], [], []
    for record, fold in zip(ds.records, ds.fold_ids):
        i = store.row[record.uid]
        probs = store.get(record.uid)
        if (store.fold_ids[i] != fold or len(probs) != record.length or not np.isfinite(probs).all()
                or np.any(probs < 0) or not np.allclose(probs.sum(1), 1, atol=1e-4)):
            raise ValueError(f"{record.uid}: invalid fold, length, or Q8 probabilities")
        pred = probs.argmax(1)
        valid = record.valid_mask
        matrix += confusion(record.q8_label[valid], pred[valid])
        rows = candidate_positions(probs)
        if not len(rows):
            continue
        emb = ds.stores[record.source].get(record.uid)
        x = features(emb, probs, rows)
        if not np.isfinite(x).all():
            raise ValueError(f"Non-finite specialist features for {record.uid}")
        inputs.append(x)
        positions.append(rows + store.offsets[i])
        truth.append(np.where(valid[rows], record.q8_label[rows], -1))
        labels.append(pred[rows])
        gaps.append(np.abs(probs[rows, G_IDX] - probs[rows, H_IDX]))
    x = np.concatenate(inputs) if inputs else np.empty((0, ds.dim + 8 * (2 * RADIUS + 1)), np.float32)
    def joined(parts, dtype):
        return np.concatenate(parts).astype(dtype, copy=False) if parts else np.empty(0, dtype)
    return Candidates(store, entries, joined(positions, np.int64), torch.from_numpy(x),
                      joined(truth, np.int64), joined(labels, np.int64), joined(gaps, np.float32), matrix)


@torch.inference_mode()
def predict_g(model, inputs):
    model.eval()
    if len(inputs) == 0:
        return np.empty(0, np.float32)
    return torch.cat([torch.sigmoid(model(batch)) for batch in inputs.split(2048)]).numpy()


def override_mask(base, gaps, p_g, rule):
    if not rule["enabled"]:
        return np.zeros(len(base), dtype=bool)
    return ((gaps <= rule["max_gap"])
            & (((base == H_IDX) & (p_g > rule["g_threshold"]))
               | ((base == G_IDX) & (p_g < 1 - rule["h_threshold"]))))


def summarize(matrix):
    result = {"valid_residues": int(matrix.sum()), "Q8": float(matrix.trace() / matrix.sum())}
    for name, i in (("G", G_IDX), ("H", H_IDX)):
        tp = int(matrix[i, i])
        actual, predicted = int(matrix[i].sum()), int(matrix[:, i].sum())
        result.update({f"{name}_precision": tp / max(predicted, 1),
                       f"{name}_recall": tp / max(actual, 1),
                       f"{name}_F1": 2 * tp / max(actual + predicted, 1)})
    return result


def audit(data, flip):
    mask = flip & (data.truth >= 0)
    t, old = data.truth[mask], data.base_labels[mask]
    new = np.where(old == H_IDX, G_IDX, H_IDX)
    matrix = data.baseline_confusion - confusion(t, old) + confusion(t, new)
    return {**summarize(matrix), "valid_overrides": int(mask.sum()), "corrected": int((new == t).sum()),
            "harmed": int((old == t).sum()), "other_state_truth": int(((t != G_IDX) & (t != H_IDX)).sum()),
            "G_as_H": int(matrix[G_IDX, H_IDX]), "H_as_G": int(matrix[H_IDX, G_IDX])}


def select_rule(validation, p_g):
    if validation.store.metadata["kind"] != "validation":
        raise ValueError("Choose routing thresholds on validation only")
    rules = [{"enabled": False}]
    for gap in (0.15, 0.35, 0.55, MAX_GAP):
        for g in (0.5, 0.6, 0.7):
            for h in (0.5, 0.6, 0.7):
                rules.append({"enabled": True, "max_gap": gap, "g_threshold": g, "h_threshold": h})
    trials = []
    for rule in rules:
        flip = override_mask(validation.base_labels, validation.gaps, p_g, rule)
        trials.append({"rule": rule, "metrics": audit(validation, flip)})
    eligible = [row for row in trials if row["metrics"]["corrected"] >= row["metrics"]["harmed"]]
    best = max(eligible, key=lambda r: (r["metrics"]["G_F1"], r["metrics"]["Q8"], -r["metrics"]["valid_overrides"]))
    return best["rule"], trials


def export_combined(data, flip, path, metadata):
    # A hard override swaps the original top two scores. Other scores and total mass
    # are preserved. These are decision scores, not calibrated specialist probabilities.
    swap = np.zeros(len(data.store.array), dtype=bool)
    swap[data.positions[flip]] = True
    writer = PredictionWriter(path, data.entries, metadata)
    for entry in data.entries:
        uid = entry["uid"]
        i = data.store.row[uid]
        changed = swap[data.store.offsets[i]:data.store.offsets[i + 1]]
        probs = data.store.get(uid).copy()
        probs[changed, :2] = probs[changed, 1::-1]
        writer.write_batch({"uid": [uid], "lengths": [entry["length"]]}, probs[None])
    writer.finish()
