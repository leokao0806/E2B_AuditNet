"""BiLSTM Q8 baseline and prediction export helpers."""

from __future__ import annotations

import gc
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence

from .constants import G_IDX, H_IDX, IGNORE_INDEX, Q8_TO_Q3


@dataclass(frozen=True)
class BaselineConfig:
    input_dim: int = 1280
    hidden_dim: int = 256
    dropout: float = 0.2


class BaselineModel(nn.Module):
    """One bidirectional LSTM layer over frozen, per-residue ESM-2 vectors."""

    def __init__(self, config: BaselineConfig):
        super().__init__()
        self.config = config
        self.lstm = nn.LSTM(config.input_dim, config.hidden_dim, batch_first=True, bidirectional=True)
        self.dropout = nn.Dropout(config.dropout)
        self.classifier = nn.Linear(2 * config.hidden_dim, 8)

    def forward(self, esm_emb: torch.Tensor, pad_mask: torch.Tensor) -> torch.Tensor:
        if esm_emb.ndim != 3 or pad_mask.shape != esm_emb.shape[:2]:
            raise ValueError("Expected embeddings [B, L, D] and pad_mask [B, L]")
        if esm_emb.shape[-1] != self.config.input_dim:
            raise ValueError(f"Expected embedding dimension {self.config.input_dim}, got {esm_emb.shape[-1]}")
        lengths = pad_mask.sum(dim=1).cpu()
        if (lengths == 0).any():
            raise ValueError("A batch contains an empty protein")
        packed = pack_padded_sequence(esm_emb, lengths, batch_first=True, enforce_sorted=False)
        encoded, _ = self.lstm(packed)
        encoded, _ = pad_packed_sequence(encoded, batch_first=True, total_length=esm_emb.shape[1])
        return self.classifier(self.dropout(encoded))


def masked_loss(logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    """Mean cross-entropy over ordered, labelled residues only."""
    return F.cross_entropy(logits.transpose(1, 2), target, ignore_index=IGNORE_INDEX)


def confusion_metrics(confusion: torch.Tensor) -> dict:
    """Metrics from a Q8 confusion matrix with truth on rows."""
    c = confusion.cpu().numpy()
    n = int(c.sum())
    if n == 0:
        raise ValueError("No valid residues to score")
    q3 = np.zeros((3, 3), dtype=np.int64)
    for truth in range(8):
        for pred in range(8):
            q3[Q8_TO_Q3[truth], Q8_TO_Q3[pred]] += c[truth, pred]
    g_total = int(c[G_IDX].sum())
    h_total = int(c[H_IDX].sum())

    def percent(numerator, denominator):
        return float(round(100 * numerator / denominator, 2)) if denominator else None

    return {
        "residues": n,
        "Q8": percent(np.trace(c), n),
        "Q3": percent(np.trace(q3), n),
        "G_recall": percent(c[G_IDX, G_IDX], g_total),
        "G_predicted_as_H": percent(c[G_IDX, H_IDX], g_total),
        "H_recall": percent(c[H_IDX, H_IDX], h_total),
    }


def _update_confusion(confusion: torch.Tensor, batch: dict, logits: torch.Tensor) -> None:
    valid = batch["valid_mask"]
    truth = batch["q8_label"][valid]
    pred = logits.argmax(dim=-1).cpu()[valid]
    confusion += torch.bincount(truth * 8 + pred, minlength=64).reshape(8, 8)


def run_epoch(model: BaselineModel, loader, device: torch.device, optimizer=None, scaler=None) -> dict:
    """Train or score one pass; loss is weighted by valid residue count."""
    training = optimizer is not None
    model.train(training)
    loss_sum, residue_count = 0.0, 0
    confusion = torch.zeros((8, 8), dtype=torch.long)
    for batch in loader:
        count = int(batch["valid_mask"].sum())
        if count == 0:
            continue
        x = batch["esm_emb"].to(device, non_blocking=True)
        pad_mask = batch["pad_mask"].to(device, non_blocking=True)
        target = batch["target"].to(device, non_blocking=True)
        with torch.set_grad_enabled(training):
            with torch.autocast(device_type=device.type, enabled=device.type == "cuda"):
                logits = model(x, pad_mask)
                loss = masked_loss(logits, target)
            if training:
                optimizer.zero_grad(set_to_none=True)
                if scaler is None:
                    loss.backward()
                    nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    optimizer.step()
                else:
                    scaler.scale(loss).backward()
                    scaler.unscale_(optimizer)
                    nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    scaler.step(optimizer)
                    scaler.update()
        loss_sum += loss.item() * count
        residue_count += count
        _update_confusion(confusion, batch, logits.detach())
    if residue_count == 0:
        raise ValueError("Loader contains no valid residues")
    result = confusion_metrics(confusion)
    result["loss"] = round(loss_sum / residue_count, 6)
    return result


class PredictionWriter:
    """Write one [L, 8] probability array per uid into an indexed .npy file."""

    def __init__(self, path: str | Path, entries: list[dict], metadata: dict):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.tmp_path = self.path.with_name(self.path.name + ".tmp")
        self.entries = entries
        self.metadata = metadata
        self.offsets = np.concatenate(([0], np.cumsum([e["length"] for e in entries], dtype=np.int64)))
        self.rows = {e["uid"]: i for i, e in enumerate(entries)}
        if len(self.rows) != len(entries):
            raise ValueError("Duplicate uid in prediction manifest")
        self.seen = set()
        self.array = np.lib.format.open_memmap(self.tmp_path, mode="w+", dtype=np.float32, shape=(int(self.offsets[-1]), 8))

    def write_batch(self, batch: dict, probs: np.ndarray) -> None:
        for j, uid in enumerate(batch["uid"]):
            if uid not in self.rows or uid in self.seen:
                raise ValueError(f"Unexpected or repeated uid: {uid}")
            i = self.rows[uid]
            length = int(batch["lengths"][j])
            if length != self.entries[i]["length"]:
                raise ValueError(f"Length mismatch for {uid}")
            p = probs[j, :length]
            if not np.isfinite(p).all() or not np.allclose(p.sum(axis=-1), 1, atol=1e-4):
                raise ValueError(f"Invalid probabilities for {uid}")
            self.array[self.offsets[i]:self.offsets[i + 1]] = p
            self.seen.add(uid)

    def finish(self) -> None:
        if len(self.seen) != len(self.entries):
            missing = set(self.rows) - self.seen
            raise ValueError(f"Missing predictions for {len(missing)} proteins")
        self.array.flush()
        self.array._mmap.close()
        del self.array
        gc.collect()
        index = {
            "format": "indexed_protein_probabilities_v1",
            "probabilities": self.path.name,
            "shape": [int(self.offsets[-1]), 8],
            "q8_order": "GHIBESTC",
            "uids": [e["uid"] for e in self.entries],
            "offsets": self.offsets.tolist(),
            "fold_ids": [e["fold_id"] for e in self.entries],
            "metadata": self.metadata,
        }
        index_path = self.path.with_suffix(".index.json")
        tmp_index = index_path.with_name(index_path.name + ".tmp")
        tmp_index.write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
        os.replace(self.tmp_path, self.path)
        os.replace(tmp_index, index_path)


class PredictionStore:
    """Read exported probabilities by unique protein uid without loading them all."""

    def __init__(self, path: str | Path):
        path = Path(path)
        index = json.loads(path.with_suffix(".index.json").read_text(encoding="utf-8"))
        if index["format"] != "indexed_protein_probabilities_v1" or index["q8_order"] != "GHIBESTC":
            raise ValueError("Unsupported prediction index")
        self.uids = index["uids"]
        self.offsets = index["offsets"]
        self.fold_ids = index["fold_ids"]
        self.metadata = index["metadata"]
        self.row = {uid: i for i, uid in enumerate(self.uids)}
        self.array = np.load(path, mmap_mode="r")
        if (len(self.row) != len(self.uids) or len(self.offsets) != len(self.uids) + 1
                or self.array.shape != tuple(index["shape"])):
            raise ValueError("Prediction index does not match the array")

    def get(self, uid: str) -> np.ndarray:
        i = self.row[uid]
        return self.array[self.offsets[i]:self.offsets[i + 1]]


def dataset_entries(dataset) -> list[dict]:
    return [
        {"uid": record.uid, "length": record.length, "fold_id": fold_id}
        for record, fold_id in zip(dataset.records, dataset.fold_ids)
    ]


@torch.no_grad()
def export_predictions(model: BaselineModel, loader, writer: PredictionWriter, device: torch.device) -> dict:
    model.eval()
    confusion = torch.zeros((8, 8), dtype=torch.long)
    for batch in loader:
        x = batch["esm_emb"].to(device, non_blocking=True)
        pad_mask = batch["pad_mask"].to(device, non_blocking=True)
        with torch.autocast(device_type=device.type, enabled=device.type == "cuda"):
            logits = model(x, pad_mask)
        probs = logits.float().softmax(dim=-1).cpu().numpy()
        writer.write_batch(batch, probs)
        _update_confusion(confusion, batch, logits)
    return confusion_metrics(confusion)


def checkpoint_payload(model: BaselineModel, *, fold_id: int | str, seed: int, epoch: int,
                       validation: dict, split_sha256: str, esm_model: str) -> dict:
    return {
        "architecture": "frozen_esm2_1layer_bilstm_q8",
        "config": asdict(model.config),
        "model_state": {key: value.cpu().clone() for key, value in model.state_dict().items()},
        "fold_id": fold_id,
        "seed": seed,
        "epoch": epoch,
        "validation": validation,
        "split_sha256": split_sha256,
        "esm_model": esm_model,
    }


def save_checkpoint(payload: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    torch.save(payload, tmp)
    os.replace(tmp, path)
