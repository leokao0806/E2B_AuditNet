"""Shared PyTorch Dataset / DataLoader (Shared Data Contract) for Roles 3-6.

Per-protein item (dict):
  pdb_id     str
  uid        str            unique key (train_NNNNN / cb513_NNNNN)
  seq        str            length L
  esm_emb    FloatTensor    [L, D] frozen ESM-2 features
  q8_label   LongTensor     [L] 0..7 = G,H,I,B,E,S,T,C  (-1 where no label; never valid)
  valid_mask BoolTensor     [L] ordered & labelled (& CB513 eval mask)
  fold_id    int            0..K-1 CV, -1 validation, -2 CB513

Batch (collate_proteins), padded to the longest protein Lmax:
  esm_emb [B, Lmax, D], q8_label [B, Lmax] (pad -1), valid_mask [B, Lmax] (pad False),
  target [B, Lmax] = q8_label where valid_mask else -100 (ready for cross_entropy),
  pad_mask [B, Lmax] True on real residues, lengths [B], fold_id [B], plus lists pdb_id/uid/seq.

Typical use:
  from e2b.dataset import make_loader
  train_dl = make_loader(folds=[1, 2, 3, 4], batch_size=16, shuffle=True)
  oof_dl   = make_loader(folds=[0], batch_size=32)
  test_dl  = make_loader(folds=[-2], batch_size=32)
  loss = F.cross_entropy(logits.transpose(1, 2), batch["target"], ignore_index=-100)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset, Sampler

from .constants import IGNORE_INDEX, LABEL_MISSING
from .esm_embed import MODELS, cache_paths
from .netsurfp import load_records
from .splits import fold_lookup, load_splits

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = "esm2_t33_650M_UR50D"


class EmbeddingStore:
    """Read-only view of one cached source. Opens the memmap lazily (safe for worker processes)."""

    def __init__(self, emb_dir: str | Path, source: str, model: str = DEFAULT_MODEL, layer: int | None = None):
        layer = MODELS[model][0] if layer is None else layer
        self.npy_path, idx_path = cache_paths(emb_dir, source, model, layer)
        if not self.npy_path.exists():
            raise FileNotFoundError(f"{self.npy_path} not found: run scripts/extract_esm.py --model {model}")
        meta = json.loads(idx_path.read_text())
        self.uids = meta["uids"]
        self.offsets = np.asarray(meta["offsets"], dtype=np.int64)
        self.dim = int(meta["dim"])
        self.row = {u: i for i, u in enumerate(self.uids)}
        self._mm = None

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_mm"] = None  # never pickle the mapped array
        return state

    def get(self, uid: str) -> np.ndarray:
        if self._mm is None:
            self._mm = np.load(self.npy_path, mmap_mode="r")
        i = self.row[uid]
        return self._mm[self.offsets[i] : self.offsets[i + 1]]


class ProteinDataset(Dataset):
    def __init__(
        self,
        folds,
        splits_path: str | Path = ROOT / "splits/splits.json",
        processed_dir: str | Path = ROOT / "data/processed",
        emb_dir: str | Path = ROOT / "data/embeddings",
        model: str = DEFAULT_MODEL,
        layer: int | None = None,
    ):
        self.folds = sorted(set(int(f) for f in folds))
        fold_of = fold_lookup(load_splits(splits_path))
        sources = {"cb513"} if self.folds == [-2] else ({"train"} if -2 not in self.folds else {"train", "cb513"})
        self.records, self.fold_ids, self.stores = [], [], {}
        for src in sorted(sources):
            store = EmbeddingStore(emb_dir, src, model, layer)
            self.stores[src] = store
            for r in load_records(Path(processed_dir) / f"{src}.npz"):
                f = fold_of.get(r.uid)
                if f is None or f not in self.folds:
                    continue
                if r.uid not in store.row:
                    raise KeyError(f"{r.uid} missing from embedding cache {store.npy_path.name}")
                i = store.row[r.uid]
                if store.offsets[i + 1] - store.offsets[i] != r.length:
                    raise ValueError(f"{r.uid}: cache length does not match sequence length")
                self.records.append(r)
                self.fold_ids.append(f)
        self.dim = next(iter(self.stores.values())).dim

    def __len__(self):
        return len(self.records)

    @property
    def lengths(self) -> list[int]:
        return [r.length for r in self.records]

    def __getitem__(self, i):
        r = self.records[i]
        emb = self.stores[r.source].get(r.uid)
        return {
            "pdb_id": r.pdb_id,
            "uid": r.uid,
            "seq": r.seq,
            "esm_emb": torch.from_numpy(np.asarray(emb, dtype=np.float32)),
            "q8_label": torch.from_numpy(r.q8_label.astype(np.int64)),
            "valid_mask": torch.from_numpy(r.valid_mask.copy()),
            "fold_id": self.fold_ids[i],
        }


def collate_proteins(batch: list[dict]) -> dict:
    lengths = torch.tensor([b["q8_label"].shape[0] for b in batch])
    B, Lmax, D = len(batch), int(lengths.max()), batch[0]["esm_emb"].shape[1]
    emb = torch.zeros(B, Lmax, D)
    q8 = torch.full((B, Lmax), LABEL_MISSING, dtype=torch.long)
    valid = torch.zeros(B, Lmax, dtype=torch.bool)
    pad = torch.zeros(B, Lmax, dtype=torch.bool)
    for i, b in enumerate(batch):
        L = lengths[i]
        emb[i, :L] = b["esm_emb"]
        q8[i, :L] = b["q8_label"]
        valid[i, :L] = b["valid_mask"]
        pad[i, :L] = True
    target = q8.masked_fill(~valid, IGNORE_INDEX)
    return {
        "esm_emb": emb,
        "q8_label": q8,
        "valid_mask": valid,
        "target": target,
        "pad_mask": pad,
        "lengths": lengths,
        "fold_id": torch.tensor([b["fold_id"] for b in batch]),
        "pdb_id": [b["pdb_id"] for b in batch],
        "uid": [b["uid"] for b in batch],
        "seq": [b["seq"] for b in batch],
    }


class LengthBucketSampler(Sampler):
    """Batches of similar-length proteins (less padding). Reshuffled every epoch when shuffle=True."""

    def __init__(self, lengths, batch_size: int, shuffle: bool = True, seed: int = 0, bucket_factor: int = 50):
        self.lengths = np.asarray(lengths)
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.seed = seed
        self.bucket = batch_size * bucket_factor
        self.epoch = 0

    def __iter__(self):
        rng = np.random.default_rng(self.seed + self.epoch)
        self.epoch += 1
        idx = rng.permutation(len(self.lengths)) if self.shuffle else np.arange(len(self.lengths))
        batches = []
        for s in range(0, len(idx), self.bucket):
            chunk = idx[s : s + self.bucket]
            chunk = chunk[np.argsort(self.lengths[chunk], kind="stable")]
            batches += [chunk[j : j + self.batch_size].tolist() for j in range(0, len(chunk), self.batch_size)]
        if self.shuffle:
            rng.shuffle(batches)
        return iter(batches)

    def __len__(self):
        n = 0
        for s in range(0, len(self.lengths), self.bucket):
            m = min(self.bucket, len(self.lengths) - s)
            n += -(-m // self.batch_size)
        return n


def make_loader(
    folds,
    batch_size: int = 16,
    shuffle: bool = False,
    num_workers: int = 0,
    bucket_by_length: bool = True,
    seed: int = 0,
    dataset: Dataset | None = None,
    **dataset_kwargs,
) -> DataLoader:
    ds = dataset if dataset is not None else ProteinDataset(folds, **dataset_kwargs)
    common = dict(collate_fn=collate_proteins, num_workers=num_workers, persistent_workers=num_workers > 0)
    if bucket_by_length:
        sampler = LengthBucketSampler(ds.lengths, batch_size, shuffle=shuffle, seed=seed)
        return DataLoader(ds, batch_sampler=sampler, **common)
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, **common)


def make_oof_loaders(k: int, batch_size: int = 16, **kw) -> tuple[DataLoader, DataLoader]:
    """(train loader on folds != k, shuffled; held-out loader on fold k, in order)."""
    splits = load_splits(kw.get("splits_path", ROOT / "splits/splits.json"))
    folds = range(splits["params"]["k_folds"])
    train = make_loader([f for f in folds if f != k], batch_size, shuffle=True, **kw)
    held = make_loader([k], batch_size, shuffle=False, **kw)
    return train, held


# ---------------------------------------------------------------------------
# Mock data with the same contract, for developing models before the cache exists
# ---------------------------------------------------------------------------


class MockProteinDataset(Dataset):
    """Synthetic proteins following the Shared Data Contract (random features and labels)."""

    def __init__(self, n: int = 64, dim: int = 1280, min_len: int = 20, max_len: int = 400, fold_id: int = 0, seed: int = 0):
        rng = np.random.default_rng(seed)
        q8_freq = np.array([3.9, 34.3, 0.6, 1.0, 21.7, 8.1, 10.9, 19.5])
        self.items = []
        for i in range(n):
            L = int(rng.integers(min_len, max_len + 1))
            q8 = rng.choice(8, size=L, p=q8_freq / q8_freq.sum())
            valid = rng.random(L) > 0.05
            self.items.append(
                {
                    "pdb_id": f"mock-{i}",
                    "uid": f"mock_{i:05d}",
                    "seq": "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=L)),
                    "esm_emb": torch.from_numpy(rng.standard_normal((L, dim)).astype(np.float32)),
                    "q8_label": torch.from_numpy(q8.astype(np.int64)),
                    "valid_mask": torch.from_numpy(valid),
                    "fold_id": fold_id,
                }
            )
        self.dim = dim

    @property
    def lengths(self):
        return [len(x["seq"]) for x in self.items]

    def __len__(self):
        return len(self.items)

    def __getitem__(self, i):
        return self.items[i]
