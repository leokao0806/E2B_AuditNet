"""Frozen ESM-2 per-residue feature extraction with long-sequence windowing.

ESM-2 was trained on crops of at most 1,024 tokens (BOS + 1,022 residues + EOS).
Sequences longer than ``max_residues`` are split into overlapping windows; each
residue takes its embedding from the window in which it sits farthest from a
window edge (i.e. with the most context on both sides). Every residue therefore
gets exactly one embedding that the model actually produced, with no averaging.

Embeddings are written to a single float16 memory-mapped .npy per source:
  <out>/<source>.<model>.L<layer>.f16.npy   shape [total_residues, D]
  <out>/<source>.<model>.L<layer>.index.json
The row range of protein i is offsets[i]:offsets[i+1] in the matching
data/processed/<source>.npz, so no per-protein files are needed.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from .constants import ESM_MAX_RESIDUES

MODELS = {
    # name: (layers, embedding dim)
    "esm2_t6_8M_UR50D": (6, 320),
    "esm2_t12_35M_UR50D": (12, 480),
    "esm2_t30_150M_UR50D": (30, 640),
    "esm2_t33_650M_UR50D": (33, 1280),
}


def plan_windows(length: int, max_len: int = ESM_MAX_RESIDUES, overlap: int = 256) -> list[tuple[int, int]]:
    """Window [start, end) ranges covering 0..length with the given overlap."""
    if length <= max_len:
        return [(0, length)]
    if not 0 <= overlap < max_len:
        raise ValueError("overlap must be in [0, max_len)")
    step = max_len - overlap
    starts = list(range(0, length - max_len, step)) + [length - max_len]
    return [(s, s + max_len) for s in sorted(set(starts))]


def assign_positions(length: int, windows: list[tuple[int, int]]) -> np.ndarray:
    """For each residue, index of the window where it has the most two-sided context."""
    best = np.full(length, -1, dtype=np.int64)
    best_ctx = np.full(length, -1, dtype=np.int64)
    pos = np.arange(length)
    for w, (s, e) in enumerate(windows):
        inside = (pos >= s) & (pos < e)
        ctx = np.minimum(pos - s, e - 1 - pos)
        better = inside & (ctx > best_ctx)
        best[better] = w
        best_ctx[better] = ctx[better]
    if (best < 0).any():
        raise AssertionError("windows do not cover the sequence")
    return best


@dataclass
class EsmEmbedder:
    model_name: str = "esm2_t33_650M_UR50D"
    layer: int | None = None  # default: last layer
    device: str = "cuda"
    half: bool = True
    max_residues: int = ESM_MAX_RESIDUES
    overlap: int = 256
    # Cap PyTorch's share of VRAM. On Windows (WDDM) allocations past the card's
    # memory silently spill into shared system RAM and run several times slower;
    # with a cap the caching allocator frees blocks instead of spilling.
    gpu_memory_fraction: float | None = 0.75

    def __post_init__(self):
        import esm  # fair-esm

        if self.device.startswith("cuda") and self.gpu_memory_fraction:
            torch.cuda.set_per_process_memory_fraction(self.gpu_memory_fraction)

        if self.model_name not in MODELS:
            raise ValueError(f"unknown model {self.model_name}; choose from {list(MODELS)}")
        n_layers, self.dim = MODELS[self.model_name]
        self.layer = n_layers if self.layer is None else self.layer
        model, alphabet = getattr(esm.pretrained, self.model_name)()
        model.eval().requires_grad_(False)  # frozen
        if self.half and self.device.startswith("cuda"):
            model = model.half()
        self.model = model.to(self.device)
        self.alphabet = alphabet
        self.batch_converter = alphabet.get_batch_converter()

    @torch.inference_mode()
    def _forward(self, seqs: list[str]) -> list[torch.Tensor]:
        """Embed sequences that each fit in one window. Returns float32 CPU [L, D] tensors."""
        _, _, tokens = self.batch_converter([(str(i), s) for i, s in enumerate(seqs)])
        out = self.model(tokens.to(self.device), repr_layers=[self.layer], return_contacts=False)
        rep = out["representations"][self.layer]
        return [rep[i, 1 : len(s) + 1].float().cpu() for i, s in enumerate(seqs)]

    def embed_long(self, seq: str) -> torch.Tensor:
        windows = plan_windows(len(seq), self.max_residues, self.overlap)
        which = assign_positions(len(seq), windows)
        out = torch.empty(len(seq), self.dim)
        for w, (s, e) in enumerate(windows):
            rep = self._forward_safe([seq[s:e]])[0]
            sel = np.nonzero(which == w)[0]
            out[sel] = rep[sel - s]
        return out

    def _forward_safe(self, seqs: list[str]) -> list[torch.Tensor]:
        """_forward, splitting the batch in half on CUDA out-of-memory."""
        try:
            return self._forward(seqs)
        except torch.OutOfMemoryError:
            if len(seqs) == 1:
                raise
            torch.cuda.empty_cache()
            half = len(seqs) // 2
            return self._forward_safe(seqs[:half]) + self._forward_safe(seqs[half:])

    @staticmethod
    def make_batches(lengths, tokens_per_batch: int = 8192, attn_per_batch: int = 4 * 1024**2) -> list[list[int]]:
        """Group indices (sorted by length) so that B*(L+2) <= tokens and B*(L+2)^2 <= attn budget.

        The attention budget matters for long proteins: fair-esm materialises
        [B, heads, L, L] attention weights (softmax in float32) in every layer.
        """
        order = sorted(range(len(lengths)), key=lambda i: lengths[i])
        batches, cur = [], []
        for i in order:
            n = len(cur) + 1
            t = lengths[i] + 2  # i is the longest so far
            if cur and (n * t > tokens_per_batch or n * t * t > attn_per_batch):
                batches.append(cur)
                cur = []
            cur.append(i)
        if cur:
            batches.append(cur)
        return batches

    def embed_many(self, seqs: list[str], tokens_per_batch: int = 8192, progress: bool = True):
        """Yield (index, [L, D] float32 tensor) for every sequence, batching by length."""
        short = [i for i in range(len(seqs)) if len(seqs[i]) <= self.max_residues]
        long_ = [i for i in range(len(seqs)) if len(seqs[i]) > self.max_residues]
        batches = [[short[j] for j in b] for b in self.make_batches([len(seqs[i]) for i in short], tokens_per_batch)]
        t0, done, res_done, next_report = time.time(), 0, 0, 500
        total_res = sum(len(s) for s in seqs)
        for batch in batches:
            for j, rep in zip(batch, self._forward_safe([seqs[j] for j in batch])):
                yield j, rep
            done += len(batch)
            res_done += sum(len(seqs[j]) for j in batch)
            if progress and done >= next_report:
                next_report += 500
                el = time.time() - t0
                print(f"  {done}/{len(seqs)} proteins, {res_done / el:,.0f} residues/s, {el:.0f}s", flush=True)
        for i in long_:
            yield i, self.embed_long(seqs[i])
            res_done += len(seqs[i])
        if progress:
            el = time.time() - t0
            print(f"  done: {len(seqs)} proteins, {total_res:,} residues in {el:.1f}s", flush=True)


def cache_paths(out_dir: str | Path, source: str, model_name: str, layer: int) -> tuple[Path, Path]:
    stem = f"{source}.{model_name}.L{layer}"
    out_dir = Path(out_dir)
    return out_dir / f"{stem}.f16.npy", out_dir / f"{stem}.index.json"


def write_cache(embedder: EsmEmbedder, uids, seqs, offsets, out_dir, source, tokens_per_batch=8192):
    """Write all embeddings for one source into a float16 memmap aligned with `offsets`."""
    npy_path, idx_path = cache_paths(out_dir, source, embedder.model_name, embedder.layer)
    npy_path.parent.mkdir(parents=True, exist_ok=True)
    total = int(offsets[-1])
    tmp_path = npy_path.with_suffix(".partial.npy")
    mm = np.lib.format.open_memmap(tmp_path, mode="w+", dtype=np.float16, shape=(total, embedder.dim))
    written = np.zeros(len(seqs), dtype=bool)
    t0 = time.time()
    for i, rep in embedder.embed_many(seqs, tokens_per_batch=tokens_per_batch):
        a, b = int(offsets[i]), int(offsets[i + 1])
        if rep.shape != (b - a, embedder.dim):
            raise AssertionError(f"{uids[i]}: got {tuple(rep.shape)}, expected {(b - a, embedder.dim)}")
        mm[a:b] = rep.numpy().astype(np.float16)
        written[i] = True
    if not written.all():
        raise AssertionError(f"{(~written).sum()} proteins were not embedded")
    mm.flush()
    del mm
    tmp_path.replace(npy_path)
    meta = {
        "source": source,
        "model": embedder.model_name,
        "layer": embedder.layer,
        "dim": embedder.dim,
        "dtype": "float16",
        "half_precision_forward": embedder.half,
        "max_residues_per_window": embedder.max_residues,
        "window_overlap": embedder.overlap,
        "n_proteins": len(seqs),
        "total_residues": total,
        "seconds": round(time.time() - t0, 1),
        "uids": list(map(str, uids)),
        "offsets": [int(x) for x in offsets],
    }
    idx_path.write_text(json.dumps(meta))
    return npy_path, idx_path
