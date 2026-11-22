"""Role 2 hand-off check: is the embedding cache aligned with the labels, and how fast is the loader?

A per-residue linear probe (one nn.Linear on frozen ESM-2 features) is trained on
the CV folds and scored on validation (-1) and CB513 (-2) with valid_mask. If
features and labels were misaligned the probe would sit near the majority-class
rate (~34% Q8); a correct cache gives a clearly higher score. This is only a
sanity check, not the project baseline (Role 3 trains the BiLSTM).

Usage:
  python scripts/check_cache.py [--model esm2_t33_650M_UR50D] [--epochs 2]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e2b.constants import G_IDX, H_IDX, IGNORE_INDEX, Q8_TO_Q3  # noqa: E402
from e2b.dataset import make_loader  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


@torch.no_grad()
def evaluate(probe, loader, device):
    conf = torch.zeros(8, 8, dtype=torch.long)
    for b in loader:
        logits = probe(b["esm_emb"].to(device))
        pred = logits.argmax(-1).cpu()
        m = b["valid_mask"]
        conf += torch.bincount(b["q8_label"][m] * 8 + pred[m], minlength=64).reshape(8, 8)
    n = conf.sum().item()
    q8 = conf.trace().item() / n
    q3_map = torch.tensor(Q8_TO_Q3)
    c3 = torch.zeros(3, 3, dtype=torch.long)
    for t in range(8):
        for p in range(8):
            c3[q3_map[t], q3_map[p]] += conf[t, p]
    g_recall = conf[G_IDX, G_IDX].item() / max(conf[G_IDX].sum().item(), 1)
    g_as_h = conf[G_IDX, H_IDX].item() / max(conf[G_IDX].sum().item(), 1)
    return {"residues": n, "Q8": round(100 * q8, 2), "Q3": round(100 * c3.trace().item() / n, 2),
            "G_recall": round(100 * g_recall, 2), "G_predicted_as_H": round(100 * g_as_h, 2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--epochs", type=int, default=2)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--num-workers", type=int, default=0)
    args = ap.parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(0)

    t0 = time.time()
    train_dl = make_loader([0, 1, 2, 3, 4], args.batch_size, shuffle=True, num_workers=args.num_workers, model=args.model)
    val_dl = make_loader([-1], 64, model=args.model)
    test_dl = make_loader([-2], 64, model=args.model)
    print(f"datasets: train {len(train_dl.dataset)}, val {len(val_dl.dataset)}, CB513 {len(test_dl.dataset)} "
          f"(built in {time.time() - t0:.1f}s), dim {train_dl.dataset.dim}", flush=True)

    probe = torch.nn.Linear(train_dl.dataset.dim, 8).to(device)
    opt = torch.optim.AdamW(probe.parameters(), lr=1e-3, weight_decay=1e-4)
    for ep in range(args.epochs):
        t0, n_res = time.time(), 0
        for b in train_dl:
            logits = probe(b["esm_emb"].to(device, non_blocking=True))
            loss = F.cross_entropy(logits.transpose(1, 2), b["target"].to(device), ignore_index=IGNORE_INDEX)
            opt.zero_grad()
            loss.backward()
            opt.step()
            n_res += int(b["lengths"].sum())
        el = time.time() - t0
        print(f"epoch {ep + 1}: loss {loss.item():.3f}, {len(train_dl)} batches in {el:.1f}s "
              f"({len(train_dl) / el:.1f} batches/s, {n_res / el:,.0f} residues/s)", flush=True)

    res = {"model": args.model, "probe": "linear", "epochs": args.epochs,
           "validation": evaluate(probe, val_dl, device), "cb513": evaluate(probe, test_dl, device)}
    print(json.dumps(res, indent=1))
    with open(ROOT / "data/embeddings/check_cache.jsonl", "a") as f:
        f.write(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
