"""Role 1, step 1-2: parse the NetSurfP-2.0 archives and audit masks/labels.

Outputs (under data/processed/):
  train.npz, cb513.npz     parsed records (see e2b/netsurfp.py: save_records)
  train.fasta, cb513.fasta sequences keyed by uid, for similarity search
  data_audit.json          mask / label / length checks

Usage:
  python scripts/prepare_data.py
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e2b.constants import ESM_MAX_RESIDUES, LABEL_MISSING, Q8_ORDER  # noqa: E402
from e2b.netsurfp import parse_netsurfp, save_records, sha256sum, write_fasta  # noqa: E402


def audit(records, source):
    lengths = np.array([r.length for r in records])
    q8 = np.concatenate([r.q8_label for r in records])
    dis = np.concatenate([r.disorder_mask for r in records])
    ev = np.concatenate([r.eval_mask for r in records])
    va = np.concatenate([r.valid_mask for r in records])
    seqs = "".join(r.seq for r in records)
    labelled = q8 != LABEL_MISSING
    counts = np.bincount(q8[va].astype(np.int64), minlength=8)
    pdb_counts = Counter(r.pdb_id for r in records)
    return {
        "source": source,
        "n_proteins": len(records),
        "n_unique_pdb_ids": len(pdb_counts),
        "length": {
            "min": int(lengths.min()),
            "median": float(np.median(lengths)),
            "mean": round(float(lengths.mean()), 1),
            "max": int(lengths.max()),
            f"n_longer_than_{ESM_MAX_RESIDUES}": int((lengths > ESM_MAX_RESIDUES).sum()),
        },
        "residues_total": int(lengths.sum()),
        "residues_unknown_X": seqs.count("X"),
        "residues_disordered": int((~dis).sum()),
        "residues_without_q8_label": int((~labelled).sum()),
        "residues_disordered_but_labelled": int((~dis & labelled).sum()),
        "residues_eval_mask_0": int((~ev).sum()),
        "residues_valid": int(va.sum()),
        "q8_counts_valid": {c: int(n) for c, n in zip(Q8_ORDER, counts)},
        "q8_percent_valid": {c: round(100 * float(n) / counts.sum(), 2) for c, n in zip(Q8_ORDER, counts)},
        "proteins_with_zero_valid_residues": int(sum(1 for r in records if not r.valid_mask.any())),
    }


def exact_overlap(train, cb513):
    train_seqs = {}
    for r in train:
        train_seqs.setdefault(r.seq, []).append(r.uid)
    dup_within_train = {s: u for s, u in train_seqs.items() if len(u) > 1}
    joined = "\n".join(r.seq for r in train)
    identical, cb_in_train = [], []
    for r in cb513:
        if r.seq in train_seqs:
            identical.append((r.uid, train_seqs[r.seq]))
        elif r.seq in joined:
            cb_in_train.append(r.uid)
    cb_joined = "\n".join(r.seq for r in cb513)
    cb_set = {r.seq for r in cb513}
    train_in_cb = [r.uid for r in train if r.seq not in cb_set and r.seq in cb_joined]
    return {
        "train_duplicate_sequences": len(dup_within_train),
        "cb513_identical_to_train": identical,
        "cb513_substring_of_train": cb_in_train,
        "train_substring_of_cb513": train_in_cb,
    }


def main():
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=root / "data/raw", type=Path)
    ap.add_argument("--out", default=root / "data/processed", type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    sources = {"train": args.raw / "Train_HHblits.npz", "cb513": args.raw / "CB513_HHblits.npz"}
    parsed, report = {}, {"files": {}}
    for name, path in sources.items():
        print(f"parsing {path.name} ...", flush=True)
        recs = parse_netsurfp(path, name)
        parsed[name] = recs
        save_records(recs, args.out / f"{name}.npz")
        write_fasta(recs, args.out / f"{name}.fasta")
        report["files"][path.name] = {"sha256": sha256sum(path), "bytes": path.stat().st_size}
        report[name] = audit(recs, name)
        print(json.dumps(report[name], indent=1), flush=True)

    report["exact_overlap"] = exact_overlap(parsed["train"], parsed["cb513"])
    print(json.dumps(report["exact_overlap"], indent=1))
    (args.out / "data_audit.json").write_text(json.dumps(report, indent=2))
    print(f"wrote {args.out / 'data_audit.json'}")


if __name__ == "__main__":
    main()
