"""Role 1, step 3: leakage removal, homology-grouped folds, and split statistics.

Inputs:
  data/processed/{train,cb513}.npz and data_audit.json   (scripts/prepare_data.py)
  data/similarity/cb513_vs_train.m8, train_vs_train.m8    (scripts/run_mmseqs.py; optional)

Outputs:
  splits/splits.json        fold_id for every usable protein + excluded list
  splits/split_stats.json   per-partition residue / Q8 statistics
  splits/split_stats.md     the same as Markdown tables (for the report / poster)

Usage:
  python scripts/make_splits.py [--k 5] [--val-fraction 0.1] [--seed 2026]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e2b.constants import ESM_MAX_RESIDUES, FOLD_CB513, FOLD_VALIDATION, G_IDX, H_IDX, Q8_ORDER  # noqa: E402
from e2b.netsurfp import load_records  # noqa: E402
from e2b.splits import FOLD_MEANING, HomologyRule, assign_folds, connected_groups, homology_edges, read_m8  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def g_next_to_h(q8: np.ndarray) -> np.ndarray:
    """True for G residues whose G segment directly touches an H residue."""
    out = np.zeros(len(q8), dtype=bool)
    i = 0
    while i < len(q8):
        if q8[i] != G_IDX:
            i += 1
            continue
        j = i
        while j < len(q8) and q8[j] == G_IDX:
            j += 1
        if (i > 0 and q8[i - 1] == H_IDX) or (j < len(q8) and q8[j] == H_IDX):
            out[i:j] = True
        i = j
    return out


def partition_stats(recs) -> dict:
    if not recs:
        return {"proteins": 0}
    lengths = np.array([r.length for r in recs])
    q8 = np.concatenate([r.q8_label for r in recs]).astype(np.int64)
    va = np.concatenate([r.valid_mask for r in recs])
    gh = np.concatenate([g_next_to_h(r.q8_label) for r in recs])
    counts = np.bincount(q8[va], minlength=8)
    n_valid = int(va.sum())
    n_g = int(counts[G_IDX])
    return {
        "proteins": len(recs),
        "residues": int(lengths.sum()),
        "valid_residues": n_valid,
        "median_length": float(np.median(lengths)),
        "max_length": int(lengths.max()),
        f"proteins_longer_than_{ESM_MAX_RESIDUES}": int((lengths > ESM_MAX_RESIDUES).sum()),
        "q8_counts": {c: int(n) for c, n in zip(Q8_ORDER, counts)},
        "q8_percent": {c: round(100 * n / max(n_valid, 1), 2) for c, n in zip(Q8_ORDER, counts)},
        "G_next_to_H": int((gh & va).sum()),
        "G_next_to_H_percent_of_G": round(100 * int((gh & va).sum()) / max(n_g, 1), 1),
    }


def stats_markdown(stats: dict, meta: dict) -> str:
    names = list(stats)
    lines = [
        "# Split statistics (Role 1)",
        "",
        f"Generated {meta['created']} by `scripts/make_splits.py`. Residue counts use `valid_mask` "
        "(ordered residues with a Q8 label; on CB513 also the evaluation mask).",
        "",
        f"- CB513 leakage rule (Train proteins removed): {meta['homology']['cb513_leakage_rule']}.",
        f"- Train grouping rule (same fold): {meta['homology']['train_grouping_rule']}.",
        f"- Status: **{meta['homology']['status']}**; {meta['homology']['train_groups']:,} groups, "
        f"largest {meta['homology']['largest_group']} proteins.",
        "",
        "| Partition | Proteins | Valid residues | Median len | Max len | >1022 | G | G % | H | H % | G next to H (% of G) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for n in names:
        s = stats[n]
        if not s["proteins"]:
            lines.append(f"| {n} | 0 | | | | | | | | | |")
            continue
        lines.append(
            f"| {n} | {s['proteins']:,} | {s['valid_residues']:,} | {s['median_length']:.0f} | {s['max_length']} | "
            f"{s[f'proteins_longer_than_{ESM_MAX_RESIDUES}']} | {s['q8_counts']['G']:,} | {s['q8_percent']['G']:.2f} | "
            f"{s['q8_counts']['H']:,} | {s['q8_percent']['H']:.2f} | {s['G_next_to_H']:,} ({s['G_next_to_H_percent_of_G']:.1f}%) |"
        )
    lines += ["", "## Q8 distribution (% of valid residues)", "",
              "| Partition | " + " | ".join(Q8_ORDER) + " |", "|---|" + "---:|" * 8]
    for n in names:
        s = stats[n]
        if s["proteins"]:
            lines.append(f"| {n} | " + " | ".join(f"{s['q8_percent'][c]:.2f}" for c in Q8_ORDER) + " |")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed", type=Path, default=ROOT / "data/processed")
    ap.add_argument("--similarity", type=Path, default=ROOT / "data/similarity")
    ap.add_argument("--out", type=Path, default=ROOT / "splits")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--val-fraction", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=2026)
    # CB513 leakage: strict, any significant local match counts (only ~3% of Train is lost).
    ap.add_argument("--test-max-evalue", type=float, default=1e-3)
    ap.add_argument("--test-min-identity", type=float, default=0.25)
    ap.add_argument("--test-min-alnlen", type=int, default=20)
    # Grouping inside Train: needs coverage, otherwise expression tags chain everything together.
    ap.add_argument("--group-max-evalue", type=float, default=1e-3)
    ap.add_argument("--group-min-identity", type=float, default=0.30)
    ap.add_argument("--group-min-coverage", type=float, default=0.5)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    train = load_records(args.processed / "train.npz")
    cb513 = load_records(args.processed / "cb513.npz")
    audit = json.loads((args.processed / "data_audit.json").read_text())
    by_uid = {r.uid: r for r in train + cb513}
    test_rule = HomologyRule(args.test_max_evalue, args.test_min_identity, args.test_min_alnlen)
    group_rule = HomologyRule(args.group_max_evalue, args.group_min_identity, 20, args.group_min_coverage)

    # 1) Train proteins that overlap CB513 are removed from training entirely.
    excluded: dict[str, str] = {}
    ov = audit["exact_overlap"]
    for cb_uid, train_uids in ov["cb513_identical_to_train"]:
        for t in train_uids:
            excluded[t] = f"identical sequence to {cb_uid} ({by_uid[cb_uid].pdb_id})"
    joined_cb = {r.uid: r.seq for r in cb513}
    for cb_uid in ov["cb513_substring_of_train"]:
        for t in train:
            if joined_cb[cb_uid] in t.seq:
                excluded.setdefault(t.uid, f"contains CB513 {cb_uid} ({by_uid[cb_uid].pdb_id}) as a substring")
    for t in ov["train_substring_of_cb513"]:
        excluded.setdefault(t, "substring of a CB513 sequence")

    cross_file = args.similarity / "cb513_vs_train.m8"
    within_file = args.similarity / "train_vs_train.m8"
    homology_status = "exact-match and substring checks only; MMseqs2 search not run yet"
    n_cross_hits = n_within_edges = None
    leakage_table: list[dict] = []
    if cross_file.exists():
        hits = [h for h in read_m8(cross_file) if h["query"].startswith("cb513_") and h["target"].startswith("train_")]
        cross = [h for h in hits if test_rule.passes(h)]
        n_cross_hits = len(cross)
        for fid in (0.25, 0.30, 0.40, 0.50, 0.90):
            s = [h for h in hits if HomologyRule(test_rule.max_evalue, fid, test_rule.min_alnlen).passes(h)]
            leakage_table.append({"min_identity": fid, "pairs": len(s), "cb513_chains": len({h["query"] for h in s}),
                                  "train_chains": len({h["target"] for h in s})})
        for h in cross:
            excluded.setdefault(
                h["target"],
                f"homolog of CB513 {h['query']} ({by_uid[h['query']].pdb_id}): "
                f"identity {h['fident']:.2f}, alignment {h['alnlen']}, E={h['evalue']:.1e}",
            )
    # 2) Homologous training proteins are grouped so they share a fold.
    edges = set()
    if within_file.exists():
        edges = homology_edges(read_m8(within_file), group_rule)
        n_within_edges = len(edges)
    if cross_file.exists() and within_file.exists():
        homology_status = "complete (exact/substring + MMseqs2 CB513-vs-Train and Train-vs-Train)"

    usable = [r for r in train if r.uid not in excluded]
    groups = connected_groups([r.uid for r in usable], edges)
    folds = assign_folds(groups, {r.uid: int(r.valid_mask.sum()) for r in usable}, args.k, args.val_fraction, args.seed)

    group_sizes = np.bincount(list(groups.values()))
    proteins = [
        {"uid": r.uid, "pdb_id": r.pdb_id, "source": "train", "fold_id": folds[r.uid], "group": groups[r.uid]}
        for r in usable
    ] + [{"uid": r.uid, "pdb_id": r.pdb_id, "source": "cb513", "fold_id": FOLD_CB513, "group": -1} for r in cb513]

    meta = {
        "version": 1,
        "created": dt.date.today().isoformat(),
        "fold_id_meaning": FOLD_MEANING,
        "params": {"k_folds": args.k, "val_fraction": args.val_fraction, "seed": args.seed},
        "homology": {
            "cb513_leakage_rule": test_rule.describe(),
            "train_grouping_rule": group_rule.describe(),
            "status": homology_status,
            "cb513_vs_train_hits_passing_rule": n_cross_hits,
            "train_vs_train_edges": n_within_edges,
            "train_groups": int(len(group_sizes)),
            "largest_group": int(group_sizes.max()),
            "groups_with_more_than_one_protein": int((group_sizes > 1).sum()),
            "cb513_vs_train_by_identity": leakage_table,
        },
        "source_files": audit["files"],
        "notes": [
            "uid is the unique key (train_NNNNN / cb513_NNNNN = index in the original .npz). "
            "pdb_id is not unique in CB513 (513 regions from 434 chains).",
            "CB513 entries cb513_00299 (1rsy-A) and cb513_00377 (2alp-A) have an all-zero evaluation mask, "
            "so they contribute no scored residues.",
        ],
    }
    splits = dict(meta, excluded=[{"uid": u, "pdb_id": by_uid[u].pdb_id, "reason": why} for u, why in sorted(excluded.items())],
                  proteins=proteins)
    (args.out / "splits.json").write_text(json.dumps(splits, indent=1))

    parts = {f"fold {k}": [by_uid[p["uid"]] for p in proteins if p["fold_id"] == k] for k in range(args.k)}
    parts["all CV folds"] = [r for k in range(args.k) for r in parts[f"fold {k}"]]
    parts["validation (-1)"] = [by_uid[p["uid"]] for p in proteins if p["fold_id"] == FOLD_VALIDATION]
    parts["CB513 test (-2)"] = cb513
    parts["excluded"] = [by_uid[u] for u in excluded]
    stats = {name: partition_stats(recs) for name, recs in parts.items()}
    (args.out / "split_stats.json").write_text(json.dumps(dict(meta=meta, partitions=stats), indent=1))
    md = stats_markdown(stats, meta)
    if leakage_table:
        md += (
            "\n## CB513 vs Train similarity (MMseqs2, E-value <= 1e-3, alignment >= 20)\n\n"
            "NetSurfP-2.0 states that no training protein has > 25% identity to CB513. "
            "MMseqs2 local alignments still find the matches below, so these Train proteins are excluded.\n\n"
            "| Identity >= | Pairs | CB513 chains | Train chains |\n|---:|---:|---:|---:|\n"
        )
        md += "".join(f"| {r['min_identity']:.0%} | {r['pairs']} | {r['cb513_chains']} | {r['train_chains']} |\n"
                      for r in leakage_table)
    if excluded:
        md += "\n## Excluded training proteins\n\n| uid | pdb_id | reason |\n|---|---|---|\n"
        md += "".join(f"| {u} | {by_uid[u].pdb_id} | {why} |\n" for u, why in sorted(excluded.items()))
    (args.out / "split_stats.md").write_text(md, encoding="utf-8")
    print(md)
    print(json.dumps(meta["homology"], indent=1))


if __name__ == "__main__":
    main()
