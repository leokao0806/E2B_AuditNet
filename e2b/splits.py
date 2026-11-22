"""Similarity-safe fold assignment and the splits.json format.

Proteins linked by a homology edge (directly or through a chain of edges) form
one group, and a whole group always goes to the same partition. Validation
groups are drawn at random; the remaining groups are spread over K folds with a
greedy "largest group first, emptiest fold" rule so folds have similar residue
counts.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .constants import FOLD_CB513, FOLD_VALIDATION

FOLD_MEANING = {
    "0..K-1": "cross-validation / out-of-fold partitions (Train_HHblits)",
    "-1": "validation partition (Train_HHblits): class-weight study, routing thresholds, early stopping",
    "-2": "CB513 held-out test set: evaluate only, never tune on it",
}


# ---------------------------------------------------------------------------
# Homology hits (MMseqs2 m8 / tab format)
# ---------------------------------------------------------------------------

M8_COLUMNS = ["query", "target", "fident", "alnlen", "qlen", "tlen", "evalue", "bits"]


@dataclass(frozen=True)
class HomologyRule:
    """A hit links two proteins when all conditions hold.

    min_cov_shorter: alignment length / length of the shorter protein. Without it,
    the shared ~20-residue His-tag/thrombin expression tags (MGSSHHHHHHSSGLVPRGSH)
    link thousands of unrelated PDB chains.
    """

    max_evalue: float = 1e-3
    min_fident: float = 0.25
    min_alnlen: int = 20
    min_cov_shorter: float = 0.0

    def passes(self, h: dict) -> bool:
        return (
            h["evalue"] <= self.max_evalue
            and h["fident"] >= self.min_fident
            and h["alnlen"] >= self.min_alnlen
            and h["alnlen"] >= self.min_cov_shorter * min(h["qlen"], h["tlen"])
        )

    def describe(self) -> str:
        s = (
            f"E-value <= {self.max_evalue:g}, identity >= {self.min_fident:.0%}, "
            f"alignment >= {self.min_alnlen} residues"
        )
        if self.min_cov_shorter:
            s += f", covering >= {self.min_cov_shorter:.0%} of the shorter protein"
        return s


def read_m8(path: str | Path) -> list[dict]:
    hits = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(M8_COLUMNS):
                continue
            q, t, fid, aln, ql, tl, ev, bits = parts[: len(M8_COLUMNS)]
            hits.append(
                dict(query=q, target=t, fident=float(fid), alnlen=int(aln), qlen=int(ql), tlen=int(tl),
                     evalue=float(ev), bits=float(bits))
            )
    return hits


def homology_edges(hits: list[dict], rule: HomologyRule) -> set[tuple[str, str]]:
    edges = set()
    for h in hits:
        if h["query"] != h["target"] and rule.passes(h):
            edges.add(tuple(sorted((h["query"], h["target"]))))
    return edges


# ---------------------------------------------------------------------------
# Grouping and fold assignment
# ---------------------------------------------------------------------------


def connected_groups(uids: list[str], edges) -> dict[str, int]:
    """Union-find over `edges`; returns uid -> group id (0..G-1, ordered by first uid)."""
    parent = {u: u for u in uids}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
    roots, out = {}, {}
    for u in uids:
        out[u] = roots.setdefault(find(u), len(roots))
    return out


def assign_folds(
    groups: dict[str, int],
    residues: dict[str, int],
    k: int = 5,
    val_fraction: float = 0.1,
    seed: int = 2026,
) -> dict[str, int]:
    """Return uid -> fold_id (0..k-1 or FOLD_VALIDATION). Groups are never split."""
    rng = np.random.default_rng(seed)
    members: dict[int, list[str]] = {}
    for u, g in groups.items():
        members.setdefault(g, []).append(u)
    gids = list(members)
    rng.shuffle(gids)

    n_total = len(groups)
    fold_of_group: dict[int, int] = {}
    n_val = 0
    rest = []
    for g in gids:
        if n_val < val_fraction * n_total and len(members[g]) <= max(1, 0.01 * n_total):
            fold_of_group[g] = FOLD_VALIDATION
            n_val += len(members[g])
        else:
            rest.append(g)

    size = {g: sum(residues[u] for u in members[g]) for g in rest}
    rest.sort(key=lambda g: -size[g])  # stable: ties keep the shuffled order
    load = np.zeros(k, dtype=np.int64)
    for g in rest:
        f = int(np.argmin(load))
        fold_of_group[g] = f
        load[f] += size[g]
    return {u: fold_of_group[g] for u, g in groups.items()}


# ---------------------------------------------------------------------------
# splits.json
# ---------------------------------------------------------------------------


def load_splits(path: str | Path) -> dict:
    with open(path) as f:
        return json.load(f)


def fold_lookup(splits: dict) -> dict[str, int]:
    """uid -> fold_id for every protein that may be used (excluded proteins are absent)."""
    return {p["uid"]: p["fold_id"] for p in splits["proteins"]}


def uids_in_folds(splits: dict, folds) -> list[str]:
    folds = set(folds)
    return [p["uid"] for p in splits["proteins"] if p["fold_id"] in folds]


def cv_folds(splits: dict) -> list[int]:
    return list(range(splits["params"]["k_folds"]))


def oof_partitions(splits: dict, k: int) -> tuple[list[str], list[str]]:
    """(train uids, held-out uids) for out-of-fold round k."""
    others = [f for f in cv_folds(splits) if f != k]
    return uids_in_folds(splits, others), uids_in_folds(splits, [k])


__all__ = [
    "FOLD_CB513",
    "FOLD_VALIDATION",
    "HomologyRule",
    "assign_folds",
    "connected_groups",
    "cv_folds",
    "fold_lookup",
    "homology_edges",
    "load_splits",
    "oof_partitions",
    "read_m8",
    "uids_in_folds",
]
