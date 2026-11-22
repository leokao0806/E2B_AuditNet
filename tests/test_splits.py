import numpy as np

from e2b.constants import FOLD_VALIDATION
from e2b.splits import HomologyRule, assign_folds, connected_groups, homology_edges, read_m8


def test_connected_groups_transitive():
    uids = ["a", "b", "c", "d", "e"]
    g = connected_groups(uids, {("a", "b"), ("b", "c")})
    assert g["a"] == g["b"] == g["c"]
    assert len({g["a"], g["d"], g["e"]}) == 3


def test_homology_rule_and_m8(tmp_path):
    p = tmp_path / "hits.m8"
    p.write_text(
        "a\tb\t0.40\t100\t120\t130\t1e-20\t150\n"  # passes
        "a\tc\t0.20\t100\t120\t90\t1e-10\t60\n"  # identity too low
        "a\td\t0.90\t10\t120\t15\t1e-5\t30\n"  # alignment too short
        "a\te\t0.50\t60\t120\t80\t0.5\t20\n"  # not significant
        "a\ta\t1.00\t120\t120\t120\t0\t300\n"  # self hit ignored
    )
    edges = homology_edges(read_m8(p), HomologyRule())
    assert edges == {("a", "b")}


def test_coverage_blocks_tag_only_matches(tmp_path):
    p = tmp_path / "tags.m8"
    p.write_text(
        "a\tb\t1.00\t22\t300\t250\t1e-6\t45\n"  # shared His-tag only: 22/250 coverage
        "a\tc\t0.35\t200\t300\t260\t1e-30\t200\n"  # real homolog: 200/260 coverage
    )
    hits = read_m8(p)
    assert homology_edges(hits, HomologyRule()) == {("a", "b"), ("a", "c")}
    assert homology_edges(hits, HomologyRule(min_fident=0.3, min_cov_shorter=0.5)) == {("a", "c")}


def test_assign_folds_keeps_groups_together_and_balances():
    rng = np.random.default_rng(0)
    uids = [f"p{i}" for i in range(2000)]
    edges = {(f"p{i}", f"p{i + 1}") for i in range(0, 300, 3)}  # 100 pairs
    groups = connected_groups(uids, edges)
    residues = {u: int(rng.integers(30, 600)) for u in uids}
    folds = assign_folds(groups, residues, k=5, val_fraction=0.1, seed=1)

    for a, b in edges:
        assert folds[a] == folds[b]
    n_val = sum(1 for f in folds.values() if f == FOLD_VALIDATION)
    assert 190 <= n_val <= 215
    load = np.zeros(5)
    for u, f in folds.items():
        if f >= 0:
            load[f] += residues[u]
    assert load.max() / load.min() < 1.01
    # deterministic for a fixed seed
    assert folds == assign_folds(groups, residues, k=5, val_fraction=0.1, seed=1)
