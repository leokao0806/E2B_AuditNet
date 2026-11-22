import numpy as np
import pytest

from e2b.esm_embed import assign_positions, plan_windows


def test_short_sequence_single_window():
    assert plan_windows(500) == [(0, 500)]
    assert plan_windows(1022) == [(0, 1022)]


@pytest.mark.parametrize("length", [1023, 1500, 1632, 2044, 3000, 5000])
def test_long_sequence_windows_cover_and_fit(length):
    wins = plan_windows(length, max_len=1022, overlap=256)
    assert wins[0][0] == 0 and wins[-1][1] == length
    assert all(e - s == 1022 for s, e in wins)
    for (s0, e0), (s1, e1) in zip(wins, wins[1:]):
        assert s1 < e0, "consecutive windows must overlap"
    which = assign_positions(length, wins)
    assert which.shape == (length,)
    # every residue is assigned to a window that contains it
    for p in range(length):
        s, e = wins[which[p]]
        assert s <= p < e


def test_make_batches_respects_budgets():
    from e2b.esm_embed import EsmEmbedder

    lengths = list(np.random.default_rng(0).integers(20, 1023, size=3000))
    batches = EsmEmbedder.make_batches(lengths, tokens_per_batch=8192, attn_per_batch=4 * 1024**2)
    assert sorted(i for b in batches for i in b) == list(range(3000))
    for b in batches:
        t = max(lengths[i] for i in b) + 2
        assert len(b) == 1 or (len(b) * t <= 8192 and len(b) * t * t <= 4 * 1024**2)


def test_assignment_prefers_central_context():
    wins = plan_windows(1632, max_len=1022, overlap=256)
    which = assign_positions(1632, wins)
    # residue 0 can only come from the first window, the last from the last window
    assert which[0] == 0 and which[-1] == len(wins) - 1
    # a residue near the end of window 0 but inside window 1 should use window 1
    s1 = wins[1][0]
    assert which[wins[0][1] - 1] == 1
    # min two-sided context of the chosen window is never worse than any alternative
    for p in range(1632):
        ctx = [min(p - s, e - 1 - p) for s, e in wins if s <= p < e]
        s, e = wins[which[p]]
        assert min(p - s, e - 1 - p) == max(ctx)
    assert s1 > 0
