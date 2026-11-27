from types import SimpleNamespace

import numpy as np
import pytest

from e2b.baseline import PredictionStore, PredictionWriter
from e2b.specialist import (audit, candidate_positions, confusion, export_combined,
                            features, override_mask, select_rule)


def probabilities():
    return np.array([[.35, .55, .10, 0, 0, 0, 0, 0],
                     [.45, .35, .20, 0, 0, 0, 0, 0],
                     [.35, .15, .50, 0, 0, 0, 0, 0],
                     [.03, .95, .02, 0, 0, 0, 0, 0]], dtype=np.float32)


def test_candidates_and_window_padding():
    p = probabilities()
    np.testing.assert_array_equal(candidate_positions(p), [0, 1])
    emb = np.arange(8).reshape(4, 2)
    x = features(emb, p, np.array([0, 3]))
    for i, pos in enumerate([0, 3]):
        assert np.array_equal(x[i, :2], emb[pos])
        neighbors = x[i, 2:].reshape(5, 8)
        for j in range(5):
            index = pos + j - 2
            np.testing.assert_array_equal(neighbors[j], p[index] if 0 <= index < 4 else np.zeros(8))


def test_eligibility_and_thresholds():
    p = probabilities()
    ix = candidate_positions(p)
    rule = {"enabled": True, "max_gap": .3, "g_threshold": .6, "h_threshold": .7}
    def flips(score):
        return override_mask(p[ix].argmax(1), abs(p[ix, 0] - p[ix, 1]), np.array(score), rule)
    np.testing.assert_array_equal(flips([.65, .25]), [True, True])
    np.testing.assert_array_equal(flips([.55, .4]), [False, False])


def test_export_and_masked_audit(tmp_path):
    p = probabilities()
    entries = [{"uid": "a", "length": 2, "fold_id": -1}, {"uid": "b", "length": 2, "fold_id": -1}]
    holder = PredictionWriter(tmp_path / "baseline.npy", entries, {"kind": "validation"})
    for i, e in enumerate(entries):
        holder.write_batch({"uid": [e["uid"]], "lengths": [2]}, p[None, 2*i:2*i + 2])
    holder.finish()
    store = PredictionStore(tmp_path / "baseline.npy")
    data = SimpleNamespace(store=store, entries=entries, positions=np.array([0, 1]),
                           truth=np.array([0, -1]), base_labels=np.array([1, 0]), gaps=np.array([.2, .1]),
                           baseline_confusion=confusion(np.array([0, 2, 1]), np.array([1, 2, 1])))
    rule, _ = select_rule(data, np.array([.9, .1]))
    assert rule["enabled"]
    flip = np.ones(2, dtype=bool)
    info = audit(data, flip)
    assert (info["corrected"], info["harmed"], info["valid_overrides"]) == (1, 0, 1)
    export_combined(data, flip, tmp_path / "combined.npy", {})
    out = PredictionStore(tmp_path / "combined.npy")
    np.testing.assert_array_equal(out.array.argmax(1), [0, 1, 2, 1])
    np.testing.assert_array_equal(out.get("b"), store.get("b"))
    np.testing.assert_array_equal(out.array[:, 2:], p[:, 2:])
    np.testing.assert_array_equal(store.array, p)
    store.metadata["kind"] = "cb513"
    with pytest.raises(ValueError, match="validation"):
        select_rule(data, np.array([.8, .3]))


def test_rule_can_keep_baseline():
    data = SimpleNamespace(store=SimpleNamespace(metadata={"kind": "validation"}),
                           truth=np.array([1]), base_labels=np.array([1]), gaps=np.array([.1]),
                           baseline_confusion=confusion(np.array([0, 1]), np.array([0, 1])))
    rule, _ = select_rule(data, np.array([.92]))
    assert not rule["enabled"]
