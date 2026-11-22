import numpy as np
import pytest

from e2b.constants import AA_ORDER, LABEL_MISSING, Q8_ORDER
from e2b.netsurfp import load_records, parse_netsurfp, save_records


def _make_sample(seq, q8, max_len=12, ordered=None, evalm=None):
    arr = np.zeros((max_len, 68))
    for i, (aa, ss) in enumerate(zip(seq, q8)):
        if aa != "X":
            arr[i, AA_ORDER.index(aa)] = 1
        arr[i, 50] = 1
        arr[i, 51] = 1 if ordered is None else ordered[i]
        arr[i, 52] = 1 if evalm is None else evalm[i]
        if ss != "-":
            arr[i, 57 + Q8_ORDER.index(ss)] = 1
    return arr


@pytest.fixture
def fake_npz(tmp_path):
    data = np.stack(
        [
            _make_sample("ACDXW", "GHH-C", ordered=[1, 1, 0, 1, 1]),
            _make_sample("MKLV", "EEBT", evalm=[0, 1, 1, 0]),
        ]
    )
    path = tmp_path / "fake.npz"
    np.savez_compressed(path, pdbids=np.array(["1abc-A", "2xyz-B"]), data=data)
    return path


def test_decode_and_labels(fake_npz):
    recs = parse_netsurfp(fake_npz, "train")
    assert [r.seq for r in recs] == ["ACDXW", "MKLV"]
    assert recs[0].q8_label.tolist() == [0, 1, 1, LABEL_MISSING, 7]
    assert recs[1].q8_label.tolist() == [4, 4, 3, 6]
    assert recs[0].uid == "train_00000" and recs[1].pdb_id == "2xyz-B"


def test_valid_mask_rules(fake_npz):
    train = parse_netsurfp(fake_npz, "train")
    # disordered residue 2 and unlabeled residue 3 are excluded
    assert train[0].valid_mask.tolist() == [True, True, False, False, True]
    # eval mask is ignored outside CB513
    assert train[1].valid_mask.all()
    cb = parse_netsurfp(fake_npz, "cb513")
    assert cb[1].valid_mask.tolist() == [False, True, True, False]


def test_roundtrip(fake_npz, tmp_path):
    recs = parse_netsurfp(fake_npz, "cb513")
    out = tmp_path / "proc.npz"
    save_records(recs, out)
    back = load_records(out)
    for a, b in zip(recs, back):
        assert a.uid == b.uid and a.seq == b.seq and a.pdb_id == b.pdb_id
        assert np.array_equal(a.q8_label, b.q8_label)
        assert np.array_equal(a.valid_mask, b.valid_mask)
        assert np.array_equal(a.eval_mask, b.eval_mask)


def test_non_contiguous_mask_rejected(tmp_path):
    arr = _make_sample("ACD", "HHH")[None]
    arr[0, 1, 50] = 0  # hole in the sequence mask
    path = tmp_path / "bad.npz"
    np.savez_compressed(path, pdbids=np.array(["1bad-A"]), data=arr)
    with pytest.raises(ValueError, match="contiguous"):
        parse_netsurfp(path, "train")
