import json
import pickle
from pathlib import Path

import numpy as np
import pytest
import torch

from e2b.constants import IGNORE_INDEX
from e2b.dataset import (
    LengthBucketSampler,
    MockProteinDataset,
    ProteinDataset,
    collate_proteins,
    make_loader,
)
from e2b.esm_embed import cache_paths
from e2b.netsurfp import ProteinRecord, save_records

ROOT = Path(__file__).resolve().parents[1]
D = 8
MODEL, LAYER = "esm2_t6_8M_UR50D", 6


def _rec(source, i, L, rng):
    valid = rng.random(L) > 0.2
    return ProteinRecord(
        uid=f"{source}_{i:05d}", pdb_id=f"{i}xyz-A", source=source, npz_index=i,
        seq="".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=L)),
        q8_label=rng.integers(0, 8, size=L).astype(np.int8),
        disorder_mask=valid.copy(), eval_mask=np.ones(L, bool), valid_mask=valid,
    )


@pytest.fixture
def fake_env(tmp_path):
    rng = np.random.default_rng(0)
    proc, emb = tmp_path / "processed", tmp_path / "emb"
    proc.mkdir()
    emb.mkdir()
    folds = {}
    for source, n in [("train", 12), ("cb513", 4)]:
        recs = [_rec(source, i, int(rng.integers(5, 40)), rng) for i in range(n)]
        save_records(recs, proc / f"{source}.npz")
        offsets = np.concatenate([[0], np.cumsum([r.length for r in recs])])
        npy, idx = cache_paths(emb, source, MODEL, LAYER)
        # row r of the cache = (protein index, position) so alignment is checkable
        arr = np.zeros((offsets[-1], D), np.float16)
        for i, r in enumerate(recs):
            arr[offsets[i]:offsets[i + 1], 0] = i
            arr[offsets[i]:offsets[i + 1], 1] = np.arange(r.length)
        np.save(npy, arr)
        idx.write_text(json.dumps({"uids": [r.uid for r in recs], "offsets": offsets.tolist(), "dim": D}))
        for i, r in enumerate(recs):
            folds[r.uid] = -2 if source == "cb513" else (-1 if i < 2 else i % 3)
    splits = {"params": {"k_folds": 3}, "proteins": [{"uid": u, "fold_id": f} for u, f in folds.items()]}
    (tmp_path / "splits.json").write_text(json.dumps(splits))
    return dict(splits_path=tmp_path / "splits.json", processed_dir=proc, emb_dir=emb, model=MODEL, layer=LAYER)


def test_dataset_items_follow_contract(fake_env):
    ds = ProteinDataset([0, 1], **fake_env)
    assert len(ds) == 6 and ds.dim == D  # train_00003..00011 with i % 3 in {0, 1}
    for k in range(len(ds)):
        it = ds[k]
        L = len(it["seq"])
        assert it["esm_emb"].shape == (L, D) and it["esm_emb"].dtype == torch.float32
        assert it["q8_label"].dtype == torch.long and it["q8_label"].shape == (L,)
        assert it["valid_mask"].dtype == torch.bool and it["fold_id"] in (0, 1)
        idx = int(it["uid"].split("_")[1])
        assert torch.all(it["esm_emb"][:, 0] == idx)  # rows belong to this protein
        assert torch.equal(it["esm_emb"][:, 1], torch.arange(L, dtype=torch.float32))


def test_fold_selection(fake_env):
    assert {x["fold_id"] for x in ProteinDataset([-1], **fake_env)} == {-1}
    test = ProteinDataset([-2], **fake_env)
    assert len(test) == 4 and all(x["uid"].startswith("cb513") for x in test)


def test_collate_padding_and_target(fake_env):
    ds = ProteinDataset([0, 1, 2], **fake_env)
    b = collate_proteins([ds[0], ds[1], ds[2]])
    Lmax = max(len(ds[i]["seq"]) for i in range(3))
    assert b["esm_emb"].shape == (3, Lmax, D)
    assert (b["target"][~b["valid_mask"]] == IGNORE_INDEX).all()
    assert torch.equal(b["target"][b["valid_mask"]], b["q8_label"][b["valid_mask"]])
    assert not b["valid_mask"][~b["pad_mask"]].any()
    assert b["pad_mask"].sum(1).tolist() == b["lengths"].tolist()
    # the loss works directly with ignore_index
    logits = torch.randn(3, Lmax, 8, requires_grad=True)
    loss = torch.nn.functional.cross_entropy(logits.transpose(1, 2), b["target"], ignore_index=IGNORE_INDEX)
    loss.backward()


def test_loader_covers_every_protein_once(fake_env):
    dl = make_loader([0, 1, 2], batch_size=3, shuffle=True, **fake_env)
    seen = [u for batch in dl for u in batch["uid"]]
    assert sorted(seen) == sorted(r.uid for r in dl.dataset.records)


def test_store_pickles_without_data(fake_env):
    ds = ProteinDataset([0], **fake_env)
    _ = ds[0]  # opens the memmap
    blob = pickle.dumps(ds)
    assert len(blob) < 200_000
    assert torch.equal(pickle.loads(blob)[0]["esm_emb"], ds[0]["esm_emb"])


def test_bucket_sampler():
    lengths = np.random.default_rng(1).integers(20, 900, size=1000)
    s = LengthBucketSampler(lengths, batch_size=16, shuffle=True)
    batches = list(s)
    assert len(batches) == len(s)
    assert sorted(i for b in batches for i in b) == list(range(1000))


def test_mock_dataset_contract():
    ds = MockProteinDataset(n=10, dim=32, max_len=50)
    b = collate_proteins([ds[i] for i in range(4)])
    assert b["esm_emb"].shape[-1] == 32 and b["target"].dtype == torch.long


@pytest.mark.skipif(not (ROOT / "splits/splits.json").exists(), reason="real splits not built")
def test_real_data_if_present():
    try:
        ds = ProteinDataset([-2])
    except FileNotFoundError:
        pytest.skip("embedding cache not built")
    assert len(ds) == 513
    item = ds[0]
    assert item["esm_emb"].shape == (len(item["seq"]), ds.dim)
    assert torch.isfinite(item["esm_emb"]).all()
