import json

import numpy as np
import pytest
import torch
import torch.nn.functional as F

from e2b.baseline import (
    BaselineConfig,
    BaselineModel,
    PredictionStore,
    PredictionWriter,
    export_predictions,
    masked_loss,
    run_epoch,
)
from e2b.dataset import MockProteinDataset, collate_proteins, make_loader


def test_padded_proteins_match_individual_inference():
    torch.manual_seed(5)
    model = BaselineModel(BaselineConfig(input_dim=12, hidden_dim=6, dropout=0)).eval()
    ds = MockProteinDataset(n=2, dim=12, min_len=5, max_len=15, seed=11)
    batch = collate_proteins([ds[0], ds[1]])
    with torch.no_grad():
        together = model(batch["esm_emb"], batch["pad_mask"])
        for i in range(2):
            item = collate_proteins([ds[i]])
            alone = model(item["esm_emb"], item["pad_mask"])
            assert torch.allclose(together[i, :item["lengths"][0]], alone[0], atol=1e-6)


def test_loss_ignores_masked_residues():
    logits = torch.randn(1, 4, 8, requires_grad=True)
    target = torch.tensor([[0, -100, 1, -100]])
    expected = F.cross_entropy(logits[0, [0, 2]], target[0, [0, 2]])
    assert torch.allclose(masked_loss(logits, target), expected)
    masked_loss(logits, target).backward()
    assert torch.count_nonzero(logits.grad[0, [1, 3]]) == 0


def test_training_and_uid_indexed_export(tmp_path):
    torch.manual_seed(7)
    ds = MockProteinDataset(n=5, dim=8, min_len=5, max_len=9, seed=3)
    loader = make_loader([0], batch_size=2, dataset=ds, bucket_by_length=False)
    model = BaselineModel(BaselineConfig(input_dim=8, hidden_dim=4, dropout=0))
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    train = run_epoch(model, loader, torch.device("cpu"), optimizer)
    assert train["residues"] == sum(int(item["valid_mask"].sum()) for item in ds)
    assert np.isfinite(train["loss"])

    entries = [{"uid": item["uid"], "length": len(item["seq"]), "fold_id": 0} for item in ds]
    writer = PredictionWriter(tmp_path / "oof.npy", entries, {"kind": "test"})
    metrics = export_predictions(model, loader, writer, torch.device("cpu"))
    writer.finish()
    index = json.loads((tmp_path / "oof.index.json").read_text())
    probs = np.load(tmp_path / "oof.npy", mmap_mode="r")
    assert index["uids"] == [item["uid"] for item in ds]
    assert index["offsets"][-1] == probs.shape[0]
    assert probs.shape[1] == 8
    assert np.allclose(probs.sum(axis=1), 1, atol=1e-5)
    assert metrics["residues"] == train["residues"]
    store = PredictionStore(tmp_path / "oof.npy")
    for i, item in enumerate(ds):
        assert store.get(item["uid"]).shape == (len(item["seq"]), 8)
        assert store.fold_ids[i] == 0


def test_export_rejects_missing_or_duplicate_uids(tmp_path):
    entries = [{"uid": "a", "length": 2, "fold_id": 0},
               {"uid": "b", "length": 2, "fold_id": 1}]
    writer = PredictionWriter(tmp_path / "oof.npy", entries, {})
    batch = {"uid": ["a"], "lengths": torch.tensor([2])}
    writer.write_batch(batch, np.full((1, 2, 8), 1 / 8, np.float32))
    with pytest.raises(ValueError, match="repeated uid"):
        writer.write_batch(batch, np.full((1, 2, 8), 1 / 8, np.float32))
    with pytest.raises(ValueError, match="Missing predictions"):
        writer.finish()
