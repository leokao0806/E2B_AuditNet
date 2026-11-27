"""Fit the H/G specialist on baseline OOF predictions, then freeze and export."""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e2b.baseline import save_checkpoint
from e2b.constants import G_IDX, H_IDX
from e2b.dataset import ROOT
from e2b.specialist import (HGSpecialist, RADIUS, MAX_GAP, audit, export_combined,
                            override_mask, predict_g, read_split, select_rule)


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=2) + "\n")
    tmp.replace(path)


def fit(train, val, seed):
    torch.manual_seed(seed)
    model = HGSpecialist()
    opt = torch.optim.AdamW(model.parameters(), lr=0.001)
    lab_train = (train.truth == G_IDX) | (train.truth == H_IDX)
    lab_val = (val.truth == G_IDX) | (val.truth == H_IDX)
    x = train.inputs[lab_train]
    y = torch.tensor(train.truth[lab_train] == G_IDX, dtype=torch.float32)
    vx = val.inputs[lab_val]
    vy = torch.tensor(val.truth[lab_val] == G_IDX, dtype=torch.float32)
    if len(y.unique()) != 2 or len(vy.unique()) != 2:
        raise ValueError("Both H and G examples are required in training and validation")
    history, best_state, best_loss, best_epoch = [], None, float("inf"), 0
    for epoch in range(1, 21):
        model.train()
        for batch in torch.randperm(len(x)).split(512):
            opt.zero_grad()
            loss = F.binary_cross_entropy_with_logits(model(x[batch]), y[batch])
            loss.backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            val_loss = F.binary_cross_entropy_with_logits(model(vx), vy).item()
        history.append(val_loss)
        if val_loss < best_loss:
            best_loss, best_epoch = val_loss, epoch
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        print(f"Epoch {epoch}: validation binary loss {val_loss:.5f}", flush=True)
        if epoch - best_epoch >= 4:
            break
    if best_state is None:
        raise ValueError("No finite validation loss")
    model.load_state_dict(best_state)
    return model, {"best_epoch": best_epoch, "validation_loss_by_epoch": history,
                   "training_proteins_in_oof_store": len(train.entries), "training_H": int((y == 0).sum()),
                   "training_G": int((y == 1).sum()), "validation_HG_for_loss": len(vy)}


def evaluate_and_export(model, data, rule, dest, kind, checkpoint_name):
    p_g = predict_g(model, data.inputs)
    flip = override_mask(data.base_labels, data.gaps, p_g, rule)
    metadata = {**data.store.metadata, "kind": kind, "architecture": "ESM2_plus_Q8_window_HG_specialist",
                "baseline_checkpoint": data.store.metadata.get("checkpoint"), "checkpoint": checkpoint_name,
                "rule": rule, "window_radius": RADIUS, "score_semantics": "baseline H/G scores swapped on override"}
    metadata["baseline_model_config"] = metadata.pop("model_config", {})
    metadata["model_config"] = {"hidden_dim": 64, "dropout": 0.2, "embedding_dim": 1280, "window_radius": RADIUS}
    export_combined(data, flip, dest / f"auditnet_{kind}.npy", metadata)
    stats = {"baseline": audit(data, np.zeros_like(flip)), "combined": audit(data, flip),
             "all_overrides": int(flip.sum()), "all_candidate_residues": len(p_g)}
    write_json(dest / f"{kind}_metrics.json", stats)
    print(kind, json.dumps(stats), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline-dir", type=Path, default=ROOT / "outputs")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "outputs/role5")
    ap.add_argument("--checkpoint", type=Path, default=ROOT / "checkpoints/hg_specialist.pt")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--export-cb513", action="store_true", help="Use the saved checkpoint and locked rule only")
    args = ap.parse_args()
    if args.export_cb513:
        ckpt = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
        data = read_split("cb513", args.baseline_dir)
        if (data.store.metadata["split_sha256"] != ckpt["split_sha256"]
                or data.store.metadata.get("checkpoint") != ckpt["baseline_checkpoint"]):
            raise ValueError("CB513 baseline does not match the validation baseline")
        model = HGSpecialist()
        model.load_state_dict(ckpt["model_state"])
        evaluate_and_export(model, data, ckpt["rule"], args.output_dir, "cb513", args.checkpoint.name)
        return
    train = read_split("oof", args.baseline_dir)
    val = read_split("val", args.baseline_dir)
    model, fit_info = fit(train, val, args.seed)
    rule, trials = select_rule(val, predict_g(model, val.inputs))
    config = {"rule": rule, "window_radius": RADIUS, "candidate_max_gap": MAX_GAP,
              "criterion": "highest validation G F1 with nonnegative net corrections; ties use Q8 then fewer overrides",
              "seed": args.seed, "split_sha256": train.store.metadata["split_sha256"],
              "baseline_checkpoint": val.store.metadata["checkpoint"]}
    save_checkpoint({**config, "model_state": model.state_dict(), "best_epoch": fit_info["best_epoch"]}, args.checkpoint)
    write_json(args.output_dir / "routing.json", config)
    write_json(args.output_dir / "validation_trials.json", {"training": fit_info, "candidates": trials})
    evaluate_and_export(model, val, rule, args.output_dir, "val", args.checkpoint.name)


if __name__ == "__main__":
    main()
