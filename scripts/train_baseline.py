"""Train the Role 3 BiLSTM baseline and export OOF, validation, and CB513 probabilities.

The validation partition selects checkpoints. CB513 is read only after the final
checkpoint is fixed. Existing completed checkpoints are reused on reruns.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e2b.baseline import (  # noqa: E402
    BaselineConfig,
    BaselineModel,
    PredictionStore,
    PredictionWriter,
    checkpoint_payload,
    confusion_metrics,
    dataset_entries,
    export_predictions,
    run_epoch,
    save_checkpoint,
)
from e2b.dataset import DEFAULT_MODEL, ProteinDataset, make_loader, make_oof_loaders  # noqa: E402
from e2b.splits import load_splits  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def fit(model, train_loader, val_loader, args, device):
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    if device.type == "cuda":
        scaler = (torch.amp.GradScaler("cuda") if hasattr(torch.amp, "GradScaler")
                  else torch.cuda.amp.GradScaler())
    else:
        scaler = None
    best_loss, best_state, best_epoch, best_validation = float("inf"), None, 0, None
    stale, history = 0, []
    for epoch in range(1, args.epochs + 1):
        start = time.perf_counter()
        train = run_epoch(model, train_loader, device, optimizer, scaler)
        with torch.no_grad():
            validation = run_epoch(model, val_loader, device)
        history.append({"epoch": epoch, "train": train, "validation": validation})
        print(f"epoch {epoch:02d} train loss={train['loss']:.4f} "
              f"val loss={validation['loss']:.4f} val Q8={validation['Q8']:.2f} "
              f"time={time.perf_counter() - start:.1f}s", flush=True)
        if validation["loss"] < best_loss - 1e-4:
            best_loss = validation["loss"]
            best_epoch = epoch
            best_validation = validation
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= args.patience:
                break
    model.load_state_dict(best_state)
    return best_epoch, best_validation, history


def load_or_fit(fold_id, train_loader, val_loader, args, device, config, split_sha256, training_settings):
    name = f"fold_{fold_id}" if isinstance(fold_id, int) else "final"
    path = args.checkpoint_dir / f"baseline_{name}.pt"
    seed = args.seed + (fold_id if isinstance(fold_id, int) else 1000)
    torch.manual_seed(seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(seed)
    model = BaselineModel(config).to(device)
    if path.exists() and not args.force:
        payload = torch.load(path, map_location="cpu", weights_only=True)
        expected = {
            "architecture": "frozen_esm2_1layer_bilstm_q8",
            "config": vars(config),
            "fold_id": fold_id,
            "seed": seed,
            "split_sha256": split_sha256,
            "esm_model": args.model,
            "training": training_settings,
        }
        for key, value in expected.items():
            if payload.get(key) != value:
                raise ValueError(f"{path} has a different {key}; use --force to retrain")
        model.load_state_dict(payload["model_state"])
        model.eval()
        print(f"reused {path}", flush=True)
    else:
        print(f"training {name}: {len(train_loader.dataset)} proteins", flush=True)
        epoch, validation, history = fit(model, train_loader, val_loader, args, device)
        payload = checkpoint_payload(model, fold_id=fold_id, seed=seed, epoch=epoch,
                                     validation=validation, split_sha256=split_sha256,
                                     esm_model=args.model)
        payload["training"] = training_settings
        payload["history"] = history
        save_checkpoint(payload, path)
        print(f"saved {path} (best epoch {epoch})", flush=True)
    return model, path, payload


def score_oof(path, records):
    store = PredictionStore(path)
    confusion = np.zeros((8, 8), dtype=np.int64)
    for record in records:
        pred = store.get(record.uid).argmax(axis=-1)
        valid = record.valid_mask
        confusion += np.bincount(
            record.q8_label[valid] * 8 + pred[valid], minlength=64).reshape(8, 8)
    return confusion_metrics(torch.from_numpy(confusion))


def args_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--splits-path", type=Path, default=ROOT / "splits/splits.json")
    parser.add_argument("--processed-dir", type=Path, default=ROOT / "data/processed")
    parser.add_argument("--emb-dir", type=Path, default=ROOT / "data/embeddings")
    parser.add_argument("--checkpoint-dir", type=Path, default=ROOT / "checkpoints")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--predict-batch-size", type=int, default=32)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--hidden-dim", type=int, default=256)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--patience", type=int, default=4)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto")
    parser.add_argument("--force", action="store_true", help="retrain and overwrite completed checkpoints")
    args = parser.parse_args()
    if args.epochs < 1 or args.patience < 1 or args.batch_size < 1 or args.predict_batch_size < 1:
        parser.error("epochs, patience, and batch sizes must be positive")
    if not 0 <= args.dropout < 1 or args.hidden_dim < 1 or args.lr <= 0 or args.weight_decay < 0:
        parser.error("invalid model or optimizer parameters")
    if args.device == "cuda" and not torch.cuda.is_available():
        parser.error("CUDA is not available")
    return args


def main():
    args = args_parser()
    device = torch.device("cuda" if args.device == "auto" and torch.cuda.is_available() else
                          "cpu" if args.device == "auto" else args.device)
    split_sha256 = hashlib.sha256(args.splits_path.read_bytes()).hexdigest()
    splits = load_splits(args.splits_path)
    k_folds = int(splits["params"]["k_folds"])
    data_args = dict(splits_path=args.splits_path, processed_dir=args.processed_dir,
                     emb_dir=args.emb_dir, model=args.model)

    # The CB513 partition is not loaded until the final model is fixed.
    all_cv = ProteinDataset(range(k_folds), **data_args)
    val = ProteinDataset([-1], **data_args)
    if all_cv.dim != val.dim:
        raise ValueError("Training and validation embedding dimensions differ")
    config = BaselineConfig(input_dim=all_cv.dim, hidden_dim=args.hidden_dim, dropout=args.dropout)
    training_settings = dict(epochs=args.epochs, patience=args.patience, lr=args.lr,
                             weight_decay=args.weight_decay, batch_size=args.batch_size)
    val_loader = make_loader([-1], args.predict_batch_size, dataset=val,
                             num_workers=args.num_workers)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    common_meta = {"architecture": "frozen_esm2_1layer_bilstm_q8", "esm_model": args.model,
                   "split_sha256": split_sha256, "model_config": vars(config)}
    oof_writer = PredictionWriter(args.output_dir / "baseline_oof.npy", dataset_entries(all_cv),
                                  {**common_meta, "kind": "out_of_fold"})
    fold_results = []
    for k in range(k_folds):
        train_loader, held_loader = make_oof_loaders(
            k, batch_size=args.batch_size, num_workers=args.num_workers,
            seed=args.seed + k, **data_args)
        model, checkpoint, payload = load_or_fit(
            k, train_loader, val_loader, args, device, config, split_sha256, training_settings)
        result = export_predictions(model, held_loader, oof_writer, device)
        fold_results.append({"fold": k, "checkpoint": checkpoint.name,
                             "best_epoch": payload["epoch"], "heldout": result})
        print(f"fold {k} held-out Q8={result['Q8']:.2f} G->H={result['G_predicted_as_H']}", flush=True)
        del model, train_loader, held_loader
        if device.type == "cuda":
            torch.cuda.empty_cache()
    oof_writer.finish()
    oof_metrics = score_oof(args.output_dir / "baseline_oof.npy", all_cv.records)

    train_loader = make_loader(range(k_folds), args.batch_size, shuffle=True,
                               num_workers=args.num_workers, seed=args.seed + 1000, **data_args)
    model, checkpoint, payload = load_or_fit(
        "final", train_loader, val_loader, args, device, config, split_sha256, training_settings)
    val_writer = PredictionWriter(args.output_dir / "baseline_val.npy", dataset_entries(val),
                                  {**common_meta, "kind": "validation", "checkpoint": checkpoint.name})
    validation = export_predictions(model, val_loader, val_writer, device)
    val_writer.finish()
    cb513 = ProteinDataset([-2], **data_args)
    if all_cv.dim != cb513.dim:
        raise ValueError("CB513 embedding dimension differs from training")
    cb_loader = make_loader([-2], args.predict_batch_size, dataset=cb513,
                            num_workers=args.num_workers)
    cb_writer = PredictionWriter(args.output_dir / "baseline_cb513.npy", dataset_entries(cb513),
                                 {**common_meta, "kind": "cb513", "checkpoint": checkpoint.name})
    cb_metrics = export_predictions(model, cb_loader, cb_writer, device)
    cb_writer.finish()
    summary = {
        **common_meta,
        "folds": fold_results,
        "oof": oof_metrics,
        "final_checkpoint": checkpoint.name,
        "final_best_epoch": payload["epoch"],
        "validation": validation,
        "cb513": cb_metrics,
        "prediction_files": ["baseline_oof.npy", "baseline_val.npy", "baseline_cb513.npy"],
    }
    (args.output_dir / "baseline_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"validation": validation, "cb513": cb_metrics}, indent=2), flush=True)


if __name__ == "__main__":
    main()
