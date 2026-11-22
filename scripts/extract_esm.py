"""Role 2: extract and cache frozen ESM-2 per-residue embeddings.

Examples:
  # timing test on 200 random training proteins (writes nothing)
  python scripts/extract_esm.py --benchmark 200

  # full cache for both sources (default model: esm2_t33_650M_UR50D, last layer)
  python scripts/extract_esm.py
  python scripts/extract_esm.py --model esm2_t30_150M_UR50D

  # check fp16 vs fp32 and batched vs single-sequence agreement
  python scripts/extract_esm.py --sanity
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from e2b.esm_embed import MODELS, EsmEmbedder, write_cache  # noqa: E402
from e2b.netsurfp import load_records  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def benchmark(emb: EsmEmbedder, recs, n: int, tokens_per_batch: int, total_residues: int):
    rng = np.random.default_rng(0)
    pick = rng.choice(len(recs), size=min(n, len(recs)), replace=False)
    seqs = [recs[i].seq for i in pick]
    list(emb.embed_many(seqs[:8], tokens_per_batch, progress=False))  # warm-up
    torch.cuda.synchronize()
    torch.cuda.reset_peak_memory_stats()
    t0 = time.time()
    for _ in emb.embed_many(seqs, tokens_per_batch, progress=False):
        pass
    torch.cuda.synchronize()
    el = time.time() - t0
    res = sum(map(len, seqs))
    rate = res / el
    out = {
        "model": emb.model_name,
        "layer": emb.layer,
        "proteins": len(seqs),
        "residues": res,
        "seconds": round(el, 2),
        "residues_per_second": round(rate),
        "peak_gpu_memory_gb": round(torch.cuda.max_memory_allocated() / 1e9, 2),
        "tokens_per_batch": tokens_per_batch,
        "estimated_minutes_all_sources": round(total_residues / rate / 60, 1),
        "estimated_cache_gb_fp16": round(total_residues * emb.dim * 2 / 1e9, 2),
        "gpu": torch.cuda.get_device_name(0),
    }
    print(json.dumps(out, indent=1))
    return out


def sanity(emb: EsmEmbedder, recs):
    """fp16 batched output vs fp32 single-sequence output on a few proteins."""
    import esm

    ref_model, _ = getattr(esm.pretrained, emb.model_name)()
    ref_model = ref_model.eval().to(emb.device)
    seqs = [r.seq for r in recs[:6]]
    batched = dict(emb.embed_many(seqs, progress=False))
    sims = []
    with torch.inference_mode():
        for i, s in enumerate(seqs):
            _, _, tok = emb.batch_converter([("x", s)])
            ref = ref_model(tok.to(emb.device), repr_layers=[emb.layer])["representations"][emb.layer][0, 1:-1].float().cpu()
            cos = torch.nn.functional.cosine_similarity(ref, batched[i], dim=-1)
            sims.append(float(cos.min()))
    print(json.dumps({"min_cosine_fp16_batched_vs_fp32_single": [round(x, 5) for x in sims]}))
    del ref_model
    return sims


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="esm2_t33_650M_UR50D", choices=list(MODELS))
    ap.add_argument("--layer", type=int, default=None, help="default: last layer")
    ap.add_argument("--sources", nargs="+", default=["train", "cb513"])
    ap.add_argument("--processed", type=Path, default=ROOT / "data/processed")
    ap.add_argument("--out", type=Path, default=ROOT / "data/embeddings")
    ap.add_argument("--tokens-per-batch", type=int, default=8192)
    ap.add_argument("--overlap", type=int, default=256)
    ap.add_argument("--fp32", action="store_true", help="run the model in float32")
    ap.add_argument("--benchmark", type=int, default=0, metavar="N")
    ap.add_argument("--sanity", action="store_true")
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    emb = EsmEmbedder(args.model, args.layer, device=device, half=not args.fp32, overlap=args.overlap)
    data = {s: load_records(args.processed / f"{s}.npz") for s in args.sources}
    total = sum(r.length for recs in data.values() for r in recs)

    if args.benchmark:
        res = benchmark(emb, data["train"], args.benchmark, args.tokens_per_batch, total)
        log = args.out / "benchmarks.jsonl"
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "a") as f:
            f.write(json.dumps(res) + "\n")
        return
    if args.sanity:
        sanity(emb, data["train"])
        return

    for source, recs in data.items():
        print(f"[{source}] {len(recs)} proteins -> {args.model} layer {emb.layer}", flush=True)
        offsets = np.concatenate([[0], np.cumsum([r.length for r in recs])])
        npy, idx = write_cache(
            emb, [r.uid for r in recs], [r.seq for r in recs], offsets, args.out, source, args.tokens_per_batch
        )
        print(f"  wrote {npy} ({npy.stat().st_size / 1e9:.2f} GB) and {idx.name}", flush=True)


if __name__ == "__main__":
    main()
