"""Role 1, step 3a: homology search with MMseqs2 (CB513 vs Train, Train vs Train).

Calls MMseqs2 modules directly (createdb / prefilter / align / convertalis), so
the Windows build does not need its BusyBox shell helpers or admin rights.

Outputs (data/similarity/):
  cb513_vs_train.m8   query=cb513 uid, target=train uid
  train_vs_train.m8   all-vs-all within Train (self hits included; ignored later)
Columns: query, target, fident, alnlen, qlen, tlen, evalue, bits

Usage:
  python scripts/run_mmseqs.py [--mmseqs tools/mmseqs/bin/mmseqs.exe] [--sensitivity 7.5]
Then rebuild the folds:
  python scripts/make_splits.py
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORMAT = "query,target,fident,alnlen,qlen,tlen,evalue,bits"


def run(cmd, cwd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), cwd=cwd, check=True)


def main():
    default_bin = ROOT / "tools/mmseqs/bin/mmseqs.exe"
    ap = argparse.ArgumentParser()
    ap.add_argument("--mmseqs", default=str(default_bin) if default_bin.exists() else (shutil.which("mmseqs") or "mmseqs"))
    ap.add_argument("--processed", type=Path, default=ROOT / "data/processed")
    ap.add_argument("--out", type=Path, default=ROOT / "data/similarity")
    ap.add_argument("--sensitivity", type=float, default=7.5)
    ap.add_argument("--max-seqs", type=int, default=2000)
    ap.add_argument("--evalue", type=float, default=1.0, help="report threshold; make_splits applies the real rule")
    args = ap.parse_args()

    out = args.out
    (out / "db").mkdir(parents=True, exist_ok=True)
    (out / "tmp").mkdir(exist_ok=True)
    mm = args.mmseqs
    for name in ("train", "cb513"):
        fasta = (args.processed / f"{name}.fasta").resolve()
        run([mm, "createdb", fasta, f"db/{name}", "--dbtype", "1"], out)

    for q, t, tag in (("cb513", "train", "cb513_vs_train"), ("train", "train", "train_vs_train")):
        run([mm, "prefilter", f"db/{q}", f"db/{t}", f"db/{tag}_pref", "-s", args.sensitivity,
             "--max-seqs", args.max_seqs], out)
        run([mm, "align", f"db/{q}", f"db/{t}", f"db/{tag}_pref", f"db/{tag}_aln", "-e", args.evalue,
             "--alignment-mode", "3"], out)
        run([mm, "convertalis", f"db/{q}", f"db/{t}", f"db/{tag}_aln", f"{tag}.m8",
             "--format-output", FORMAT], out)
        print(f"wrote {out / (tag + '.m8')}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
