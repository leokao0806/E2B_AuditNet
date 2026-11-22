# E2B-AuditNet: Protein Secondary Structure H/G Correction

## Team Setup & Workflow (Start Here — Roles 3–6)

GitHub is the single source of truth for code, splits (`splits/`), and processed targets (`data/processed/`). Only the large ESM-2 embedding files (`data/embeddings/*.npy` + `*.index.json`) are stored on Google Drive.

### 1. Setup (~10 minutes + Drive download)

1. **Clone this repository:**
   ```bash
   git clone https://github.com/leokao0806/E2B_AuditNet.git
   cd E2B_AuditNet
   ```
2. **Download ESM-2 embeddings from Google Drive:**
   - Download **only the `embeddings_650M/` folder** (4 files, ~7.48 GB; default model `D = 1280`).
   - *(Optional)* Download `embeddings_150M_optional/` (~3.74 GB, `D = 640`) only if we run a 150M ablation.
   - **Do not** download `E2B_AuditNet_code.zip` from Drive — use this GitHub repo for code.
   - Place the downloaded files directly into `data/embeddings/` (keep the original filenames).
   - *(Colab users)* Copy `data/embeddings/` from Drive to local Colab disk (e.g., `/content/`) before training; reading a 7 GB memmap directly from a mounted Drive is much slower.
3. **Install dependencies** (Python 3.10+; install PyTorch for your hardware first from https://pytorch.org/get-started/locally/):
   ```bash
   pip install -r requirements.txt
   ```
4. **Verify installation:**
   ```bash
   python -m pytest
   ```
   - Expected: **`25 passed`**. If you see `24 passed, 1 skipped`, the embedding files are not in `data/embeddings/`.

### 2. What Each Role Does Next

- **Role 3 (Baseline & OOF):**
  - Implement **Frozen ESM-2 (1280-dim) -> 1-layer BiLSTM -> Linear Q8 classifier** using `make_oof_loaders(k=0..4)` and `F.cross_entropy(..., ignore_index=-100)`.
  - Train 5 folds (`k=0..4`, early-stopping on `val_dl = make_loader([-1])`) and save **Out-of-Fold (OOF) Softmax probabilities** (`[L, 8]`) for every training protein keyed by **`uid`** (e.g., `train_00001`, **not** `pdb_id`).
  - Also export final Baseline predictions on `val (-1)` and `CB513 (-2)` for Roles 5 & 6.
- **Role 4 (Class-Weighted Baseline):**
  - Reuse Role 3's exact architecture and folds, apply class-weighted cross-entropy loss tuned on `val (-1)`, and export predictions on `CB513 (-2)` for Role 6.
- **Role 5 (H/G Specialist & Routing):**
  - Train the H/G Specialist on frozen ESM-2 embeddings + Role 3's OOF probabilities + local sequence window features; tune ambiguous H/G routing thresholds strictly on `val (-1)` and export combined **E2B-AuditNet** predictions on `CB513 (-2)`.
- **Role 6 (Evaluation & Audit):**
  - Evaluate all three systems (*Baseline* vs. *Class-Weighted* vs. *E2B-AuditNet*) on `CB513 (-2)` using `valid_mask` (already applies disorder + CB513 eval mask): overall Q8/Q3, G/H Precision/Recall/F1, 8x8 confusion matrices, and SecNet-style H/G boundary & net-correction audit.
- **Git hygiene:** Save large model weights and OOF/prediction dumps under `checkpoints/` or `outputs/` (both are gitignored).

---

## Data Pipeline, ESM-2 Cache & Shared DataLoader (Roles 1 and 2)

This folder holds the Role 1 (Data & Partitioning) and Role 2 (ESM-2 & Infrastructure)
deliverables from the task allocation proposal:

| Deliverable | Where |
|---|---|
| Dataset parser (NetSurfP-2.0 `.npz` → per-protein records) | `e2b/netsurfp.py`, `scripts/prepare_data.py` |
| Mask / label / length audit | `data/processed/data_audit.json` |
| Similarity-safe fold file | **`splits/splits.json`** |
| Residue and class statistics per split (incl. H/G) | **`splits/split_stats.md`** / `.json` |
| Frozen ESM-2 extractor with >1,022-residue windowing | `e2b/esm_embed.py`, `scripts/extract_esm.py` |
| Embedding cache (float16 memmap) | `data/embeddings/` (shared on the drive, not in git) |
| Shared `Dataset` / padded `DataLoader` / mock data | **`e2b/dataset.py`** |
| Tests | `tests/` (`python -m pytest`) |

## 1. Quick start for Roles 3–6

```python
from e2b.dataset import make_loader, make_oof_loaders, MockProteinDataset
import torch.nn.functional as F

# Out-of-fold round k (Role 3): train on the other 4 folds, predict fold k
train_dl, heldout_dl = make_oof_loaders(k=0, batch_size=16)

val_dl  = make_loader([-1], batch_size=32)   # validation: class weights (Role 4), routing thresholds (Role 5)
test_dl = make_loader([-2], batch_size=32)   # CB513: final evaluation only

for batch in train_dl:
    x = batch["esm_emb"]                     # [B, Lmax, 1280] float32
    logits = model(x, batch["pad_mask"])     # [B, Lmax, 8]
    loss = F.cross_entropy(logits.transpose(1, 2), batch["target"], ignore_index=-100)
```

Until you have the cache, develop against `MockProteinDataset(n=64, dim=1280)`; it returns the
same dictionaries. `make_loader(..., dataset=MockProteinDataset(...))` also works.

### Shared Data Contract (as implemented)

Per protein (`ProteinDataset[i]`):

| Key | Type / shape | Meaning |
|---|---|---|
| `pdb_id` | str | chain id from `pdbids` (**not unique in CB513**: 513 regions from 434 chains) |
| `uid` | str | unique key: `train_NNNNN` / `cb513_NNNNN` = row in the original `.npz` |
| `seq` | str, length L | decoded from channels [0:20]; unknown residues = `X` |
| `esm_emb` | FloatTensor [L, D] | frozen ESM-2 last layer (D = 1280 for 650M, 640 for 150M) |
| `q8_label` | LongTensor [L] | 0..7 = G, H, I, B, E, S, T, C (0 = G, 1 = H) |
| `valid_mask` | BoolTensor [L] | ordered (`[51]`) and labelled; on CB513 also the evaluation mask (`[52]`) |
| `fold_id` | int | 0..4 CV/OOF folds, −1 validation, −2 CB513 |

A batch (`collate_proteins`) adds padding to the longest protein, plus:
`target` (= `q8_label` where `valid_mask`, else −100), `pad_mask` (True on real residues),
`lengths`, and lists `pdb_id` / `uid` / `seq`. Padding is −1 in `q8_label` and False in `valid_mask`.
Always score with `valid_mask` (Role 6). Batches are bucketed by length to reduce padding;
pass `bucket_by_length=False` for plain batching.

Q8 → Q3 mapping (NetSurfP-2.0): G, H, I → H; B, E → E; S, T, C → C (`e2b.constants.Q8_TO_Q3`).

## 2. Data (Role 1)

Source: NetSurfP-2.0 "Training/Test Data" tab — `Train_HHblits.npz` (433 MB; float64 array
10,848 × 1,632 × 68, 9.6 GB uncompressed) and `CB513_HHblits.npz` (23 MB; 513 × 874 × 68).
SHA-256 values are in `data/processed/data_audit.json` and `splits/splits.json`.

Channels used: [0:20] one-hot residues, [50] sequence mask, [51] disorder mask, [52] evaluation
mask, [57:65] Q8 labels (GHIBESTC). HMM profiles and the structural channels are ignored.

Checks done by `scripts/prepare_data.py`:

- **Residue order** of [0:20] is `ACDEFGHIKLMNPQRSTVWY` (not stated on the NetSurfP page). Verified
  by decoding CB513 entries 154l-A and 1a45-A and matching the RCSB FASTA exactly.
- Sequence mask is always a contiguous prefix; every residue has exactly one Q8 bit; no residue
  has two amino-acid bits.
- Train: 2,778,438 residues, 16,799 `X`, 179,591 disordered (they still carry labels, so the
  disorder mask must be applied), 36 proteins longer than 1,022 (max 1,632).
- CB513: 144,360 residues, of which 84,336 are inside the evaluation mask; 83,772 scored residues
  after the disorder mask. Entries `cb513_00299` (1rsy-A) and `cb513_00377` (2alp-A) have an
  all-zero evaluation mask, so they contribute no scored residues.

### Leakage check and folds

NetSurfP-2.0 states that no training protein has more than 25% identity to CB513. We checked this
ourselves with exact/substring matching and MMseqs2 (`-s 7.5`, local alignments):

| CB513 vs Train (E-value ≤ 1e-3, alignment ≥ 20) | Pairs | CB513 chains | Train chains |
|---|---:|---:|---:|
| identity ≥ 25% | 470 | 203 | 323 |
| identity ≥ 50% | 72 | 27 | 45 |
| identity ≥ 90% | 42 | 13 | 29 |

Two CB513 chains are identical to Train chains (1tii-C; 1wfb-B vs 1wfb-A), and one is a substring
of a Train chain (1mrt-A in 4mt2-A). **All 324 Train proteins that match CB513 under the 25% rule
are excluded** (listed with reasons in `splits/split_stats.md`), so CB513 stays held out. This is
stricter than E2BNet, which removed exact duplicates only.

Inside Train, homologous proteins are put in the same partition (connected components of the
homology graph). The rule is E-value ≤ 1e-3, identity ≥ 30%, and an alignment covering ≥ 50% of
the shorter protein. The coverage condition matters: without it, ~72,000 "high-identity" edges
come from shared 20–22 residue expression tags (e.g. `MGSSHHHHHHSSGLVPRGSH`) and chain unrelated
proteins into one 2,138-protein group. With it there are 9,376 groups (largest: 713 proteins).

Validation (−1) is ~10% of Train, drawn as whole groups at random; the rest is spread over five
folds with equal numbers of scored residues (seed 2026).

| Partition | Proteins | Scored residues | G % | H % |
|---|---:|---:|---:|---:|
| folds 0–4 (each) | 1,831–2,124 | 450,320–450,321 | 3.78–3.91 | 33.8–34.8 |
| validation (−1) | 1,053 | 245,405 | 3.86 | 34.88 |
| CB513 (−2) | 513 | 83,772 | 3.70 | 30.69 |
| excluded | 324 | — | — | — |

About 21% of G residues sit directly next to an H segment in every partition (18% on CB513);
full tables are in `splits/split_stats.md`.

`splits.json` layout: `params`, `homology` (rules and counts), `excluded` (uid, pdb_id, reason),
and `proteins` (uid, pdb_id, source, fold_id, group). Helpers: `e2b.splits.load_splits`,
`fold_lookup`, `uids_in_folds`, `oof_partitions`.

## 3. ESM-2 features (Role 2)

Model: `esm2_t33_650M_UR50D` (fair-esm 2.0.0), weights frozen, last layer (33), D = 1280.
`esm2_t30_150M_UR50D` (the E2BNet backbone, D = 640) is cached too, so the team can compare them.
Pass `model="esm2_t30_150M_UR50D"` to `make_loader` to use it.

- **Full-length chains.** E2BNet truncates to 128 residues; we keep every residue. Chains of up to
  1,022 residues are encoded in one pass. The 36 longer chains (max 1,632) are split into
  1,022-residue windows overlapping by 256. Each residue takes its vector from the window where it
  has the most context on both sides, so no averaging is involved (`plan_windows`,
  `assign_positions`, tested in `tests/test_windows.py`).
- **Precision.** The forward pass runs in fp16 and the cache is stored in fp16. Compared with fp32
  single-sequence inference, the minimum per-residue cosine similarity is 0.9997.
- **Memory.** Batches are capped by tokens (8,192) and by attention size (B·L² ≤ 4·1024²), because
  fair-esm materialises float32 attention maps. On out-of-memory errors the batch is split in
  half. Peak use is 2.1 GB on the longest chains.
- **Format.** One memory-mapped array per source and model, e.g.
  `data/embeddings/train.esm2_t33_650M_UR50D.L33.f16.npy` with shape [2,778,438, 1280], plus
  `*.index.json` (uids and row offsets). The rows of protein i are `offsets[i]:offsets[i+1]`.
  Re-running is byte-identical.

| | 650M (default) | 150M |
|---|---:|---:|
| Extraction, Train + CB513 (RTX 5060 Laptop 8 GB) | 5.5 min | 2.8 min |
| Throughput | ~8,900 residues/s | ~17,000 residues/s |
| Cache size (fp16) | 7.48 GB | 3.74 GB |
| Linear probe, validation Q8 / Q3 | 73.51 / 84.36 | 70.53 / 82.24 |
| Linear probe, CB513 Q8 / Q3 | **70.97 / 83.73** | 67.45 / 81.18 |
| Linear probe, CB513 G recall | 24.1% | 13.1% |
| Linear probe, CB513 G predicted as H | 28.9% | 32.4% |

The linear probe (`scripts/check_cache.py`: one `nn.Linear` on frozen features, 2 epochs) is a
sanity check, not the baseline. Misaligned features would score near the majority-class rate
(~31% on CB513). Even this probe sends ~29% of G residues to H. That is close to SecNet's 27.87%
and supports the project's motivation; Role 3 should confirm it with the BiLSTM.

Loader speed: `make_loader` with batch size 32 and `num_workers=0` delivers ~150,000–200,000
residues/s from the memmap, so data loading is unlikely to be the bottleneck.
On Windows, `num_workers > 0` needs the usual `if __name__ == "__main__":` guard.

**Sharing the cache.** `data/embeddings/*.npy`, `data/embeddings/*.json`, and `data/raw/` are gitignored. Upload the embedding files to the team Drive (`data/processed/` and `splits/` are already tracked in Git), or teammates can regenerate them with the commands below (needs a CUDA GPU for reasonable speed).

## 4. Reproduce from scratch

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/Scripts/python.exe "torch==2.14.1+cu130" --index-url https://download.pytorch.org/whl/cu130
uv pip install --python .venv/Scripts/python.exe -r requirements.txt

# 1. download into data/raw/ (see links on the NetSurfP-2.0 page)
#    https://services.healthtech.dtu.dk/services/NetSurfP-2.0/training_data/Train_HHblits.npz
#    https://services.healthtech.dtu.dk/services/NetSurfP-2.0/training_data/CB513_HHblits.npz
python scripts/prepare_data.py        # parse + audit       (~20 s)
python scripts/run_mmseqs.py          # MMseqs2 searches    (~1 min; needs MMseqs2 in tools/ or on PATH)
python scripts/make_splits.py         # splits.json + stats (~10 s)
python scripts/extract_esm.py         # ESM-2 650M cache    (see timing above)
python scripts/check_cache.py         # linear-probe sanity check of the cache
python -m pytest                      # unit tests
```

On Windows, MMseqs2 comes from the official build (`mmseqs-win64.zip` from mmseqs.com, unzipped
into `tools/`). `run_mmseqs.py` calls its modules directly, so no BusyBox install or admin rights
are needed. On Linux/macOS use `conda install -c bioconda mmseqs2`.
