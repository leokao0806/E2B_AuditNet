# Draft text for the final report — Method: Data and Features (Roles 1 and 2)

Drop-in paragraphs for the "Method" section (IEEE template). Numbers come from
`splits/split_stats.json`, `data/processed/data_audit.json` and `data/embeddings/*.index.json`;
re-check them if the pipeline is re-run.

---

## A. Dataset

We use the NetSurfP-2.0 release [b6]: `Train_HHblits.npz` (10,848 PDB chains) for training and
validation, and `CB513_HHblits.npz` (513 regions from 434 chains) as the held-out test set. Each
file stores, per residue, a one-hot amino-acid encoding, a sequence mask, a disorder mask, an
evaluation mask (CB513 only) and the eight-state DSSP label (G, H, I, B, E, S, T, C) [b1]. We
discard the HHblits profile channels because frozen ESM-2 representations replace them. The
one-hot order (ACDEFGHIKLMNPQRSTVWY) is not documented by the release; we confirmed it by matching
decoded CB513 chains against their RCSB sequences. A residue is scored only if it is ordered and
labelled and, on CB513, inside the evaluation mask. This leaves 2,598,847 scored residues in the
training file and 83,772 in CB513. G makes up 3.87% and H 34.31% of the scored training residues
(3.70% and 30.69% in CB513).

## B. Leakage control and partitioning

The release states that no training chain shares more than 25% identity with CB513. Because
SecNet [b2] showed that such overlap inflates scores, we re-checked it with exact-sequence
matching and an MMseqs2 search (sensitivity 7.5). Two CB513 chains were identical to training
chains and one was contained in a training chain. With local alignments (E-value ≤ 10⁻³, at least
20 aligned residues), 203 of the 513 CB513 regions matched 323 training chains at ≥ 25% identity,
and 13 matched training chains at ≥ 90% identity. We removed every training chain that matched
CB513 under this rule (324 chains, 3.0% of the file). This is stricter than E2BNet [b5], which
removed exact duplicates only.

The remaining 10,524 chains are grouped so that homologous chains never fall in different
partitions. Two chains are linked if MMseqs2 finds an alignment with E-value ≤ 10⁻³, at least 30%
identity, and coverage of at least half of the shorter chain; groups are the connected components
of these links (9,376 groups, the largest with 713 chains). The coverage condition is necessary:
without it, short identical expression tags (e.g. a 20-residue His-tag/thrombin site) linked
thousands of unrelated chains. We randomly assign whole groups to a validation partition (1,053
chains, 245,405 scored residues), used only for model selection, class-weight tuning and routing
thresholds. The other groups are distributed over five cross-validation folds with equal numbers
of scored residues (450,320–450,321 per fold; G 3.78–3.91%). These folds produce the out-of-fold
baseline predictions used to train the H/G specialist.

## C. Frozen ESM-2 features

Each chain is encoded with the pretrained ESM-2 model esm2\_t33\_650M\_UR50D (650M parameters) [b3, b4], with all weights
frozen. We take the final-layer representation of every residue (dimension 1,280) and cache
it on disk in half precision, so that all models read identical features. E2BNet truncated chains
to 128 residues [b5]; we instead encode full chains up to ESM-2's 1,022-residue limit. The 36
longer chains (up to 1,632 residues) are split into 1,022-residue windows that overlap by 256
residues, and each residue takes its representation from the window in which it has the most
sequence context on both sides. Half-precision inference matched full precision closely (minimum
per-residue cosine similarity 0.9997 on sampled chains). Extraction of all 11,361 chains took
5.5 minutes on one RTX 5060 Laptop GPU (8 GB), and the cache occupies 7.5 GB.

As a check that features and labels are aligned, a single linear layer trained on the frozen
features for two epochs reaches 70.97% Q8 and 83.73% Q3 on CB513. Even this simple probe
labels 28.9% of true G residues as H, consistent with the 27.87% reported by SecNet [b2].

## Table (optional): partitions

| Partition | Chains | Scored residues | G (%) | H (%) |
|---|---:|---:|---:|---:|
| CV folds 0–4 (total) | 9,471 | 2,251,604 | 3.84 | 34.36 |
| Validation | 1,053 | 245,405 | 3.86 | 34.88 |
| CB513 (test) | 513 | 83,772 | 3.70 | 30.69 |
| Removed (CB513 homologs) | 324 | 101,838 | 4.57 | 31.72 |
