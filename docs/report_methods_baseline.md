# Baseline method and handoff (Role 3)

## Method

We train a Q8 secondary-structure baseline on the cached, frozen final-layer
ESM-2 representations (1,280 features per residue). The classifier has one
bidirectional LSTM layer with 256 hidden units per direction, followed by 0.2
dropout and a linear layer with eight outputs in G, H, I, B, E, S, T, C order.
Sequences are packed by their true lengths before the LSTM, so batch padding
does not affect neighboring residue representations. Cross-entropy is computed
only where the shared data loader marks a residue as valid; masked and padded
positions have target value -100.

We use the five similarity-controlled folds defined in `splits/splits.json`. For
each fold, a model trains on the other four folds and predicts probabilities for
every residue of the held-out fold. The separate validation partition selects
the best checkpoint by residue-weighted cross-entropy, with patience of four
epochs and a maximum of 20 epochs. The resulting out-of-fold probabilities are
keyed by unique protein uid, never by PDB chain ID. A final model trains on all
five folds under the same validation procedure and predicts the validation and
CB513 partitions. CB513 is read only after that final checkpoint is fixed.

The optimizer is AdamW with learning rate 0.001 and weight decay 0.0001.
Training clips gradient norm at 1.0 and uses mixed precision on CUDA. The
default batch size is 16, with a prediction batch size of 32. Record any changes
to these settings when reporting results.

## Results

| Partition | Scored residues | Q8 (%) | Q3 (%) | G recall (%) | G predicted as H (%) |
|---|---:|---:|---:|---:|---:|
| OOF, folds 0–4 | 2,251,604 | 77.57 | 87.25 | 40.85 | 21.18 |
| Validation, final model | 245,405 | 77.63 | 86.94 | 45.51 | 19.13 |
| CB513, final model | 83,772 | 75.35 | 86.44 | 47.64 | 17.95 |

The five individual held-out Q8 scores range from 77.16% to 78.52%.
The final baseline still confuses 3-10 helices with alpha helices: 17.95% of
scored CB513 G residues are predicted as H. All metrics apply `valid_mask`,
which includes the CB513 evaluation and disorder masks. The linear-probe
figures in the main README are a separate cache check. Full fold results and
checkpoint epochs are in `outputs/baseline_summary.json`.

## Prediction handoff

Each `outputs/baseline_*.npy` file stores concatenated float32 Q8 softmax
probabilities, with one row per residue. Its `.index.json` file gives parallel
`uids`, `fold_ids`, and `offsets` arrays. For protein `uids[i]`, read rows
`offsets[i]:offsets[i + 1]`. All residues are present, including masked ones;
apply `valid_mask` when computing loss or metrics. The OOF index covers only
training folds 0–4, and each protein appears once. Checkpoints and output arrays
are gitignored, so share them with Roles 5 and 6 outside Git.
