# Role 5 handoff for evaluation

## What to give Role 6

Extract `role5_handoff.zip` from the team file share at the repository root. Role 6 can score the saved decisions immediately. The two files needed for CB513 are `outputs/role5/auditnet_cb513.npy` and `outputs/role5/auditnet_cb513.index.json`. For comparison, use Sujal’s `outputs/baseline_cb513.npy` and its index, and the ground truth and masks in `data/processed/cb513.npz`.

Read each array with `e2b.baseline.PredictionStore`, match chains by `uid`, take `argmax(axis=1)` for the Q8 labels in `GHIBESTC` order, and score only positions where `valid_mask` is true. The export includes masked residues so that row positions match the original chains. Role 4’s class-weighted predictions are a separate comparison.

The archive also contains the validation prediction pair, the fitted head in `checkpoints/hg_specialist.pt`, the frozen rule in `outputs/role5/routing.json`, the validation selection results, and metrics for each split. This package was reproduced locally with the original training code and seed. Its Q8, G F1 and H F1 match the saved Colab run at the precision reported there.

## How the specialist makes a decision

The eight-state baseline supplies one probability vector per residue. The specialist sees the frozen ESM-2 embedding at the central residue and five baseline vectors, one at each position from two residues before to two after. A 64-unit hidden layer estimates the probability of G given H or G. Training used true H and G examples among eligible positions, with out-of-fold baseline predictions from all five training folds. Validation selected the epoch and the rule; we then froze both before evaluating CB513.

For the selected rule, H and G must be the two highest baseline scores, with distinct values and a gap of at most 0.35. A baseline H changes to G when the specialist’s G probability is above 0.50. A baseline G changes to H when that probability is below 0.30. In all other cases, the baseline decision stands. On a change, the export swaps just the two baseline scores. Use the exported values for decisions. They are not calibrated estimates of the specialist’s confidence.

## Results to check

| Partition | Scored residues | Baseline G F1 | With specialist | Errors fixed | New errors | Wrong before and after |
|---|---:|---:|---:|---:|---:|---:|
| Validation | 245,405 | 48.95% | 49.12% | 95 | 92 | 86 |
| CB513 | 83,772 | 47.00% | 47.12% | 28 | 33 | 23 |

On CB513, 84 scored residues changed label. The specialist identified 15 more true G residues and 20 fewer true H residues, for five fewer correct Q8 predictions overall. Q8 accuracy went from 75.355% to 75.349%. The 23 changes that stayed wrong had a true label outside H/G. Q3 decisions are identical because H and G both map to helix. For the poster, the supported conclusion is a small gain in G detection with a cost in H accuracy.

## Reproducing the exports

With Leo’s four 650M embedding files and Sujal’s three baseline prediction pairs in the locations specified in the README, run the following at the repository root.

```sh
python scripts/train_specialist.py
python scripts/train_specialist.py --export-cb513
```

To regenerate only the CB513 decisions with the packaged weights, extract the archive and run the second command. It reads the saved checkpoint and its frozen rule.
