"""Shared constants for the E2B-AuditNet data pipeline.

Channel layout of the NetSurfP-2.0 ``.npz`` files (``data[samples, length, 68]``),
as documented on the NetSurfP-2.0 "Training/Test Data" tab.
"""

# One-hot order of channels [0:20]. The NetSurfP page does not state it; we
# verified it by decoding CB513 entries 154l-A and 1a45-A and comparing them
# with the RCSB FASTA records (exact match). Unknown residues are all-zero -> X.
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
UNKNOWN_AA = "X"

# Q8 labels, channels [57:65], in this order. 0 = G (3-10 helix), 1 = H (alpha helix).
Q8_ORDER = "GHIBESTC"
Q8_INDEX = {c: i for i, c in enumerate(Q8_ORDER)}
G_IDX = Q8_INDEX["G"]
H_IDX = Q8_INDEX["H"]

# NetSurfP-2.0: "Q8 -> Q3: HHHEECCC"
Q3_ORDER = "HEC"
Q8_TO_Q3 = [0, 0, 0, 1, 1, 2, 2, 2]

# Channel indices
CH_AA = slice(0, 20)
CH_HMM = slice(20, 50)  # not used: replaced by frozen ESM-2
CH_SEQ_MASK = 50  # 1 = residue, 0 = padding
CH_DISORDER = 51  # 1 = ordered, 0 = disordered
CH_EVAL = 52  # CB513 only: 1 = evaluate, 0 = ignore
CH_Q8 = slice(57, 65)

# fold_id values from the shared data contract
FOLD_VALIDATION = -1
FOLD_CB513 = -2

# Label value stored for residues that have no Q8 label (never valid).
LABEL_MISSING = -1
# Target value used in batched tensors for positions that must not be scored.
IGNORE_INDEX = -100

# ESM-2 positional limit: 1024 tokens including BOS/EOS.
ESM_MAX_RESIDUES = 1022
