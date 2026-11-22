"""Parser for the NetSurfP-2.0 ``*_HHblits.npz`` archives.

The Train archive holds a float64 array of shape (10848, 1632, 68), about 9.6 GB
once decompressed. ``iter_samples`` streams it one protein at a time from the zip
member, so parsing needs well under 1 GB of RAM.

Only channels [0:20], [50], [51], [52] and [57:65] are read (see constants.py).
"""

from __future__ import annotations

import hashlib
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np

from .constants import (
    AA_ORDER,
    CH_AA,
    CH_DISORDER,
    CH_EVAL,
    CH_Q8,
    CH_SEQ_MASK,
    LABEL_MISSING,
    UNKNOWN_AA,
)

_AA_ARRAY = np.array(list(AA_ORDER))


@dataclass
class ProteinRecord:
    uid: str  # unique key, e.g. "train_00042" / "cb513_00007"
    pdb_id: str  # identifier from `pdbids` (not unique in CB513: 513 regions, 434 chains)
    source: str  # "train" or "cb513"
    npz_index: int
    seq: str
    q8_label: np.ndarray  # int8 [L]; 0..7 (GHIBESTC) or LABEL_MISSING
    disorder_mask: np.ndarray  # bool [L]; True = ordered
    eval_mask: np.ndarray  # bool [L]; True = evaluate (all True outside CB513)
    valid_mask: np.ndarray  # bool [L]; residue may be used in loss/metrics

    @property
    def length(self) -> int:
        return len(self.seq)


def sha256sum(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def _read_header(f):
    version = np.lib.format.read_magic(f)
    if version == (1, 0):
        return np.lib.format.read_array_header_1_0(f)
    if version == (2, 0):
        return np.lib.format.read_array_header_2_0(f)
    raise ValueError(f"unsupported .npy version {version}")


def read_pdbids(path: str | Path) -> np.ndarray:
    with np.load(path, allow_pickle=False) as z:
        return z["pdbids"]


def iter_samples(path: str | Path) -> Iterator[tuple[int, np.ndarray]]:
    """Yield (index, array[length, 68]) without loading the whole archive."""
    with zipfile.ZipFile(path) as zf, zf.open("data.npy") as f:
        shape, fortran, dtype = _read_header(f)
        if fortran or len(shape) != 3:
            raise ValueError(f"unexpected data layout: shape={shape}, fortran={fortran}")
        n, max_len, n_ch = shape
        per_sample = max_len * n_ch * dtype.itemsize
        for i in range(n):
            buf = f.read(per_sample)
            if len(buf) != per_sample:
                raise IOError(f"truncated archive at sample {i}")
            yield i, np.frombuffer(buf, dtype=dtype).reshape(max_len, n_ch)


def decode_sequence(onehot: np.ndarray) -> str:
    """[L, 20] one-hot -> amino-acid string; all-zero rows become X."""
    present = onehot.sum(axis=1) > 0
    letters = _AA_ARRAY[onehot.argmax(axis=1)]
    letters[~present] = UNKNOWN_AA
    return "".join(letters.tolist())


def parse_sample(
    sample: np.ndarray, pdb_id: str, source: str, index: int, use_eval_mask: bool
) -> ProteinRecord:
    seq_mask = sample[:, CH_SEQ_MASK] > 0.5
    length = int(seq_mask.sum())
    if not seq_mask[:length].all():
        raise ValueError(f"{source}[{index}] {pdb_id}: sequence mask is not a contiguous prefix")
    s = sample[:length]

    onehot = s[:, CH_AA]
    if (onehot.sum(axis=1) > 1.5).any():
        raise ValueError(f"{source}[{index}] {pdb_id}: residue with more than one amino-acid bit")
    seq = decode_sequence(onehot)

    q8 = s[:, CH_Q8]
    has_label = q8.sum(axis=1) > 0.5
    if (q8.sum(axis=1) > 1.5).any():
        raise ValueError(f"{source}[{index}] {pdb_id}: residue with more than one Q8 bit")
    q8_label = np.where(has_label, q8.argmax(axis=1), LABEL_MISSING).astype(np.int8)

    disorder_mask = s[:, CH_DISORDER] > 0.5
    eval_mask = s[:, CH_EVAL] > 0.5 if use_eval_mask else np.ones(length, dtype=bool)
    valid_mask = disorder_mask & eval_mask & has_label

    return ProteinRecord(
        uid=f"{source}_{index:05d}",
        pdb_id=str(pdb_id),
        source=source,
        npz_index=index,
        seq=seq,
        q8_label=q8_label,
        disorder_mask=disorder_mask,
        eval_mask=eval_mask,
        valid_mask=valid_mask,
    )


def parse_netsurfp(path: str | Path, source: str, use_eval_mask: bool | None = None) -> list[ProteinRecord]:
    """Parse one NetSurfP-2.0 archive. The evaluation mask is applied only to CB513 by default."""
    if use_eval_mask is None:
        use_eval_mask = source == "cb513"
    ids = read_pdbids(path)
    return [parse_sample(arr, ids[i], source, i, use_eval_mask) for i, arr in iter_samples(path)]


# ---------------------------------------------------------------------------
# Compact on-disk format: one .npz per source with residues concatenated.
# No pickled objects, so teammates can load it with allow_pickle=False.
# ---------------------------------------------------------------------------


def save_records(records: list[ProteinRecord], path: str | Path) -> None:
    lengths = np.array([r.length for r in records], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    np.savez_compressed(
        path,
        uid=np.array([r.uid for r in records]),
        pdb_id=np.array([r.pdb_id for r in records]),
        source=np.array([r.source for r in records]),
        npz_index=np.array([r.npz_index for r in records], dtype=np.int64),
        offsets=offsets,
        seq_concat=np.frombuffer("".join(r.seq for r in records).encode("ascii"), dtype=np.uint8),
        q8_label=np.concatenate([r.q8_label for r in records]),
        disorder_mask=np.concatenate([r.disorder_mask for r in records]),
        eval_mask=np.concatenate([r.eval_mask for r in records]),
        valid_mask=np.concatenate([r.valid_mask for r in records]),
    )


def load_records(path: str | Path) -> list[ProteinRecord]:
    with np.load(path, allow_pickle=False) as z:
        off = z["offsets"]
        seq_all = z["seq_concat"].tobytes().decode("ascii")
        q8, dis, ev, va = z["q8_label"], z["disorder_mask"], z["eval_mask"], z["valid_mask"]
        out = []
        for i, uid in enumerate(z["uid"]):
            a, b = int(off[i]), int(off[i + 1])
            out.append(
                ProteinRecord(
                    uid=str(uid),
                    pdb_id=str(z["pdb_id"][i]),
                    source=str(z["source"][i]),
                    npz_index=int(z["npz_index"][i]),
                    seq=seq_all[a:b],
                    q8_label=q8[a:b],
                    disorder_mask=dis[a:b],
                    eval_mask=ev[a:b],
                    valid_mask=va[a:b],
                )
            )
    return out


def write_fasta(records: list[ProteinRecord], path: str | Path) -> None:
    with open(path, "w", encoding="ascii", newline="\n") as f:
        for r in records:
            f.write(f">{r.uid} {r.pdb_id}\n{r.seq}\n")
