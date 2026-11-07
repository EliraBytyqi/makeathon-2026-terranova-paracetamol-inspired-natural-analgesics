"""Transparently gzipped model serialisation.

Trained bundles carry two heavy things: the fitted ensemble, and an
applicability-domain fingerprint matrix with one 2048-bit row per training
compound. Both are highly repetitive, so the pickles compress about fourfold
- the hepatic panel goes from 128 MB to 32 MB.

That matters for distribution rather than for speed. GitHub refuses any file
over 100 MB, and the uncompressed panel exceeded it, which would have forced
either Git LFS (extra tooling on the reviewer's machine, and a bandwidth quota
that can expire mid-review) or shipping the repository without its models, so
a fresh clone could not run the live demo at all. Compressing instead keeps
every model in the repository and needs nothing installed.

Files keep the `.pkl` extension: the contents are gzipped, not the name. Every
existing glob, path constant and pipeline stage therefore keeps working, and
`load` sniffs the gzip magic number so uncompressed bundles written by an
earlier version still load.

Decompression costs well under a second per model, against a load time already
dominated by RDKit import.
"""

from __future__ import annotations

import gzip
import pickle
from pathlib import Path
from typing import Any

GZIP_MAGIC = b"\x1f\x8b"
COMPRESS_LEVEL = 6          # 4.0x on our bundles; level 9 buys ~2% for 3x the time


def save(obj: Any, path: str | Path) -> Path:
    """Pickle `obj` to `path`, gzipped."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb", compresslevel=COMPRESS_LEVEL) as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def load(path: str | Path) -> Any:
    """Unpickle from `path`, whether or not it is gzipped."""
    path = Path(path)
    with open(path, "rb") as fh:
        compressed = fh.read(2) == GZIP_MAGIC
    opener = gzip.open if compressed else open
    with opener(path, "rb") as fh:
        return pickle.load(fh)


def compress_in_place(path: str | Path) -> tuple[int, int]:
    """Rewrite an uncompressed bundle as a gzipped one. Returns (before, after).

    Writes to a temporary file first and replaces only on success, so an
    interrupted run cannot leave a truncated model behind.
    """
    path = Path(path)
    before = path.stat().st_size
    with open(path, "rb") as fh:
        if fh.read(2) == GZIP_MAGIC:
            return before, before
    obj = load(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    save(obj, tmp)
    tmp.replace(path)
    return before, path.stat().st_size


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1] / "models"
    total_before = total_after = 0
    for p in sorted(root.glob("*.pkl")):
        b, a = compress_in_place(p)
        total_before += b
        total_after += a
        note = "already compressed" if a == b else f"{b / a:.1f}x"
        print(f"  {p.name:34s} {b / 1e6:7.1f} MB -> {a / 1e6:6.1f} MB   {note}")
    print(f"\n  total {total_before / 1e6:.0f} MB -> {total_after / 1e6:.0f} MB")
