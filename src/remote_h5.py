"""Read the metadata of a remote .h5ad over HTTP range requests.

Why this exists: an h5ad's `obs`, `var` and the *structure* of `X` are a few MB
inside a file that may be tens of GB.  Deciding whether a dataset can answer a
question does not require the matrix, and downloading it to find out is the
expensive mistake this module avoids.

curl is used rather than a Python HTTP client because the container's outbound
path is a configured proxy with its own CA bundle; curl already has it.

No scanpy (R10).  h5py + numpy only.
"""

from __future__ import annotations

import io
import subprocess
import time

import h5py
import numpy as np


class HTTPRangeFile(io.RawIOBase):
    """Minimal seekable read-only file object backed by HTTP range requests."""

    def __init__(self, url: str, block: int = 1 << 20, max_cache: int = 8):
        self.url = url
        self.block = block
        self.max_cache = max_cache
        self._pos = 0
        self._cache: dict[int, bytes] = {}
        self.bytes_fetched = 0
        self.size = self._content_length()

    # -- transport ---------------------------------------------------------
    def _content_length(self) -> int:
        out = subprocess.run(
            ["curl", "-sSL", "-o", "/dev/null", "-w", "%{size_download}",
             "-r", "0-0", self.url],
            capture_output=True, text=True, check=True,
        )
        # a 1-byte range reply tells us little; ask for the header instead
        hdr = subprocess.run(
            ["curl", "-sSLI", self.url], capture_output=True, text=True, check=True,
        ).stdout
        length = None
        for line in hdr.splitlines():
            # mawk-style case bugs are why this is an explicit lower() compare
            if line.lower().startswith("content-length:"):
                length = int(line.split(":", 1)[1].strip())
        if length is None:
            raise RuntimeError(f"no Content-Length for {self.url}\n{hdr}\n{out.stdout}")
        return length

    def _fetch(self, index: int) -> bytes:
        if index in self._cache:
            return self._cache[index]
        start = index * self.block
        end = min(start + self.block, self.size) - 1
        want = end - start + 1
        data = b""
        last = ""
        for attempt in range(5):
            proc = subprocess.run(
                ["curl", "-sSL", "--retry", "2", "-r", f"{start}-{end}", self.url],
                capture_output=True,
            )
            data = proc.stdout
            # the failure mode that cost a day: a range reply that is silently
            # the whole file, or empty, with exit status 0.  Every part is
            # size-verified, never trusted on exit status.
            if proc.returncode == 0 and len(data) == want:
                break
            last = f"rc={proc.returncode} got={len(data)} want={want}"
            time.sleep(2 ** attempt)
        else:
            raise RuntimeError(f"range {start}-{end} failed after 5 tries: {last}")
        self.bytes_fetched += len(data)
        # Access is sequential but alternates between two dataset regions
        # (indices and data), so evict FIFO by insertion order rather than by
        # block index — dropping the lowest index would evict the region that
        # is still being walked.
        while len(self._cache) >= self.max_cache:
            del self._cache[next(iter(self._cache))]
        self._cache[index] = data
        return data

    # -- file protocol -----------------------------------------------------
    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            self._pos = offset
        elif whence == io.SEEK_CUR:
            self._pos += offset
        elif whence == io.SEEK_END:
            self._pos = self.size + offset
        return self._pos

    def tell(self) -> int:
        return self._pos

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self.size - self._pos
        n = max(0, min(n, self.size - self._pos))
        out = bytearray()
        pos = self._pos
        while len(out) < n:
            idx = pos // self.block
            off = pos - idx * self.block
            chunk = self._fetch(idx)[off:]
            take = min(len(chunk), n - len(out))
            out += chunk[:take]
            pos += take
        self._pos = pos
        return bytes(out)

    def readinto(self, b) -> int:
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def open_remote(url: str, block: int = 1 << 20, max_cache: int = 8):
    """Open a remote HDF5 file. Returns (h5py.File, HTTPRangeFile)."""
    fobj = HTTPRangeFile(url, block=block, max_cache=max_cache)
    return h5py.File(fobj, "r"), fobj


def read_str_array(group) -> np.ndarray:
    """Read an anndata string array (categorical or plain) as python strings."""
    if isinstance(group, h5py.Group) and "categories" in group:
        cats = group["categories"][:]
        codes = group["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
        return out
    arr = group[:]
    return np.array([a.decode() if isinstance(a, bytes) else str(a) for a in arr])
