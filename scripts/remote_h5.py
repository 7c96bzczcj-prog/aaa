"""Read-only, block-cached HTTP range file object so h5py can open a remote
.h5ad and read only the groups it touches (e.g. obs), without downloading X."""
from __future__ import annotations
import io, requests

class RangeFile(io.RawIOBase):
    def __init__(self, url: str, block: int = 4 << 20):
        self.url, self.block, self.pos, self.cache = url, block, 0, {}
        self.s = requests.Session()
        r = self.s.get(url, headers={"Range": "bytes=0-0"}, timeout=60, allow_redirects=True)
        assert r.status_code == 206, r.status_code
        self.size = int(r.headers["Content-Range"].split("/")[1])
        self.fetched = 0
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos
    def _blk(self, i):
        if i not in self.cache:
            s = i * self.block; e = min(s + self.block, self.size) - 1
            for a in range(5):
                try:
                    r = self.s.get(self.url, headers={"Range": f"bytes={s}-{e}"}, timeout=120)
                    if r.status_code == 206 and len(r.content) == e - s + 1: break
                except requests.RequestException: pass
            else: raise IOError(f"range {s}-{e} failed")
            self.cache[i] = r.content; self.fetched += len(r.content)
        return self.cache[i]
    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0: return 0
        out, p, end = bytearray(), self.pos, self.pos + n
        while p < end:
            i = p // self.block; d = self._blk(i); o = p - i * self.block
            take = min(len(d) - o, end - p); out += d[o:o + take]; p += take
        b[:n] = out; self.pos += n; return n
