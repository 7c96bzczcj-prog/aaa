#!/usr/bin/env python3
"""Chunked, resumable, offset-verified download of the Zenodo NK atlas objects.

The proxy drops long single-stream transfers, so the file is pulled in bounded byte
ranges. Each response's Content-Range is checked and the payload is written at the
exact offset (never appended blindly), so a server that ignores Range can't corrupt
or inflate the file. Progress is tracked in a sidecar so restarts resume cleanly.
"""
import json, os, sys, time
import requests

URL = "https://zenodo.org/api/records/8275845/files/%s/content"
CHUNK = 32 * 1024 * 1024


def download(fname, dest, expected):
    url = URL % fname
    state_path = dest + ".state"
    done = 0
    if os.path.exists(state_path) and os.path.exists(dest):
        try:
            done = int(json.load(open(state_path))["done"])
        except Exception:  # noqa: BLE001
            done = 0
    # preallocate
    if not os.path.exists(dest) or os.path.getsize(dest) != expected:
        with open(dest, "wb") as f:
            f.truncate(expected)
        if done and os.path.getsize(dest) != expected:
            done = 0
    fails = 0
    while done < expected:
        end = min(done + CHUNK, expected) - 1
        hdr = {"Range": "bytes=%d-%d" % (done, end), "User-Agent": "npc-sc-audit/1.0"}
        try:
            r = requests.get(url, headers=hdr, stream=True, timeout=300)
            if r.status_code != 206:
                raise RuntimeError("no partial content: HTTP %s" % r.status_code)
            cr = r.headers.get("Content-Range", "")
            # expect: bytes <start>-<end>/<total>
            start = int(cr.split()[1].split("-")[0])
            if start != done:
                raise RuntimeError("server returned wrong offset %d (wanted %d)" % (start, done))
            want = end - done + 1
            buf = bytearray()
            for blk in r.iter_content(1024 * 512):
                if blk:
                    buf.extend(blk)
                    if len(buf) > want:
                        break
            if len(buf) < want:
                raise RuntimeError("short chunk %d/%d" % (len(buf), want))
            with open(dest, "r+b") as f:
                f.seek(done)
                f.write(bytes(buf[:want]))
            done += want
            json.dump({"done": done}, open(state_path, "w"))
            sys.stderr.write("  %d/%d (%.1f%%)\n" % (done, expected, 100.0 * done / expected))
            fails = 0
        except Exception as e:  # noqa: BLE001
            fails += 1
            sys.stderr.write("  offset %d failed (%s): %s\n" % (done, type(e).__name__, str(e)[:130]))
            if fails >= 8:
                raise SystemExit("stalled at offset %d" % done)
            time.sleep(min(2 ** fails, 30))
    size = os.path.getsize(dest)
    ok = size == expected
    print("%s -> %s (%d bytes, expected %d, match=%s)" % (fname, dest, size, expected, ok))
    return ok


if __name__ == "__main__":
    fname, dest, expected = sys.argv[1], sys.argv[2], int(sys.argv[3])
    sys.exit(0 if download(fname, dest, expected) else 1)
