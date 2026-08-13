#!/usr/bin/env bash
# Fetch VT2018 (E-MTAB-6701) 10x processed matrix + metadata.
#
# Byte-range download with per-part size verification.
#
# Header parsing here deliberately avoids awk's IGNORECASE: the system awk is
# mawk, which accepts `BEGIN{IGNORECASE=1}` and then ignores it, so a
# case-insensitive header match silently yields an empty size -- which turns
# every computed byte range into nonsense while curl still exits 0. Per-part
# size verification is what catches that class of failure.
set -euo pipefail

BASE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/701/E-MTAB-6701/Files"
DEST="${1:-data/VT2018}"
NPARTS="${NPARTS:-4}"
TRIES="${TRIES:-5}"

mkdir -p "$DEST/.parts"

for f in meta_10x.txt E-MTAB-6701.sdrf.txt E-MTAB-6701.idf.txt; do
  [ -s "$DEST/$f" ] || curl -sSL --retry 4 --retry-delay 2 -m 600 -o "$DEST/$f" "$BASE/$f"
done

# Total size. HEAD against this endpoint is flaky (it intermittently answers
# with no Content-Length), and an empty TOTAL silently turns every subsequent
# range into nonsense -- so retry, then fall back to the Content-Range of a
# 1-byte GET, and refuse to proceed without a plausible number.
get_total() {
  local v=""
  for _ in 1 2 3; do
    v=$(curl -sSLI -m 120 "$BASE/raw_data_10x.txt" \
        | grep -i '^content-length:' | tail -1 | tr -dc '0-9')
    [[ "$v" =~ ^[0-9]+$ ]] && { echo "$v"; return 0; }
    sleep 3
  done
  v=$(curl -sS -L -m 120 -r 0-0 -D - -o /dev/null "$BASE/raw_data_10x.txt" \
      | grep -i '^content-range:' | tail -1 | sed 's|.*/||' | tr -dc '0-9')
  [[ "$v" =~ ^[0-9]+$ ]] && { echo "$v"; return 0; }
  return 1
}

TOTAL=$(get_total) || { echo "could not determine remote size; aborting"; exit 1; }
if ! [[ "$TOTAL" =~ ^[0-9]+$ ]] || [ "$TOTAL" -lt 1000000 ]; then
  echo "implausible remote size: '$TOTAL'"; exit 1
fi
echo "raw_data_10x.txt total bytes: $TOTAL"

if [ -s "$DEST/raw_data_10x.txt" ] && [ "$(stat -c%s "$DEST/raw_data_10x.txt")" = "$TOTAL" ]; then
  echo "already complete"; exit 0
fi

CHUNK=$(( (TOTAL + NPARTS - 1) / NPARTS ))

fetch_part() {           # $1 = part index
  local i="$1"
  local start=$((i * CHUNK))
  local end=$((start + CHUNK - 1))
  [ "$end" -ge "$TOTAL" ] && end=$((TOTAL - 1))
  local want=$((end - start + 1))
  local part="$DEST/.parts/part.$i"

  for try in $(seq 1 "$TRIES"); do
    if [ -s "$part" ] && [ "$(stat -c%s "$part")" = "$want" ]; then
      return 0
    fi
    curl -sS -L --fail --retry 3 --retry-delay 3 -m 7200 \
         -H 'Accept-Encoding: identity' \
         -r "${start}-${end}" -o "$part" "$BASE/raw_data_10x.txt" || true
    local got=0; [ -f "$part" ] && got=$(stat -c%s "$part")
    if [ "$got" = "$want" ]; then return 0; fi
    echo "part $i attempt $try: want=$want got=$got, retrying" >&2
    rm -f "$part"
    sleep $((try * 3))
  done
  echo "part $i FAILED after $TRIES attempts" >&2
  return 1
}

pids=()
for i in $(seq 0 $((NPARTS - 1))); do
  fetch_part "$i" &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ "$fail" = 0 ] || { echo "one or more parts failed"; exit 1; }

cat "$DEST"/.parts/part.* > "$DEST/raw_data_10x.txt"
GOT=$(stat -c%s "$DEST/raw_data_10x.txt")
[ "$GOT" = "$TOTAL" ] || { echo "FINAL SIZE MISMATCH want=$TOTAL got=$GOT"; exit 1; }
rm -rf "$DEST/.parts"
echo "OK $GOT bytes"
