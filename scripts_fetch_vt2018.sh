#!/usr/bin/env bash
# Fetch VT2018 (E-MTAB-6701) 10x processed matrix + metadata.
# Parallel byte-range download; verifies total size before concatenating.
set -euo pipefail

BASE="https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/701/E-MTAB-6701/Files"
DEST="${1:-data/VT2018}"
NPARTS="${NPARTS:-8}"

mkdir -p "$DEST/.parts"

for f in meta_10x.txt E-MTAB-6701.sdrf.txt E-MTAB-6701.idf.txt; do
  [ -s "$DEST/$f" ] || curl -sSL --retry 4 --retry-delay 2 -m 600 -o "$DEST/$f" "$BASE/$f"
done

TOTAL=$(curl -sSLI -m 120 "$BASE/raw_data_10x.txt" | awk 'BEGIN{IGNORECASE=1}/^content-length/{v=$2}END{gsub(/\r/,"",v);print v}')
echo "raw_data_10x.txt total bytes: $TOTAL"

if [ -s "$DEST/raw_data_10x.txt" ] && [ "$(stat -c%s "$DEST/raw_data_10x.txt")" = "$TOTAL" ]; then
  echo "already complete"; exit 0
fi

CHUNK=$(( (TOTAL + NPARTS - 1) / NPARTS ))
pids=()
for i in $(seq 0 $((NPARTS - 1))); do
  start=$((i * CHUNK))
  end=$((start + CHUNK - 1))
  [ "$end" -ge "$TOTAL" ] && end=$((TOTAL - 1))
  part="$DEST/.parts/part.$i"
  want=$((end - start + 1))
  if [ -s "$part" ] && [ "$(stat -c%s "$part")" = "$want" ]; then continue; fi
  ( curl -sSL --retry 6 --retry-delay 3 -m 5400 -r "${start}-${end}" -o "$part" "$BASE/raw_data_10x.txt" ) &
  pids+=($!)
done
for p in "${pids[@]:-}"; do [ -n "$p" ] && wait "$p"; done

for i in $(seq 0 $((NPARTS - 1))); do
  start=$((i * CHUNK)); end=$((start + CHUNK - 1))
  [ "$end" -ge "$TOTAL" ] && end=$((TOTAL - 1))
  want=$((end - start + 1)); got=$(stat -c%s "$DEST/.parts/part.$i")
  [ "$got" = "$want" ] || { echo "PART $i SIZE MISMATCH want=$want got=$got"; exit 1; }
done

cat "$DEST"/.parts/part.* > "$DEST/raw_data_10x.txt"
GOT=$(stat -c%s "$DEST/raw_data_10x.txt")
[ "$GOT" = "$TOTAL" ] || { echo "FINAL SIZE MISMATCH want=$TOTAL got=$GOT"; exit 1; }
rm -rf "$DEST/.parts"
echo "OK $GOT bytes"
