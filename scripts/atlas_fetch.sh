#!/bin/bash
# Parallel ranged download of Zenodo files with md5 verification against the Zenodo API.
# usage: atlas_fetch.sh <record_id> <out_dir> <file_key> [file_key ...]
set -u
REC=$1; OUT=$2; shift 2
mkdir -p "$OUT"; cd "$OUT"
API=$(curl -sS -L --max-time 60 -H "Accept: application/json" "https://zenodo.org/api/records/$REC")
for KEY in "$@"; do
  read SIZE MD5 < <(python3 -c "import json,sys;d=json.loads(sys.argv[1]);f=[x for x in d['files'] if x['key']==sys.argv[2]][0];print(f['size'],f['checksum'].split(':')[1])" "$API" "$KEY")
  if [ -s "$KEY" ] && [ "$(md5sum "$KEY" | cut -d' ' -f1)" = "$MD5" ]; then echo "OK (cached) $KEY"; continue; fi
  URL="https://zenodo.org/records/$REC/files/$(python3 -c "import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))" "$KEY")?download=1"
  CH=25000000; N=$(( (SIZE + CH - 1) / CH )); P=".parts_${KEY// /_}"; mkdir -p "$P"
  seq 0 $((N-1)) | xargs -P 12 -I{} bash -c '
    i={}; s=$((i*'$CH')); e=$((s+'$CH'-1)); [ $e -ge '$SIZE' ] && e=$(('$SIZE'-1)); want=$((e-s+1)); f="'$P'/$(printf %06d $i)"
    [ -s "$f" ] && [ $(stat -c %s "$f") -eq $want ] && exit 0
    for a in 1 2 3 4 5; do curl -sS -L --max-time 300 -r $s-$e -o "$f" "'"$URL"'" && [ $(stat -c %s "$f") -eq $want ] && exit 0; sleep $((2**a)); done
    echo "CHUNK FAIL $i" >&2; exit 1'
  cat "$P"/* > "$KEY"
  GOT=$(md5sum "$KEY" | cut -d' ' -f1)
  if [ "$GOT" = "$MD5" ]; then echo "OK $KEY $SIZE md5=$MD5"; rm -rf "$P"; else echo "MD5 MISMATCH $KEY got=$GOT want=$MD5"; fi
done
