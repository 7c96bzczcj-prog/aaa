#!/bin/bash
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
q="$1"; out="$2"
ids=$(curl -s --max-time 120 "$E/esearch.fcgi?db=gds&retmax=120&term=$(python3 -c "import urllib.parse,sys;print(urllib.parse.quote(sys.argv[1]))" "$q")" | grep -o '<Id>[0-9]*</Id>' | sed 's/<[^>]*>//g' | tr '\n' ',')
[ -z "$ids" ] && { echo "no hits: $q"; exit 0; }
curl -s --max-time 180 "$E/esummary.fcgi?db=gds&id=${ids%,}" > "$out"
echo "$(grep -c '<DocSum>' "$out") hits -> $out"
