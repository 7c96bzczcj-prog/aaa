#!/bin/bash
cd /home/user/aaa/data/raw
B=/home/user/aaa/phase1_registry/batches_to_fetch.txt
fetch(){ f="GSE154826_amp_batch_ID_$1.tar.gz"
  [ -s "$f" ] && return 0
  for a in 1 2 3 4; do
    curl -sS --max-time 900 -o "$f" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154826/suppl/$f" && [ -s "$f" ] && return 0
    sleep $((2**a))
  done
  echo "FAILED $1" >> /home/user/aaa/data/raw/failed.txt; }
export -f fetch
cat $B | xargs -P 4 -I{} bash -c 'fetch {}'
echo "DONE downloads: $(ls -1 *.tar.gz 2>/dev/null|wc -l)"
