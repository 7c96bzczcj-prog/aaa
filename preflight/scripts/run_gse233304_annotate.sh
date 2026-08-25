#!/bin/bash
set -u
PY=/home/user/aaa/.venv/bin/python
cd /home/user/aaa/preflight
# wait for the HCA pass to finish so the two never contend for the 4 cores
while ! grep -q ALL_DONE scan/annotate_hca.log 2>/dev/null; do sleep 20; done
D=/home/user/aaa/data/gse233304
declare -A PAT=(
 [S1]=P1 [S3]=P1 [S4]=P2 [S6]=P2 [S7]=P7 [S8]=P8 [S9]=P9 [S11]=P9
 [S12]=c4 [S13]=c4 [S14]=c5 [S15]=c5 [S16]=c8 [S17]=c8 [S18]=c9 [S19]=c9
 [S20]=c10 [S21]=c10 [S22]=P4 [S26]=P15 [S30]=P15 [S32]=P16 [S36]=P16
 [S38]=P21 [S42]=P21 [S44]=P22 [S48]=P22 [S50]=P24 [S54]=P24 )
for p in P1 P2 P9 P15 P16 c4 c5 c8 c9 c10 P4 P7 P8 P21 P22 P24; do
  spec="["
  first=1
  for f in $D/*_barcodes.tsv.gz; do
    base=$(basename $f); gsm=${base%%_*}; sid=$(echo $base | cut -d_ -f2)
    [ "${PAT[$sid]:-}" = "$p" ] || continue
    b=$D/${gsm}_${sid}_barcodes.tsv.gz; ft=$D/${gsm}_${sid}_features.tsv.gz; mt=$D/${gsm}_${sid}_matrix.mtx.gz
    [ $first -eq 0 ] && spec="$spec,"
    spec="$spec[\"$sid\",\"$b|$ft|$mt\",\"mtx\"]"
    first=0
  done
  spec="$spec]"
  $PY scripts/annotate.py --donor $p --files "$spec" --out labels/GSE233304_$p.npz || echo "FAIL $p"
done
echo ALL_DONE
