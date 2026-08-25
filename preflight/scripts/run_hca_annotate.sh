#!/bin/bash
set -u
PY=/home/user/aaa/.venv/bin/python
cd /home/user/aaa/preflight
for d in MantonBM1 MantonBM2 MantonBM3 MantonBM4 MantonBM5 MantonBM6 MantonBM7 MantonBM8 \
         MantonBL1 MantonBL2 MantonBL3 MantonBL4 MantonBL5 MantonBL6 MantonBL7 MantonBL8; do
  if [ -s labels/$d.npz ] && [ -s labels/${d}_clusters.json ] && [ "$d" != "MantonBM1" ] && [ "$d" != "MantonBL1" ]; then :; fi
  $PY -c "
import json,glob,os
fs=sorted(glob.glob('/home/user/aaa/data/hca_ica/${d}_*.h5'))
json.dump([[os.path.basename(f).replace('_raw_feature_bc_matrix.h5',''),f,'h5'] for f in fs],open('/tmp/${d}.json','w'))"
  $PY scripts/annotate.py --donor $d --files "$(cat /tmp/$d.json)" --out labels/$d.npz || echo "FAIL $d"
done
echo ALL_DONE
