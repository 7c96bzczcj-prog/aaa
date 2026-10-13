"""Cell types the screen's WITHIN-STRATUM relative rule excluded, but that the reference layer supports.

The CASP1 bottom-tercile test is relative: a cell type with genuinely low CASP1 still fails if a third of
its stratum is lower. This finds IL18-top cell types that were excluded on CASP1 rank alone yet show a
favourable IL18:CASP1 ratio in matched sorted/atlas references. These are POTENTIAL candidates, reported
alongside the screen's own output rather than instead of it.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
from common import write_tsv
from evidence_integrate import ratios_for, norm_source, RATIO_STRONG

ROOT = Path(__file__).resolve().parents[1]
ref = pd.read_csv(ROOT/"reference"/"ext"/"unified_reference.tsv", sep="\t")
ref = ref[ref.match_quality.isin(["exact","close"])]

# explicit screen-label -> reference-label map for non-candidate cell types worth checking
EXTRA = {
 ("human","lung","pulmonary alveolar type 1 cell"):"lung AT1",
 ("human","lung","pulmonary alveolar type 2 cell"):"lung AT2",
 ("human","lung","alveolar macrophage"):"lung alveolar macrophage",
 ("human","lung","respiratory basal cell"):"lung multiciliated/ciliated epithelial cell",
 ("human","liver","hepatocyte"):"liver hepatocyte",
 ("human","liver","cholangiocyte"):"liver cholangiocyte",
 ("human","pancreas","pancreatic acinar cell"):"pancreas acinar cell",
 ("human","pancreas","pancreatic ductal cell"):"pancreas ductal cell",
 ("human","pancreas","type B pancreatic cell"):"islet beta cell",
 ("mouse","liver","hepatocyte"):"liver hepatocyte",
 ("mouse","liver","Kupffer cell"):"liver Kupffer cell",
 ("mouse","pancreas","pancreatic acinar cell"):"pancreas acinar cell",
 ("mouse","pancreas","type B pancreatic cell"):"islet beta cell",
 ("mouse","pancreas","pancreatic A cell"):"islet alpha cell",
 ("mouse","brain","neuron"):"brain neuron",
 ("mouse","lung","lung macrophage"):"lung interstitial macrophage",
}
s = pd.read_csv(ROOT/"results"/"il18_casp_summary.tsv", sep="\t", comment="#")
miss = s[(s.IL18_top == True) & (s.is_candidate != True) & (s.n_donors >= 3) & (s.level == "author")]
rows=[]
for _, r in miss.iterrows():
    lab = EXTRA.get((r.species, r.tissue, str(r.celltype)))
    if not lab: continue
    sub = ref[(ref.species==r.species)&(ref.cell_type_screen==lab)]
    rs = ratios_for(sub) if len(sub) else []
    steady=[x for x in rs if not any(k in str(x["state"]).lower() for k in ("challenge","disease","tumour","tumor","psoria"))] or rs
    good=[x for x in steady if x["casp1_absent"] or x["ratio"]>=RATIO_STRONG]
    if not good: continue
    best=max(steady,key=lambda x:x["ratio"])
    rows.append({"dataset_id":r.dataset_id,"species":r.species,"tissue":r.tissue,
                 "suspension_type":r.suspension_type,"celltype":r.celltype,"n_donors":int(r.n_donors),
                 "pct_IL18":r.pct_IL18,"pct_CASP1":r.pct_CASP1,
                 "rank_IL18":r.rank_IL18,"rank_CASP1":r.rank_CASP1,"n_ct_stratum":r.n_ct_stratum,
                 "excluded_because":"CASP1 未进分层内最低三分位（相对规则）",
                 "ref_label":lab,"n_ref_sources_supporting":len({g['source'] for g in good}),
                 "best_ref":f"{best['source']}[{best['unit'][:18]}] IL18={best['il18']:.3g}/CASP1={best['casp1']:.3g} ratio={best['ratio']:.3g}"})
out=pd.DataFrame(rows).sort_values(["species","tissue"])
write_tsv(out, ROOT/"results"/"screen_missed_but_reference_supported.tsv")
print(out[["species","tissue","celltype","n_donors","pct_IL18","pct_CASP1","rank_IL18","rank_CASP1","n_ct_stratum","n_ref_sources_supporting"]].to_string(index=False) if len(out) else "none")
