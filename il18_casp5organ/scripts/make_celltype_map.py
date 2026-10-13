"""Generate celltype_map.tsv — explicit per-label lineage assignment (§6.2, §7.1). No regex.
Every distinct author label is hand-assigned below; the script errors if any label in the data is unmapped.
"""
import glob, os, sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]

EPI = "parenchymal_epithelial"; NEU="neuron"; GLIA="glia"; RMAC="resident_macrophage"
MMY="monocyte_derived_myeloid"; DC="dendritic_cell"; LYM="lymphoid"; GRAN="granulocyte"
ENDO="endothelial"; MES="mesenchymal"; OTH="other"

M = {
 # epithelial — lung airway/alveolar
 "pulmonary alveolar type 2 cell":EPI,"pulmonary alveolar type 1 cell":EPI,"pulmonary alveolar epithelial cell":EPI,
 "multiciliated columnar cell of tracheobronchial tree":EPI,"multiciliated epithelial cell":EPI,"ciliated cell":EPI,
 "respiratory basal cell":EPI,"basal cell":EPI,"club cell":EPI,"goblet cell":EPI,"bronchial goblet cell":EPI,
 "tracheobronchial goblet cell":EPI,"nasal mucosa goblet cell":EPI,"mucus secreting cell":EPI,"tracheobronchial serous cell":EPI,
 "serous secreting cell":EPI,"respiratory tract hillock cell":EPI,"brush cell of tracheobronchial tree":EPI,"ionocyte":EPI,
 "pulmonary neuroendocrine cell":EPI,"epithelial cell of lower respiratory tract":EPI,"epithelial cell":EPI,"acinar cell":EPI,
 # epithelial — liver/pancreas parenchyma
 "hepatocyte":EPI,"periportal region hepatocyte":EPI,"centrilobular region hepatocyte":EPI,"midzonal region hepatocyte":EPI,
 "cholangiocyte":EPI,"duct epithelial cell":EPI,"pancreatic acinar cell":EPI,"pancreatic ductal cell":EPI,
 "type B pancreatic cell":EPI,"pancreatic A cell":EPI,"pancreatic D cell":EPI,"pancreatic PP cell":EPI,"pancreatic endocrine cell":EPI,
 # neuron
 "neuron":NEU,"neurons":NEU,"pyramidal neuron":NEU,"hippocampal neuron":NEU,"glutamatergic neuron":NEU,"cortical interneuron":NEU,
 "L4/5 intratelencephalic projecting glutamatergic neuron of the primary motor cortex":NEU,
 "L2/3-6 intratelencephalic projecting glutamatergic neuron":NEU,"L6 corticothalamic-projecting glutamatergic cortical neuron":NEU,
 "L6 intratelencephalic projecting glutamatergic neuron of the primary motor cortex":NEU,"sst GABAergic cortical interneuron":NEU,
 "VIP GABAergic cortical interneuron":NEU,"lamp5 GABAergic cortical interneuron":NEU,
 "L5/6 near-projecting glutamatergic neuron of the primary motor cortex":NEU,"L6b glutamatergic cortical neuron":NEU,
 "pvalb GABAergic cortical interneuron":NEU,"L5 extratelencephalic projecting glutamatergic cortical neuron":NEU,
 "corticothalamic-projecting glutamatergic cortical neuron":NEU,"sncg GABAergic cortical interneuron":NEU,
 "L2/3 intratelencephalic projecting glutamatergic neuron of the primary motor cortex":NEU,"sst chodl GABAergic cortical interneuron":NEU,
 "Cajal-Retzius cell":NEU,"meis2 expressing cortical GABAergic cell":NEU,
 "L2/3 intratelencephalic projecting glutamatergic neuron":NEU,"L4/5 intratelencephalic projecting glutamatergic neuron":NEU,
 "L5 intratelencephalic projecting glutamatergic neuron":NEU,"L5/6 near-projecting glutamatergic neuron":NEU,
 "L6 intratelencephalic projecting glutamatergic neuron":NEU,
 # glia
 "astrocyte":GLIA,"astrocytes":GLIA,"oligodendrocyte":GLIA,"Oligos":GLIA,"oligodendrocyte precursor cell":GLIA,"OPC":GLIA,
 "ependymal cell":GLIA,"Bergmann glial cell":GLIA,"Schwann cell":GLIA,
 # resident macrophage
 "alveolar macrophage":RMAC,"Kupffer cell":RMAC,"central nervous system macrophage":RMAC,"microglial cell":RMAC,"microglia":RMAC,
 "perivascular macrophage":RMAC,
 # monocyte-derived / monocyte
 "classical monocyte":MMY,"intermediate monocyte":MMY,"non-classical monocyte":MMY,"monocyte":MMY,"elicited macrophage":MMY,
 "inflammatory macrophage":MMY,"lung macrophage":MMY,"macrophage":MMY,"myeloid leukocyte":MMY,
 # dendritic
 "dendritic cell":DC,"conventional dendritic cell":DC,"plasmacytoid dendritic cell":DC,"CD1c-positive myeloid dendritic cell":DC,
 "myeloid dendritic cell":DC,
 # lymphoid
 "T cell":LYM,"CD4-positive, alpha-beta T cell":LYM,"CD8-positive, alpha-beta T cell":LYM,"regulatory T cell":LYM,
 "mature NK T cell":LYM,"natural killer cell":LYM,"B cell":LYM,"plasma cell":LYM,
 # granulocyte
 "neutrophil":GRAN,"mast cell":GRAN,"basophil":GRAN,
 # endothelial
 "endothelial cell":ENDO,"capillary endothelial cell":ENDO,"vein endothelial cell":ENDO,"pulmonary artery endothelial cell":ENDO,
 "endothelial cell of lymphatic vessel":ENDO,"blood vessel endothelial cell":ENDO,
 "endothelial cell of pericentral hepatic sinusoid":ENDO,"endothelial cell of hepatic sinusoid":ENDO,"endo/vascular":ENDO,"endothelial cell of artery":ENDO,
 # mesenchymal
 "fibroblast":MES,"fibroblast of lung":MES,"bronchus fibroblast of lung":MES,"pulmonary interstitial fibroblast":MES,
 "alveolar adventitial fibroblast":MES,"alveolar type 1 fibroblast cell":MES,"adventitial cell":MES,"myofibroblast cell":MES,
 "hepatic stellate cell":MES,"pancreatic stellate cell":MES,"pericyte":MES,"lung pericyte":MES,
 "smooth muscle cell":MES,"tracheobronchial smooth muscle cell":MES,"bronchial smooth muscle cell":MES,
 "smooth muscle cell of the pulmonary artery":MES,"vascular associated smooth muscle cell":MES,"stromal cell":MES,
 "chondrocyte":MES,"vascular leptomeningeal cell":MES,"vascular leptomeningeal cell (Mmus)":MES,
 # other / unresolved
 "unknown":OTH,"leukocyte":OTH,"hematopoietic stem cell":OTH,"hematopoietic cell":OTH,"erythroid lineage cell":OTH,
 "enucleate erythrocyte":OTH,"mesothelial cell":OTH,"choroid plexus epithelial cell":OTH,"meninges/schwann":OTH,"cell":OTH,
}
# GSE103892 fine neuron subtypes DE-/DI-/VE-/VI-/VC-/VM-/ME-/MI- (+ 'DE12') -> neuron (each listed explicitly below)
for pre,rng in [("DE",range(1,17)),("DI",range(1,10)),("VE",range(1,5)),("VI",range(1,6)),("VC",range(1,3)),("VM",range(1,3)),("ME",range(1,2)),("MI",range(1,5))]:
    for i in rng: M[f"{pre}-{i}"]=NEU
M["DE12"]=NEU

rows=[]; unmapped=set()
for f in sorted(glob.glob(str(ROOT/"data/cells/*.parquet"))):
    eid=os.path.basename(f)[:-8]
    for lab in pd.read_parquet(f).celltype_original.unique():
        if lab not in M: unmapped.add(lab)
        else: rows.append((eid,lab,M[lab]))
if unmapped:
    print("UNMAPPED:",sorted(unmapped)); sys.exit(1)
df=pd.DataFrame(rows,columns=["entry_id","celltype_original","lineage"]).sort_values(["entry_id","lineage","celltype_original"])
with open(ROOT/"celltype_map.tsv","w") as fh:
    fh.write("# celltype_map.tsv — explicit per-label lineage mapping (§6.2, §7.1); no regex.\n")
    fh.write("# resident_macrophage lumps alveolar mac / Kupffer / microglia; cite level-2 (celltype_original) for any claim about these (§7.1).\n")
    df.to_csv(fh,sep="\t",index=False)
print("wrote",len(df),"rows,",df.entry_id.nunique(),"entries")
