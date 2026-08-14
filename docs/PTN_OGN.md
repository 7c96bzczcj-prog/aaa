# Do decidual NK cells transcribe PTN and OGN?

Three checks, ordered by cost. Checks 2 and 3 are complete and are reported
here; check 1 is running and its cost premise needed correcting first.

The short version: **on every readout reachable from here, PTN and OGN in dNK
sit inside the carryover envelope, while SPP1/OPN — the third factor from the
same paper — sits three orders of magnitude above it.** The claim is not
uniformly wrong; it separates by gene.

One paper that could overturn this could not be read. It is named in §5.

---

## 0. What each check was supposed to decide, and what it decided

| check | question | cost as expected | cost as measured | verdict |
|---|---|---|---|---|
| 1 | does NK PTN/OGN detection track tissue stromal content? | free, data on disk | **not on disk**; ~560 MB stream, one ruler only | running |
| 2 | do stromal genes appear in their sorted dNK RNA? | one paper read | one paper read **+ their own deposited RNA-seq** | **yes — contamination demonstrated** |
| 3 | four methodological gates on the induced-NK paper | one paper read | one paper read + their own deposited RNA-seq | **two gates fail, one is weak, one passes** |

---

## 1. Check 1 — the cost premise was wrong, twice

**"Already in your archive" is not true.** `data/` holds VT2018 and nothing
else (4.2 GB). Only the atlas *`obs` table* was ever read, remotely, as
`D1` records — 8.8 MB of a 1,765 MB file, enough to establish that the atlas
contains no decidual sample and to close it as a dataset. No expression value
from it has ever been on this disk.

**And the two rulers cannot both be carried.** `all_nk_cells.h5ad` is
NK-only: 89,216 cells, `obs` columns `sample`/`dataset`/`source`, no stromal
compartment anywhere in the file. So:

- **Ruler B is unavailable.** Its denominator is the source lineage's own
  expression. There are no fibroblasts in the file to measure.
- **Tissue stromal content is unavailable as an in-dataset covariate**, for
  the same reason — there are no stromal cells to count. The per-tissue
  all-cells h5ads on the same Zenodo record do carry them, but total 41 GB
  against 25 GB of free disk.
- **Ruler A survives, in its prior-selected form only.** COL1A1/COL1A2/
  COL3A1/DCN/LUM cannot be transcribed by an NK cell, so their entire signal
  inside an NK barcode is pickup. That is an ambient index measured in the NK
  cells themselves, and it is admissible for the reason set out in
  `TRANSFERABLE.md` §5: the controls come from prior exclusivity, not from a
  statistic, and they cross a boundary the target cannot cross.

So check 1 runs with **one ruler, not two**, and that limitation is a property
of the file rather than a choice. What it can still deliver is the designed
test in its usable form: across 61 samples in 14 tissues, does PTN/OGN
detection in NK track the collagen-pickup index in the same cells?

**Cost as measured:** the counts layer is CSR over 89,216 × 11,866 with
70,466,919 non-zeros, so two gene columns cannot be sliced out — the row index
must be walked. `src/netskar_probe_stream.py` walks it once over HTTP range
requests, keeping per-cell totals and 16 probe genes, ~560 MB of transfer and
no disk. All 16 panel genes are present (R11 gate passed at 100%).

---

## 2. Check 2 — their own sorted dNK libraries carry stroma

GSE184719 is the bulk RNA-seq behind the 2022 induced-NK paper, from the same
group. Twenty libraries: `pNK1-4` (blood NK), `dNK1-4` (sorted decidual NK),
and twelve induced-NK cultures. The sorted dNK libraries answer check 2
directly.

### Stromal transcripts in sorted dNK, CPM

| gene | dNK1 | dNK2 | dNK3 | dNK4 | pNK mean |
|---|---|---|---|---|---|
| COL1A1 | 149.6 | 388.4 | 43.2 | 28.7 | 0.73 |
| COL1A2 | 46.3 | 223.0 | 12.6 | 13.6 | 0.07 |
| COL3A1 | 33.7 | 493.9 | 32.2 | 45.5 | 0.21 |
| DCN | 108.5 | 1496.3 | 95.0 | 44.4 | 0.90 |
| LUM | 66.6 | 1807.7 | 140.8 | 57.7 | 1.12 |
| **sum** | **404.7** | **4409.3** | **323.7** | **189.9** | **3.04** |
| PAEP (glandular epithelium) | 1866.2 | 59.9 | 366.3 | 138.0 | 0.76 |

Every one of the four sorted dNK libraries carries stromal transcript at
**60–1450× the blood-NK baseline**, and a second contaminating compartment
(PAEP, glandular epithelium) besides. Surface-marker purity did not exclude
either, which is the point: protein sticks to cells and free RNA rides along,
and neither is visible to the sorter.

**Contamination is therefore not a hypothesis about this data. It is
measured, in their own deposit.** What remains is whether it is *enough* to
account for the PTN/OGN they report.

### Is it enough? — ruler B, platform-matched, inside VT2018 alone

The GSE184719 libraries have no stromal arm of their own, so the source
lineage has to come from somewhere. VT2018 has it: 14,872 decidual stromal
cells alongside 11,198 dNK, same tissue, same platform, one matrix.

Ruler B: a gene's level in the NK gate as a percentage of its level in
decidual stroma. Genes NK cannot transcribe fix the envelope; a gene NK
really transcribes must sit above it.

| gene | dNK1 | dNK2 | dNK3 | reading |
|---|---|---|---|---|
| COL1A1 | 2.00 | 1.06 | 2.67 | envelope |
| COL1A2 | 1.92 | 1.25 | 3.64 | envelope |
| COL3A1 | 2.20 | 1.11 | 2.85 | envelope |
| DCN | 0.61 | 1.54 | 0.71 | envelope |
| LUM | 0.54 | 1.42 | 0.54 | envelope |
| **PTN** | **0.67** | **0.97** | **0.21** | **inside — no excess to explain** |
| **OGN** | **1.52** | **1.61** | **1.83** | **inside — no excess to explain** |
| **SPP1** (OPN) | **8593** | **5437** | **11350** | **far above — genuinely transcribed** |
| VIM | 22.9 | 35.9 | 35.7 | above — which is why it is not a control |
| PTPRC | 65324 | 77657 | 90626 | reciprocal control, as expected |

Carryover envelope: **0.54 % – 3.64 %** of stromal level. PTN and OGN fall
inside it in every arm; PTN/dNK3 falls *below* it.

**SPP1 is the positive control and it works.** The same instrument, on the
same cells, puts the third factor of Fu 2017 three orders of magnitude above
the envelope. This is not a method that returns "contamination" for
everything asked of it — it separates the three genes, and it separates them
the same way twice.

**VIM is the negative control for control choice, and it also works.** At
23–36 % it is far outside the envelope, confirming that leukocytes transcribe
vimentin and that it must never be used as a contamination marker — which is
exactly how it is treated in the source scripts.

### The primary readout (R3), for completeness

Detection rate, raw counts — an upper bound on the depth-matched rate:

| gene | dNK1 | dNK2 | dNK3 | Stromal |
|---|---|---|---|---|
| PTN | 0.34 % | 0.39 % | 0.10 % | 58.0 % |
| OGN | 0.68 % | 0.50 % | 0.65 % | 52.3 % |
| DCN | 13.5 % | 18.3 % | 12.7 % | 99.8 % |
| LUM | 16.3 % | 20.9 % | 13.5 % | 99.6 % |
| SPP1 | 60.5 % | 29.1 % | 62.6 % | 3.5 % |

PTN and OGN sit in the **floor band** (< 5 %), the same band as CCR5 at
0.67 % — so on the primary scale they have no power at all in this design,
and that is stated rather than discovered later (v1.4 A11/A12).

The internally consistent part: DCN and LUM are detected 20–50× more often
than PTN and OGN in the same cells, and in stroma they are 50–60× more
abundant (5535 / 7460 CPM vs 118 / 90). The ratio matches. Pure carryover
predicts exactly this spread, and it is what is there.

### What does not carry a verdict

`out/GSE184719/ptn_ogn_carryover_prediction.csv` extends the VT2018 stromal
ratios to predict PTN/OGN in the GSE184719 bulk libraries. It is reported and
**adjudicates nothing**, because the implied stromal fraction moves from
0.77 % to 154 % depending only on which anchor gene is used — and a fraction
above 100 % is not a quantity. 10x 3′ and bulk poly-A do not measure long
collagens on the same scale. This is the `TRANSFERABLE.md` §1 failure form
caught before it was used, not after: the output was a function of the anchor
choice rather than of the contamination.

The verdict above rests on the within-VT2018 envelope, where no cross-platform
step exists. No count matrix was merged (R9).

### The mouse arm

From the abstracts of Fu 2017 and of the companion Sci Transl Med 2020
(Zhou, Fu, …, Wei; PMID 32238574), the mouse genetics are **`Nfil3`^−/−,
`Tbx21`^−/−, adoptive transfer of induced CD49a⁺Eomes⁺ NK, and "inactivation
of `Pbx1` in mouse dNK cells"** — every one of them at the level of the
transcription factor or the subset. A literature search returns no `Ptn`^fl/fl
or `Ogn`^fl/fl NK-conditional allele anywhere.

So the answer to the question as posed is: **subset- and TF-level
manipulation, not conditional deletion of the two factors in NK.** That
design can establish that the subset matters for fetal growth; it cannot
establish that these two factors come from NK cells.

**Marked uncertainty.** Both papers are paywalled with no PMC deposit and
cell.com returns 403, so the Methods sections could not be read. The statement
above is inference from abstracts plus a negative literature search. It should
be checked against the STAR Methods by anyone with access before it is
relied on.

---

## 3. Check 3 — the induced-NK paper, four gates

Du et al., *Front Immunol* 2022;13:823227 (PMC8854499), "Human-Induced CD49a⁺
NK Cells Promote Fetal Growth".

| gate | finding | verdict |
|---|---|---|
| transport inhibitor | *"the cells were cultured for 4 h in the presence of monensin (2.5 μg/mL; Sigma)"* | **passes** |
| antibody validation | anti-PTN LS-C162291, anti-OGN LS-B10948 (LifeSpan), 1:50; isotype control rabbit IgG (CST 3900S). No knockout, no blocking peptide | **weak** |
| supernatant ELISA | none. The claim is *secretion*; every measurement is intracellular | **fails** |
| mRNA in stroma-free culture | measurable in their own deposit — see below | **fails** |

### The fourth gate, measured rather than looked up

The twelve induced-NK libraries in GSE184719 are feeder-free cultures from
cord-blood HSC, bone-marrow HSC, and blood NK. **There is no stroma in them.**
PTN or OGN mRNA there could not be contamination, so this is the decisive
experiment, and it was already deposited.

**PTN, raw reads across all twelve stroma-free libraries: 0, 0, 0, 0, 1, 0, 0,
0, 0, 0, 2, 0 — three reads in total.** In CPM: 0.00–0.11, against 1.72–15.08
in the four sorted dNK libraries.

**OGN: 0.00–0.47 CPM, indistinguishable from peripheral-blood NK (0.16–0.32
CPM)** — which is the group the paper contrasts the induced cells *against*.

This is a **powered** negative, not a floor artefact. Libraries are 12.6–24.7 M
reads, so one read is ≈ 0.04–0.08 CPM; the assay resolves two orders of
magnitude below the dNK level, and it sees PTN perfectly well in the dNK
libraries from the same submission. Absence here is a measurement.

### The contradiction this creates

The paper reports PTN and OGN **protein** in these cultures by intracellular
staining (Fig. 4A–D). Its own RNA-seq of the same cultures shows **no PTN
message and OGN at blood-NK background**. Protein without message, with an
antibody validated by isotype control alone and no knockout or blocking
peptide — the antibody is the unsupported step, and gate 2 is where the claim
breaks rather than gate 1.

Consistent with this, the paper's own growth-factor heatmap (Fig. 4E) shows
VEGFA, LIF, IL-32 and CSF2 — **not PTN and not OGN.** The transcriptional
evidence for the two headline factors is absent from the figure that would
have carried it.

---

## 4. What the three checks together support

**Separating by gene, which the original claim does not:**

- **SPP1 / osteopontin** — genuinely transcribed by dNK. 5,000–11,000 % of
  stromal level, detected in 29–63 % of dNK cells against 3.5 % of stromal
  cells. Nothing here touches it.
- **PTN and OGN** — every reachable measurement puts them inside the carryover
  envelope, in the floor band on the primary readout, and absent from
  stroma-free cultures of the cells claimed to make them.

**What this does not establish.** That dNK never transcribe PTN or OGN. A
floor-band readout has no power in points (v1.4 A12), so "inside the envelope"
means *not distinguishable from carryover*, not *demonstrated absent*. The
honest label is the same one CCL3 got: **neither established nor excluded** —
except that here, unlike CCL3, there is a positive control in the same table
showing the instrument can see the thing when it is there.

---

## 5. The open item that could overturn this

Zhou Y, Fu B, et al. **PBX1 expression in uterine natural killer cells drives
fetal growth.** *Sci Transl Med* 2020;12(537):eaax1798 (PMID 32238574).

Its abstract states: *"PBX1 drives pleiotrophin and osteoglycin transcription
in dNK cells."* If that rests on ChIP, a promoter–reporter assay, or PBX1
overexpression in a stroma-free NK line with an mRNA readout, it is direct
transcription-level evidence and it **outranks everything above** — a
contamination argument cannot survive a promoter binding to a gene the cell
supposedly does not transcribe.

**It is paywalled, has no PMC deposit, and science.org returns 403. It could
not be read.** Anyone with access should read that paper's transcription
experiments before acting on this file. It is named here rather than left out
because leaving out the strongest contrary evidence is how a review becomes an
argument.

---

## 6. Provenance

| output | script | input |
|---|---|---|
| `out/GSE184719/panel_cpm_per_library.csv` | `src/gse184719_ptn_ogn.py` | GSE184719 `gene_counts.txt.gz` |
| `out/GSE184719/panel_rawcounts_per_library.csv` | same | same |
| `out/GSE184719/vt2018_lineage_cpm.csv` | `src/stromal_ruler.py` | VT2018 canonical h5ad |
| `out/GSE184719/vt2018_nk_over_stroma_pct.csv` | same | same |
| `out/GSE184719/ptn_ogn_carryover_prediction.csv` | same | **reported, adjudicates nothing** |
| `out/NETSKAR2024/cell_probe_counts.csv.gz` | `src/netskar_probe_stream.py` | remote, streamed |

R10 held throughout: anndata + numpy + scipy + pandas + h5py, no scanpy.
