# GSE180759 — database retrieval record

Date: 2026-09-19. Retrieved directly from NCBI (GEO SOFT via FTP, E-utilities
for PubMed / PMC / SRA). Raw output saved under `raw_meta/`:

- `raw_meta/GSE180759_self.txt` — series block of the GEO SOFT family file
- `raw_meta/GSE180759_gsm.txt` — all 20 sample blocks
- `raw_meta/GSE180759_annotation.txt.gz` — the authors' per-nucleus annotation
  (462 KB, verbatim copy of the GEO supplementary file)

## 1. The citation and the accession are two different papers

The request paired *Nat Neurosci* 2023;26(10):1701–1712
(doi 10.1038/s41593-023-01435-z) with GSE180759. They do not belong to the
same study:

| | paper | relation to GSE180759 |
|---|---|---|
| doi 10.1038/s41593-023-01435-z | **Andreadou M, … Mundt S, Becher B.** *IL-12 sensing in neurons induces neuroprotective CNS tissue adaptation and attenuates neuroinflammation in mice.* Nat Neurosci 2023;26(10):1701–1712. PMID 37749256, PMC10545539. University of Zurich. | **Re-analysis only.** Methods section "Reanalysis of single-nucleus RNA sequencing of multiple sclerosis tissue": the public expression matrix was downloaded from GSE180759 and re-clustered in Seurat to score `IL12RB1` / `IL12RB2` (their Fig. with n = 66,432 nuclei). The paper's own data are mouse EAE (flow, bulk and snRNA-seq). |
| doi 10.1038/s41586-021-03892-7 | **Absinta M, Maric D, … Calabresi PA, Reich DS.** *A lymphocyte-microglia-astrocyte axis in chronic active multiple sclerosis.* Nature 2021;597:709–714. PMID 34497421. NIH/NINDS. | **The depositing study.** GEO lists PMID 34497421 as the series' publication. |

The 2023 paper's data-availability statement (PMC10545539) confirms the split:
its own sequencing is deposited as **GSE236464** (bulk RNA-seq) and
**GSE236540** (mouse snRNA-seq), with count files on Mendeley
(10.17632/zgr9bj57r4.2); "Data from ref. 13 were accessed at GSE180759" and
ref. 14 (Schirmer 2019) at PRJNA544731. Its human re-analysis used the
**filtered count matrix as deposited**, log-normalised per nucleus, Seurat v4,
integrated by individual, 30 PCs.

So GSE180759 is the Absinta 2021 human MS snRNA-seq dataset; the 2023 paper is
one of its ~36 re-users (section 5).

## 2. What GSE180759 is (from the SOFT record)

| field | value |
|---|---|
| title | A lymphocyte-microglia-astrocyte axis in chronic active multiple sclerosis |
| organism / tissue | *Homo sapiens*, autopsy brain white matter |
| assay | **single-nucleus** RNA-seq, 10x Chromium Single Cell 3′ **v3**, 5,000 nuclei loaded per sample; sucrose-gradient nuclei prep from snap-frozen tissue |
| design | 5 progressive-MS autopsy brains (lesion core, chronic-active edge, chronic-inactive edge, periplaque WM) + 3 non-neurological control WM |
| samples | **20 GSMs**, two submission batches: 12 on 24 Jul 2021 (GSM5470483–94), 8 added 16 Feb 2022 (GSM5903107–14) |
| platforms | GPL16791 HiSeq 2500 (14 libraries), GPL24676 NovaSeq 6000 (6 libraries) |
| processing | Cell Ranger, pre-mRNA GRCh38-3.0 |
| status | public 1 Sep 2021, last updated 17 Feb 2022 |
| raw reads | SRA SRP329679 / BioProject PRJNA749443, 38 runs, paired-end, ≈ 264 GB of `.lite` SRA objects |
| processed | `GSE180759_annotation.txt.gz` (462 KB) and `GSE180759_expression_matrix.csv.gz` (113 MB gz; dense CSV, genes × 66,432 nuclei, ≈ 3.5 GB decompressed) |
| contact | Martina Absinta, NIH |

Sample-level characteristics carry only two fields, `pathological_tissue_location`
and `clinical phenotype`; the donor is **not** in the GEO sample metadata. Donor
identity comes from the annotation file instead.

### Processed-file structure (verified)

- Annotation: 66,432 rows, columns `nucleus_barcode, NBB_case, pathology,
  seurat_cluster, cell_type`. Barcodes are unique; the `-N` suffix (N = 1…20)
  is the library index, and each suffix maps to exactly one donor × region.
- Expression matrix: header = the same 66,432 barcodes (no duplicates), rows
  are gene symbols, integer counts. It is the **filtered, annotated** matrix
  — not raw droplets — so empty-droplet-based ambient estimation would need
  the SRA FASTQs re-run through Cell Ranger.
- The GSM titles (`s1`…`s23`, with gaps) cannot be mapped to the barcode
  suffixes from GEO alone; the library table below is from the annotation
  file.

## 3. Composition (from the annotation file)

Cell types, all nuclei:

| cell_type | nuclei | Seurat clusters |
|---|---|---|
| oligodendrocytes | 45,697 | 0,1,2,3,6 |
| astrocytes | 8,211 | 4,9,12,14 |
| immune (microglia / myeloid) | 5,432 | 5,10,17 |
| neurons | 2,803 | 7,15 |
| OPC | 2,175 | 8 |
| vascular cells | 1,708 | 11,13 |
| **lymphocytes** | **406** | 16 |

Libraries (barcode suffix → donor × region), with the two immune labels:

| lib | nuclei | donor | region | immune | lymphocytes |
|---|---|---|---|---|---|
| 1 | 4,597 | 13_047 | chronic_active edge | 596 | 45 |
| 2 | 2,180 | 09_034 | chronic_active edge | 312 | 51 |
| 3 | 3,980 | 13_047 | chronic_active edge | 738 | 27 |
| 4 | 5,830 | 12_078 | chronic_active edge | 742 | 48 |
| 5 | 439 | 09_067 | chronic_active edge | 12 | 2 |
| 6 | 4,496 | 13_047 | chronic_active edge | 473 | 46 |
| 7 | 1,384 | 09_067 | lesion core | 285 | 12 |
| 8 | 1,640 | 13_015 | chronic_inactive edge | 79 | 1 |
| 9 | 3,717 | 13_047 | chronic_inactive edge | 128 | 18 |
| 10 | 3,987 | 13_015 | chronic_inactive edge | 266 | 31 |
| 11 | 4,369 | 09_034 | chronic_inactive edge | 290 | 10 |
| 12 | 1,797 | 12_078 | chronic_inactive edge | 84 | 6 |
| 13 | 1,795 | 12_002 | control WM | 142 | 8 |
| 14 | 2,269 | 14_043 | control WM | 137 | 5 |
| 15 | 1,227 | 11_69 | control WM | 10 | 3 |
| 16 | 2,955 | 13_047 | lesion core | 366 | 12 |
| 17 | 6,194 | 09_034 | periplaque WM | 79 | 0 |
| 18 | 5,710 | 09_067 | periplaque WM | 175 | 13 |
| 19 | 2,729 | 13_015 | periplaque WM | 150 | 14 |
| 20 | 5,137 | 12_078 | periplaque WM | 368 | 54 |

Donors: 5 MS (09_034, 09_067, 12_078, 13_015, 13_047; 2–3 regions each, so
within-donor region contrasts exist) and 3 controls (11_69, 12_002, 14_043;
one region each). Lymphocytes per donor: 148, 108, 61, 46, 27 (MS) and 8, 5, 3
(controls). Donor 13_047 alone contributes 3 of the 6 chronic-active-edge
libraries and 30% of all nuclei.

## 4. Against this project's admission criteria (A1–A4)

Recorded for completeness; the registry (`registry.csv`) is **not** modified,
because the dataset is not a tumour / adjacent-normal design and its intended
use here has not been stated.

| criterion | reading |
|---|---|
| A1 paired individuals | Not applicable as written (no tumour). Within-donor pairing of lesion regions exists for 5 donors; the 3 controls are unpaired. |
| A2 condition ⟂ library | **Fails.** Every library is one donor × one region; no hashing, no shared emulsion. |
| A3 all lineages in one droplet pool | Unsorted nuclei, so glia and immune share a pool — but the protocol's five *immune* lineages are not resolvable: the authors annotate one `immune` (myeloid/microglia) and one `lymphocytes` cluster, with no NK / CD8 / CD4 / B split. |
| A4 NK ≥ 30 per group | **Fails by a wide margin.** All lymphocytes together are ≤ 54 per library and 406 in total; NK is an unannotated fraction of that. |
| raw droplet matrices | No (filtered matrix only); FASTQs in SRA. |

Relevance to the README's "suggested next step 2" (a dissociation-free readout
that still resolves NK): this is the canonical public snRNA-seq of inflamed
human white matter and it is dissociation-free, but it carries lymphocytes at
~0.6% of nuclei and NK at an unmeasurable fraction of those. It cannot play
that role. It does show what snRNA-seq of tissue does to lymphocyte recovery —
the same point as "enrichment, not cohort size" from the GEO rescan.

## 5. Who has re-used GSE180759

PubMed full-text (PMC) search for the accession returns **36** articles
(2021–2026). Those most relevant to the immune side:

- Absinta 2021 *Nature* (the deposit); Andreadou 2023 *Nat Neurosci*
  (IL-12R expression, the cited paper)
- Lerma-Martin 2024 *Nat Neurosci* — spatial + snRNA-seq of subcortical MS
  lesions (tissue niches)
- IMSGC 2023 *Nature* — MS severity locus, CNS resilience
- Maggi 2023 *EBioMedicine* — B-cell depletion does not resolve chronic
  active lesions
- Fagiani 2025 *Nat Commun*, Hyvärinen 2025 *J Neuroinflammation*, Du 2025
  *Nat Immunol* (CSF1R / microglia), Gonzalez Cruz 2026 *J Neuroinflammation*
  (FKBP5 / myeloid)
- Pernin 2024 *Nat Commun* is the only one PubMed indexes with the accession
  in its abstract/metadata; the rest mention it in full text only.

None of these re-analyses isolates NK cells from this dataset.
