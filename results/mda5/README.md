# MDA5 / IFIH1 result tables

Computed tables backing `docs/MDA5_IFIH1_CNS_AND_MS.md`.

## Baseline (normal brain)
| file | contents |
|---|---|
| `out_census_baseline_nucleus.csv` | CELLxGENE Census, normal human brain, single-nucleus, per cell type |
| `out_census_baseline_all.csv` | the same including dissociated-cell suspensions |
| `out_census_per_dataset.csv` | per-dataset values, used to show endothelium ranks top in 14/15 datasets |
| `out_hpa_brain_ifih1.csv` | Human Protein Atlas single-nuclei brain panel |
| `out_gtex_ifih1_tissues.csv` | GTEx v8 median TPM across 54 tissues |
| `out_baseline_consolidated.csv` | Census and HPA merged into coarse cell classes |
| `out_180759_baseline.csv`, `out_schirmer_baseline.csv` | per-cohort baselines |

## Disease effect
| file | contents |
|---|---|
| `out_meta_donor_result.csv` | donor-level meta-analysis across 21 brains |
| `out_meta_donor_cp10k.csv` | the per-donor values behind it |
| `out_meta_nb_glm.csv` | negative-binomial pseudobulk model, the conservative estimate |
| `out_lerma_*.csv` | GSE279180 cohort: cell-level, donor-level, depth-matched, microglial subtypes |
| `out_jakel_ms_vs_ctrl.csv` | GSE118257 cohort |
| `out_schirmer_ms_vs_ctrl.csv` | Schirmer 2019 cortex |
| `out_lerma_isg_panel_donor.csv` | the RLR/ISG panel, depth-normalised, donor-level |

## Grey versus white matter
| file | contents |
|---|---|
| `out_gm_wm_depthmatched.csv` | the same depth-matched statistic in both compartments |
| `out_gm_vs_wm_baseline.csv` | baseline hierarchy, each cell type relative to astrocytes in its own dataset |

## Literature
| file | contents |
|---|---|
| `literature_claims_verified.csv` | 145 literature claims across 8 research angles, each adversarially re-verified: 84 confirmed, 57 overstated, 3 refuted, 1 unverifiable |
| `literature_workflow_raw.json` | the full raw output including evidence and verifier reasoning |

## Figure
`mda5_figure.html` (interactive, light and dark), with `mda5_figure_light.png` and
`mda5_figure_dark.png` as static renders.
