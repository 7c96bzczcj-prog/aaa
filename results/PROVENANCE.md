# Provenance

Generated 2026-08-25 by `preflight/scripts/make_provenance.py`.

Frozen constants for this run live in `preflight/lib/pf_core.py` (`SPEC_VERSION`); the verdict lookup tables are dumped next to the results as `task_[abc]_verdict_table.json` so they can be diffed against the spec.

| deliverable | produced by | inputs |
|---|---|---|
| `manifest.yaml` | `preflight/scripts/geo_scan.py + gsm_meta.py + manual verification (barcode line counts, GEO SOFT, HCA Azul API)` | GEO esearch/esummary/filelist.txt; HCA Azul index/projects and index/files |
| `results/soup_by_compartment.tsv` | `preflight/scripts/task_a.py -> task_a_stats.py` | HCA ICA *_raw_feature_bc_matrix.h5 (dcp60); GSE233304 *_matrix.mtx.gz; GSE181543 *_raw_feature_bc_matrix.h5; cell labels from preflight/scripts/annotate.py |
| `results/soup_by_donor.tsv` | `preflight/scripts/task_a_stats.py` | as above |
| `results/target_gene_detection.tsv` | `preflight/scripts/task_a_stats.py` | admission criterion C2, per dataset x compartment |
| `results/early_nk_stages.tsv` | `preflight/scripts/task_c.py -> task_c_stats.py -> task_c_supp.py -> assemble_task_c.py` | HCA ICA MantonBM1-8 (71 libraries) and GSE233304 marrow samples, unfiltered |
| `results/ltnk_age.tsv` | `preflight/scripts/task_b.py -> task_b_stats.py` | GSE120221 (20 donors, ages 24-84) and HCA ICA MantonBM1-8 (ages 26-52) |
| `results/rho_by_lineage.tsv` | `preflight/scripts/rho_by_lineage.py` | cross-lineage rho, to check whether the NK-specific soup excess reproduces here |
| `VERDICT.md` | `hand-written from the tsv files; preflight/scripts/verify_verdict.py re-checks every number in it` | - |

## Data versions

* HCA Census of Immune Cells, project `cc95ff89-2e68-4a08-a234-480eca21ce79`, catalog `dcp60`; 127 raw h5 files admitted, 3.22 GB. Per-file sha256 in `preflight/scan/hca_ica_files.json`.
* GEO series downloaded from `ftp.ncbi.nlm.nih.gov/geo`; the `filelist.txt` of every scanned series is summarised in `preflight/scan/bm_scan.json` and `preflight/scan/pb_scan.json`.

## Environment

```
python 3.11.15
numpy 2.4.6
scipy 1.17.1
pandas 2.3.3
scanpy 1.11.5
anndata 0.12.19
h5py 3.16.0
statsmodels 0.14.6
scikit-image 0.26.0
```

## Output checksums

| file | sha256 (first 16) | bytes |
|---|---|---|
| `manifest.yaml` | `69ef3acae7c60d0a` | 14471 |
| `results/soup_by_compartment.tsv` | `e15a1f1897375330` | 8402 |
| `results/soup_by_donor.tsv` | `47ce3c402b11eab3` | 52099 |
| `results/target_gene_detection.tsv` | `8311bae42cf83148` | 2381 |
| `results/early_nk_stages.tsv` | `e80f1ec997723a96` | 5205 |
| `results/ltnk_age.tsv` | `d718b8620ca6f820` | 15969 |
| `results/rho_by_lineage.tsv` | `fbfc91fb5d9cc350` | 62387 |
| `VERDICT.md` | `61fb91d12062d7c2` | 12594 |
