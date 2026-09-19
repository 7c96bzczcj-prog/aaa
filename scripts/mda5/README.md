# MDA5 / IFIH1 analysis scripts

Run order (all expect the GEO/Census downloads described in
`docs/MDA5_IFIH1_CNS_AND_MS.md` under Methods):

| script | what it does |
|---|---|
| `libsize.py` | one full pass over the GSE180759 dense matrix to get per-nucleus total counts (needed for CP10k; the matrix ships no library sizes) |
| `final180759.py` | GSE180759: depth check, baseline per cell type, donor-level CP10k, depth-matched re-test, lesion-region gradient, RLR/ISG specificity, P2RY12+ composition control, ADAR panel |
| `lerma_analyze.py` | GSE279180: chunked extraction of the gene panel plus per-nucleus totals from the nine per-cell-type h5ad files |
| `census_ifih1.py` | CELLxGENE Census pull of IFIH1 across all human brain cells (normal + MS) |

Gene panel extraction from the GSE180759 matrix is a one-line awk filter on the
first CSV field; see the Methods section. Outputs land in `results/mda5/`.
