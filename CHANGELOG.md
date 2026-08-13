# CHANGELOG

## v1.0 — 2026-08-13 — pre-registration frozen

- `docs/PREREGISTRATION.md` committed **before** any expression value was
  read from any count matrix. Only cell-annotation metadata had been read
  at that point (see its §8).
- `panel/chemokine_panel_v1.tsv` frozen at 68 genes. Ensembl IDs resolved
  against Ensembl REST (GRCh38), independent of any dataset.
  - `CCL3L1` resolves to three Ensembl gene IDs
    (ENSG00000276085 / ENSG00000277336 / ENSG00000277768). The primary is
    recorded; the ambiguity is documented in the panel `notes` because it
    *is* the hazard the specification names. A miss on this gene must not
    be repaired by folding it into `CCL3`.
- Panel changes require a version bump (`chemokine_panel_v2.tsv`) and an
  entry here.

### Pre-analysis correction to the specification's recalled values

Measured from `meta_10x.txt` of E-MTAB-6701, before analysis:

| quantity | spec recalled | measured |
|---|---|---|
| total cells | ~70,000 | 64,734 |
| dNK1–3 cells | 11,927 | 11,202 (11,932 with dNKp) |
| usable donors | 11 | 7 total, 6 with decidua, 4 with paired blood |

The donor-count difference is the consequential one: it sets the ceiling
on every test in the study. Recorded in `docs/PREREGISTRATION.md` §2 and
`docs/DECISIONS.md`.
