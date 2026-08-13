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

## v1.0 — 2026-08-13 — VT2018 executed

Pipeline built (`src/`, manifest-driven, no dataset branches) and run on the
anchor dataset. Outputs T1–T6 in `out/VT2018/`, robustness arm in
`out/VT2018/robustness_tpurged/`.

- Ingest: 64,734 cells × 31,764 genes, 694,826,706 counts, verified integer,
  round-tripped through the canonical h5ad.
- Panel match 98.5% (67/68). `PAEP` is absent from the published matrix
  entirely (D11); `CCL3L1` matched via Ensembl only (D9).
- **The purity stop fired** (dNK1 triple-positive 17.35% > 15%). Diagnosed
  before continuing: depth confounding plus NK–T doublets, not mislabelled
  T cells (D10). Both arms run and retained.
- **R6 sensitivity became mandatory** — retention differs by 43.1 pp between
  dNK1 and dNK3 in donor D8 (D13). All comparisons repeated at three depth
  floors.
- Rulers re-calibrated in-dataset: ambient ceiling 0.092 (prior 0.402),
  true-NK floor 0.021 (prior 0.012–0.013), pickup band 1.253 (prior 0.118).
  Ruler B is shown to have little power in decidua; the conservative band was
  kept (D12).
- Empirical null measured, not assumed: +0.47 ± 4.05 pp over 405
  expression-matched genes.
- **No gene reaches q ≤ 0.05**, as the pre-registration established was
  near-impossible at n = 6 (D14). Results are reported as effect sizes,
  sign concordance and empirical-null z.

No change was made to `docs/PREREGISTRATION.md` at any point after freezing.

### Pipeline additions beyond the original six scripts

- `src/ingest.py` — source → canonical h5ad, so analysis code never parses a
  source format and every script sees identical input.
- `src/purity_diagnosis.py` — what the purity stop requires: diagnose before
  proceeding. Reports only; re-labelling stays barred.
- `--drop-triple-positive` / `--outdir` on the analysis scripts, for the
  both-ways robustness run the specification's step 4 requires.
