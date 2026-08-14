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

## v1.1 — 2026-08-13 — three rule corrections, applied after the v1.0 run

`docs/PREREGISTRATION.md` (v1.0) is untouched and still has exactly one commit.
Amendments are in `docs/PREREGISTRATION_v1.1.md`, each with **its direction of
effect on the results disclosed**, because two of the three increase the
number of positive calls.

- **A1 — R5 becomes per-arm.** `uninterpretable_ceiling` now requires BOTH
  arms above 0.85. A one-sided ceiling compresses the difference toward zero,
  so the effect is a lower bound, not an unreadable number. New columns
  `ceiling_arm_a`, `ceiling_arm_b`, `magnitude_is_lower_bound`.
  Effect: excluded contrasts 86 → 5; 81 retained as direction-only.
- **A2 — inference moves to the empirical null.** BH over the floor-limited
  signed-rank p can never reject at n ≤ 6, as v1.0 §2.2 established before the
  data. BH is now applied to the empirical p from the 405-gene matched null
  (`q_bh_empirical` in T4); `q_bh` is retained and reported but adjudicates
  nothing. α and family size unchanged.
  Effect: 0 → 8 rows at q ≤ 0.05, all at |z| ≥ 5.6.
- **A3 — ruler B is per-gene.** The denominator is the lineage with the
  highest measured CPM for that gene, not a fixed myeloid reference. A fixed
  denominator only tests pickup when that lineage is the gene's source; for
  tissue-sourced genes both gates pick up equally and the ratio is ~1
  regardless of truth. Decidual band 1.253 → **0.0441**; every ambient control
  now lands low against a sensible source. v1.0 verdicts retained as
  `scale_b_verdict_v10_myeloid_only`.

T2's v1.0 columns are frozen in place and order; v1.1 only appends. The test
suite asserts this, plus the per-arm ceiling semantics and the per-gene ruler-B
source assignment.

### Corrections to the v1.0 report

- **D16 retracts the doublet explanation in D10.** It was an inference from a
  pattern that depth-driven pickup produces equally well. Measured: UMI ratio
  flagged/unflagged is **1.31–1.67**, not the ~2.0 doublets require, and the
  NK:T ratio is **unimodal** (bimodality coefficient 0.24–0.47, threshold
  0.555). The evidence supports depth-driven ambient pickup. The
  dNK1-specific residual after depth matching remains unexplained.
  New: `src/doublet_evidence.py`, `out/VT2018/doublet_evidence.{tsv,png}`.
- **D17 verifies PAEP's absence against the 4.1 GB source file**, not just the
  ingested h5ad. Absent at source; CSH1/DCN/IGFBP1/GNLY present; no alias, no
  `ENSG00000122133`; chromosome-9 neighbours all present.
- **D15 audits the n = 11 conflict.** E-MTAB-6701 has one donor column with 7
  values; zero libraries and zero matrix runs map to more than one donor.
  Merging the Smart-seq2 companion (E-MTAB-6678: D3, D5, D6–D9) reaches 9
  donors, not 11. No donor-level statistic from this accession can have
  n > 7, or n > 6 for decidua.
- `docs/RESULTS.md` rewritten under v1.1, and now reports dNK3 in full.

## v1.2 — 2026-08-13 — circularity test, and four corrections that tighten the reading

v1.0's preregistration is still untouched (one commit). v1.1 stands as written.
Amendments in `docs/PREREGISTRATION_v1.2.md`. Unlike v1.1, **three of these
four make the results weaker or narrower**; none was chosen to promote a gene.

- **A4 — authorised exception to the no-re-clustering bar**, for a circularity
  test only. `src/circularity_test.py` re-groups decidual NK on CD39/CD103
  equivalents (ENTPD1/ITGAE); CD160 and KLRB1 are barred as gating genes,
  being dNK3 characterisation genes. No output feeds the main analysis.
- **A5 — one BH family of 81** (27 genes × 3 contrasts), not three of 27. All
  8 rows still clear, two at **q = 0.0499**.
- **A6 — the empirical p is floored at 1/(N+1)**. `q < 0.0001` was
  unsupportable over 405 nulls; six rows now read "≤ 0.0025, at the resolution
  limit". z is kept as an effect size and never converted to a p.
- **A7 — ruler B ranks sources by total pool contribution**, not per-cell CPM,
  and NK being the largest contributor forces `unmeasurable`. XCL1, XCL2,
  CCL5, CXCR4 and CCL4 all leave ruler-B positive; **no headline gene meets
  the "both rulers positive" bar any more.** Band 1.253 → 0.0441 → **0.0538**.

### The circularity result (D18)

Vento-Tormo's subsets are unsupervised whole-transcriptome clusters
characterised by differential expression, and CCL5/CXCR4/XCL1 are among the
characterising genes — so the contamination is partial. Marker-based
regrouping retains 0.62× (XCL1), 0.53× (XCL2), 0.46× (CCL5), 0.30× (CXCR4).
Dilution cannot explain a 2-fold spread ordered that way, and the within-label
test — a gate contrast inside one fixed cluster, where circularity is
impossible — separates them:

- XCL1/XCL2 survive inside dNK1 (+6.7/+5.7 pp, 3/3 donors)
- CCL5 survives inside dNK2 along CD103 (+6.0 pp, 3/3)
- **CXCR4 reverses inside dNK1** (−3.5 pp, 3/3)

**The study's largest z (CXCR4, +14.56) is its least trustworthy row; XCL1
dNK2 > dNK1 is its strongest, and is a cross-modality reproduction of Huhn
2020.** All within-label n = 3 (`min_achievable_p` 0.25), so only direction
and magnitude are claimed.

### Also

- **D19** — the dNK1 triple-positive residual is not library T-cell content:
  Pearson r = +0.061 (p = 0.90) over 7 libraries. Candidate excluded; the
  other two remain undecidable here.
- **D15 addendum** — no library-level count equals 11, but 11 falls between
  "runs containing any dNK1-3" (12) and "runs with ≥30 dNK1-3" (8), so a
  different cell-count threshold reaches it. The libraries-as-donors reading
  stays live and remains the most economical explanation of a clean 11/11.
  **The OXPHOS result should not be cited until the donor column is printed.**

## v1.3 — 2026-08-13 — the FDR criterion is withdrawn; z is shown never to be comparable

- **A8 — no FDR, no replacement rule.** With N nulls per target the rank floor
  is 1/(N+1), and BH at m = 81 needs 15 rows tied at it when N = 112, 6 when
  N = 317, 5 when N = 371 — against 8 candidate rows in the whole study. The
  same data gave 6, then 0, then 8 passes as N alone changed. Rows carry
  `fdr_estimable = FALSE`. Descriptive quantities replace it: effect in
  points, z as a standardised effect size never converted to a p, donor sign
  concordance, ruler-A margin, and the observed effect's rank in its null set.
  Keeping BH as a "reported but non-adjudicating" column (v1.1 A2) was a
  half-measure — a number printed beside a result is read as a verdict.
- **A9 — nulls matched on the reference arm's baseline detection rate**, not
  CPM, computed over the whole transcriptome at matched depth per contrast.
  **It did not fix the S1PR5 problem; it sharpened it** (z +5.01 → +10.34 on
  the same +1.5 pp). A z of 5 buys 32.1 points at >30% baseline detection and
  1.63 points at <0.5% — 20-fold. Recorded as a limit of z on a bounded scale,
  not patched, and no effect-size filter was added afterwards.
- **D23 — ruler A splits by lineage class**, tested on data-derived markers:
  haematopoietic median 0.250 / 1.8% saturated, non-haematopoietic 1.000 / 90%
  (p = 3.1e-19). So PAEP, being glandular epithelial, would be expected to
  return 1.000 and leave the ceiling untouched — **the CXCR4 risk raised in
  D21 is downgraded**. A first pass got this wrong by ranking markers on
  lineage CPM alone, which selects housekeeping genes owned by the largest
  lineage; 5× specificity fixed it.
- Attribution in D21/RESULTS made factual: the retracted PAEP argument is
  attributed to the v1.2 CHANGELOG entry, not to a person.

**The pattern across v1.1–v1.3.** Four instruments withdrawn rather than
tuned — BH over the signed-rank p, ruler B's fixed denominator, the FDR
criterion, z as a standalone ranking — each because its output was set by a
nuisance parameter instead of by the effect. **What survived is exactly what
never needed a null distribution.**
