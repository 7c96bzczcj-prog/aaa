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

## v1.4 — 2026-08-13 — the effect scale is fixed before any further result

- **A10 — percentage points are primary; z is auxiliary and never read alone.**
  Settled now, before more results exist, because with two scales and no rule
  every later result could be reported on whichever flatters it. No
  transformation fixes this: a logit would push the reading *further* toward
  low-baseline genes. No table in RESULTS is sorted by z.
- **A11 — the measurable band is declared.** 14 of 27 primary targets sit below
  5% detection (no power in points), 12 in 5–85%, 1 above 85%. **The two
  failure modes close on the same genes:** four of the 14 floor-band genes —
  S1PR5, SELL, CCR9, CXCR6 — are exactly those that produced v1.3's largest and
  most misleading z values. No power in points and false power in z, together.
- **A12 — a flat result now carries its band**: untested (<5%), a real null
  (5–85%), or uninterpretable (>85%).
- **Correction.** The reading that CCL3's flat result "may be no power rather
  than no difference" is not supported and not adopted: CCL3 is at 0.817,
  inside the measurable band. CCL4 (0.944) is ceiling-limited, a different
  problem whose remedy is opposite (more depth worsens it).
- **D26 — CXCR4's ceiling margin quantified.** 3 of 56 haematopoietic controls
  fall below its 0.048; adding 5 more has a ~24% chance of flipping its ruler-A
  verdict. 0 of 56 fall below XCL1's 0.021. CXCR4 is now negative on all three
  measured axes (regrouping 0.30×, within-label reversal, ceiling margin);
  XCL1 passes all three.

## v1.4.1 — 2026-08-13 — transferable findings extracted; two corrections

- **`docs/TRANSFERABLE.md`** added: the failure form that recurred four times
  (output determined by a nuisance parameter, not by the effect), the
  bounded-readout pathology and its three-step pre-emption, the same pathology
  in flow-cytometry %positive, and the one-number band check.
- **D28 — CCR5 is in the floor band** at 0.67% (68/10,130 decidual dNK1–3
  cells). Flagged because any archived analysis resting on a CCR5 detection
  rate from a droplet platform was operating in the band that has no power in
  points and false power in z. **The number is explicitly not transferable
  across species/tissue**; the check is. MC38 data is not in this workspace.
- **Correction to D26.** XCL1's "0.0% chance" was a point estimate presented
  as zero. Rule of three gives ~5.4% per control, ~24% over five — the same as
  CXCR4. XCL1's real defence is the 4.4× gap to the observed minimum (0.043)
  versus CXCR4's 1.1×. Also noted: treating "add a control" as random sampling
  is a communication device, not a fix; the lower-prediction-bound
  construction remains the repair.

## v1.4.2 — 2026-08-13 — the ceiling definition set is declared; D26's probabilities withdrawn

- **D29.** Two ruler-A margins had been quoted against two different minima in
  one sentence (0.092 for XCL1, 0.043 for CXCR4). The **definition set is now
  declared: the 9 pre-registered panel ambient controls, ceiling 0.0920.**
  Both margins recomputed against it: **XCL1 4.4×, CXCR4 1.9×**.
  The 56 data-derived haematopoietic markers are a *diagnostic* set and may not
  define the ceiling — the three of them below CXCR4 include **`CCL3`, a primary
  target of this study**, and a gene under test cannot set the threshold that
  judges genes under test. Note the direction: this is the choice that lets
  CXCR4 pass.
  **D26's probability table is withdrawn** — controls are specific genes, not
  exchangeable draws, and the draws it imagined adding already existed. It is
  replaced by a sharper statement: CXCR4's ruler-A verdict is
  definition-dependent (1.9× under set A, 0.89× under set B) and XCL1's is not
  (4.4× and 2.1×). That is a fourth axis on which CXCR4 is negative.
- **D30.** On the decidual side CCR5 is not merely floor-band: every dNK subset
  (0.15–1.43%) sits below both T (4.00%) and Myeloid (2.51%), its dominant pool
  contributor is Myeloid at 57.5%, and NK is not the largest contributor. The
  supported statement is **CCR5 expression is not detectable in dNK at
  transcript level in this dataset**. Explicitly not transferable to mouse
  tumour NK.

## v1.4.3 — 2026-08-13 — the fourth axis against CXCR4 is retracted

- **D31.** All three set-B markers below CXCR4 — `CCL3` (90.9% detected in NK),
  `IFNG`, `LINC00861` — are genes **NK expresses itself**. Their low soup
  fractions are correct readings carrying no ceiling information, so set B never
  delivered an adverse verdict on CXCR4 at all. Confirmed by recomputation:
  screening on NK/source ≤ 0.05 (the criterion IGKC/LYZ/C1QA meet) raises set
  B's minimum from 0.043 to **0.0821, above CXCR4's 0.0484, with zero markers
  below it**.
  **CXCR4 is negative on two axes, not four** — circularity retention 0.30× and
  the within-label sign reversal. The conclusion is unchanged (the reversal
  alone is decisive) but "negative on every axis" overstated it.
  Also: D23's haematopoietic quantification is qualified (insensitive to the
  screen, 0.250 → 0.270, but not clean); `HBB` should not have counted as a
  control (no erythroid population in decidua, myeloid CPM 1.9); the PAEP
  downgrade stands, resting on the untouched non-haematopoietic side.
  **This correction moves in CXCR4's favour and is made for that reason** — the
  review had been accumulating evidence against it, and this piece does not
  survive checking.
- **D32.** Two CCR5 statements corrected. The rulers do **not** disagree — both
  point away from pure ambient; `unmeasurable` comes from the power gate
  (myeloid detects CCR5 at 4.0%, below the 5% gate). And **"CCR5 is not
  detectable in dNK" is withdrawn as overstated**: the floor band means the
  readout has no comparative power, not that the gene is absent. The
  defensible version establishes neither expression nor its absence.
- **`TRANSFERABLE.md` §6** added: ambient controls cannot be chosen from data on
  the target's own side of the lineage tree, because separating "expresses" from
  "picks up" is the quantity the estimator measures. Prior exclusivity knowledge
  is required. Wet-lab analogue: isotype/FMO controls are chosen by biology, not
  by which channel was quiet.

## v1.4.4 — 2026-08-13 — scope correction to D31; numeric provenance stamped

- **D33 corrects D31's scope.** What fails is *auto-expanding the control set
  with a statistic*, not the haematopoietic side as such. IGKC/LYZ/C1QA are
  chosen from prior exclusivity knowledge — "NK does not transcribe
  immunoglobulin/lysozyme/C1q" is as hard as "…collagen" — so **ruler A remains
  usable for haematopoietic target genes**, ceiling 0.0920, and XCL1's 0.021
  reading is not withdrawn. XCL1 never depended on it, but a layer of standing
  evidence should not be discarded by over-correction.
- **HBB, secondary.** "No source lineage" voids **ruler B** (confirmed: current
  calibration assigns HBB to Stromal at 8.3 CPM and excludes it below the
  10 CPM gate) but does not automatically void **ruler A** — free haemoglobin
  from never-captured red cells is a different contamination route from RNA
  shed by captured cells, and 0.408 may be reporting it. Low priority.
- **Two quantities stamped with their provenance**, after each appeared at two
  values: CCL3 detection **0.817** (depth-matched, dNK1–3, per-subset mean — the
  band assignment) vs 0.909 (raw, incl. dNKp); ruler B band **0.0538** (v1.2,
  source by total counts — current) vs 0.0441 (v1.1, source by per-cell CPM).
  Neither changes a conclusion. `TRANSFERABLE.md` §8 carries the table.
- **`TRANSFERABLE.md` §7** added: negative controls are chosen by biology in
  every setting. Naming "whichever readout stayed flat" as the negative control
  is the wet-lab form of §6's circularity — the selection statistic is the
  measurement.

## v1.4.5 — 2026-08-13 — "CCL3 does not reproduce" split and half withdrawn

- **D34.** v1.4 A11 assigned CCL3 to the measurable band on a **cross-subset
  mean (0.817)** — a statistic the per-arm ceiling rule does not use. Per arm:
  dNK1 0.800, dNK2 0.777, **dNK3 0.884, above the 0.85 ceiling.**
  The pipeline had already flagged both dNK3 contrasts
  `magnitude_is_lower_bound`; the prose had not. **The claim is split:**
  dNK1-vs-dNK2 (−2.3 pp, both arms measurable) is a powered true negative;
  the dNK3 contrasts (+8.8 and +10.1 pp, saturated dNK3 arm) are **not nulls**
  — direction readable, magnitude a lower bound, "no difference" not
  adjudicable. Not a revival of the rejected floor-band argument: that was the
  opposite end of the scale and was refuted by 0.817.
  **Process fix:** band assignment is reported per arm, matching the rule that
  consumes it.
- **D35.** The `>= 10 CPM` source-lineage power gate independently caught two
  unrelated failures — CCR5 (myeloid source at 9.37 CPM, a genuinely
  low-abundance receptor) and HBB (source assigned to Stromal at 8.3 CPM, a
  gene whose real source lineage is absent from the tissue). One untuned rule,
  two different scenarios.
- **`TRANSFERABLE.md` §7 extended to five settings**, adding Western/qPCR
  loading controls and assay normalisation bases. `GAPDH` is a glycolytic
  enzyme: normalising a glycolysis/OXPHOS experiment to it is this error in its
  purest form — chosen because it looks stable, when its stability is the
  perturbed quantity. Same for choosing cells vs protein vs DNA after seeing
  which looks cleanest.
