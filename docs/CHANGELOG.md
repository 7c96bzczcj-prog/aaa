# CHANGELOG — NK lineage-assertion evidence screen

## v1.0 — 2026-08-14

- `docs/PREREGISTRATION.md` created and committed **before any
  literature retrieval**. Contains spec §5 admission criteria, the §2
  evidence taxonomy, the §3 decidability gates, the §4/C.2 power and
  chrM-coverage thresholds, and the §6 hard rules.
- No amendments. Any amendment to the preregistration requires a v1.1
  entry here carrying the reason and a timestamp.

## v1.0 — run complete, 2026-08-14

Preregistration **unchanged** throughout. No amendment was made, so no v1.1 was
needed. The §5 decision rule applied is byte-identical to the one committed at
`9692011` before the first search.

- Phase A: 25 audits + 25 independent adversarial re-examinations (48 agents,
  0 errors). 9 of 25 classifications changed; 3 more re-examined without change.
- Phase B: 17 `INF_*` claims × (decidability + occupancy) = 34 agents, 0 errors.
- Phase C: 38 GSE302113 mtscATAC libraries downloaded and parsed; NK typing
  validated against the dataset's own CD56+/CD56- sorted control.
- Outputs: `out/T1_evidence_audit.tsv`, `T2_decidability.tsv`, `T3_datasets.tsv`,
  `T4_power.tsv`, `T4_power_GSE302113_libraries.tsv`, `T5_shortlist.tsv`;
  `docs/RESULTS.md`; primary-source texts in `refs/`.
- `scripts/verify_results.py` re-derives every quantitative claim in RESULTS.md
  from the result files. It caught three wrong counts in the first draft
  (adversarial changes 10→9, fulltext TRUE 5→10, PARTIAL 19→14); those are
  corrected and all checks now pass.
