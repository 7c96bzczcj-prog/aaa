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

## v1.1 — 2026-08-15 — bidirectional load-bearing

**Amendment**, authorised by the spec's author after the v1.0 run.

**Reason.** The v1.0 `load_bearing` definition asked only what presupposes a
claim is TRUE. A claim is equally load-bearing when the field presupposes its
NEGATION: work built on the claim being false must be re-read if it turns out
true. The v1.0 definition was written too narrowly; this is a defect in the
criterion, not in its execution.

**Change.** `load_bearing = HIGH` if EITHER a positive dependency OR a negation
dependency exists (`docs/PREREGISTRATION.md` §A1.2).

**Scope.** Re-applied to **X1 and S1 only**, deliberately and on instruction.
All other rows retain their v1.0 grade and are marked
`load_bearing_basis = v1.0_positive_only`.

**Effect.** X1 MEDIUM → HIGH, and X1 enters the shortlist (it already passed
criteria 2, 3 and 4). S1 MEDIUM → HIGH with no change to its status, because it
fails criterion 3 on both decidability gates independently of load-bearing.

**Not changed.** The four shortlist criteria, the evidence taxonomy, the
decidability gates, and every Phase C threshold including the 20× chrM floor —
which is separately known to be uncalibrated (DEC-14) and whose recalibration is
deferred to a future round so that it is not adjusted alongside a result it
would affect.

**The v1.0 result remains on the record** in `docs/RESULTS.md` next to the v1.1
result; it is not overwritten.

## v1.1 — second pass, 2026-08-15

- **T1 re-graded MEDIUM → HIGH** under the bidirectional criterion, carried by the
  negation side (Sojka 2014 PMID 24714492, whose title asserts T1's negation;
  Klose 2014 PMID 24725403; Nixon 2022 PMID 35394814). **Status unchanged** — T1
  fails criterion 4 independently.
- **Verified mechanically that X1 was the only row a load-bearing re-grade could
  ever have moved.** All 13 other sub-HIGH rows are blocked by an independent
  criterion. Encoded as a check in `scripts/verify_results.py`.
- **T3, T5, T6, R1, R3 opened but not graded** — no negation dependency could be
  established from primary sources, so under G3/G7 their v1.0 grades stand.
  Declared in `phaseA/load_bearing_v11.json` under `_unfinished`.
- `docs/X1_DESIGN.md` added, deliberately **outside** the screening record
  (spec §0 fixed "no new directions" as a non-goal of the round).
