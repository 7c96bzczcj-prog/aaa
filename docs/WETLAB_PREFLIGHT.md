# Wet-lab pre-flight constraints — write into the protocol BEFORE the first run

**Hard date: experiments start next week.** Every constraint below is free to
impose now and expensive-to-impossible to impose afterwards. That asymmetry is
the only reason this document exists.

**Provenance note.** Constraints 1–6 are the spec author's, recorded here as
given, with the failure mode and the compliance check filled in. Constraint 7 is
theirs too (flagged for the human-NK arm). Nothing here is derived from data —
this environment has no access to that project's workspace.

**How to use it.** Each item has: the rule, the failure it prevents, **the exact
sentence to put in the protocol**, and **how to verify compliance after the fact**.
The last column matters — a constraint nobody can audit later is a wish.

---

## 1. No overnight rest, no IL-15/IL-2 rescue — fresh cells straight on

**Rule.** Cells go onto the assay fresh. No overnight hold, no cytokine rescue
step before readout.

**Failure prevented.** Time in culture and cytokine exposure are each sufficient
to move the readout on their own, so a rested/rescued arm cannot attribute an
effect to the manipulation. The stated basis is the 19% → 66% shift.

**Protocol sentence.**
> Cells are assayed within [X] hours of isolation, without overnight rest and
> without IL-15 or IL-2 supplementation at any point before readout. Where a hold
> is unavoidable, hold time is recorded per sample and reported.

**Compliance check.** Record isolation time and assay start time per sample. If
hold time varies, plot the readout against hold time before interpreting anything
— if it correlates, the experiment is measuring the hold.

**Note.** If any arm *requires* a rescue condition, it is a separate arm with its
own control, not a preprocessing step applied to everything.

---

## 2. Declare floor/ceiling-limited readouts BEFORE the first run

**Rule.** Before running, determine which flow readouts sit near 0% or 100% in the
control condition, and declare those **unmeasurable in advance**.

**Failure prevented.** A readout whose control is at 2% cannot go down and a
readout at 97% cannot go up. Post-hoc, a floored readout looks like "no effect in
this arm" — which is then reported as a negative result rather than as an
instrument that could not have detected one.

**Protocol sentence.**
> Readouts whose control-condition value falls below 5% or above 95% are declared
> **not measurable** for directional change in this design, before data
> collection. Such readouts are reported as floor/ceiling-limited, not as null.

**Compliance check.** The declared list is timestamped and committed before the
first run. Any readout later reported as "unchanged" is cross-checked against that
list.

**How to get the control values without burning the experiment.** Use pilot data,
prior runs, or a single control-only plate. This is the one item that may need a
small pre-run — worth it, because it is what separates "no effect" from "no
sensitivity".

---

## 3. Fix the effect scale before seeing the first batch

**Rule.** Decide now whether the effect is **percentage points**, **fold change**,
or **MFI**, per readout. Written down before any data is seen.

**Failure prevented.** At low baselines these disagree violently: 2% → 6% is +4
points or 3×, and the choice made after seeing the data is the choice that makes
the result look best. This is the most common silent p-hack in flow work.

**Protocol sentence.**
> For each readout, the effect scale (percentage-point difference / fold change /
> MFI ratio) and the summary statistic are fixed in this document prior to data
> collection and are not changed afterwards.

**Compliance check.** One table, per readout, committed before run 1. Any analysis
on a different scale is reported as a secondary, clearly-labelled analysis.

**Recommendation where you have a free choice.** For frequencies with low
baselines, percentage points understate and fold overstates; pre-specify **both**
and nominate one as primary. Pre-specifying both is legitimate; choosing after is
not.

---

## 4. n is animals or donors, never wells

**Rule.** The unit of replication is the animal or the donor. Wells are technical
replicates.

**Failure prevented.** Pseudoreplication — treating 3 wells from 1 mouse as n = 3
inflates significance by roughly the square root of the technical replicate count,
and is invisible in the final figure.

**Protocol sentence.**
> The experimental unit is the individual animal (or human donor). Technical
> replicates are averaged within animal/donor before any statistical test, and
> every reported n is stated with its unit.

**Compliance check.** Every figure legend states "n = X animals" or "n = X
donors". A legend that says only "n = 6" fails the check. Where wells are pooled,
say so.

**Corollary.** Power calculations use the animal/donor n, not the well count.

---

## 5. Negative and loading controls chosen by biology, not by flatness

**Rule.** Controls are selected because the biology says they should not move —
never because they happened to come out flat.

**The specific trap, stated because it is the purest form of it.** **GAPDH is a
glycolytic enzyme.** Using it as a loading control in an experiment that
manipulates the glycolysis/OXPHOS balance means normalising the signal to
something the manipulation is expected to change. **ACTB is not a safe fallback**
— actin moves under metabolic stress and hypoxia, both of which are in play here.

**What to use instead.** For metabolic experiments the correct answer is not a
different housekeeping gene but **total-protein normalisation** — stain-free gels,
REVERT or Ponceau — because it does not assume any gene is invariant. If a
housekeeping protein must be used, validate it *in these conditions* and show the
validation.

**Protocol sentence.**
> Western blot signal is normalised to total protein (stain-free / REVERT).
> GAPDH and ACTB are not used as loading controls in metabolic-manipulation arms.
> Any housekeeping control used elsewhere is validated under the same conditions
> and the validation is shown.

**Compliance check.** Every blot figure shows the total-protein image. Reviewers
ask for this now anyway.

---

## 6. Fix the Seahorse normalisation basis before reading plate 1

**Rule.** Choose cell number, protein, or DNA **before** the first plate is read.

**Failure prevented.** OCR/ECAR are per-well quantities and mean nothing
unnormalised. Choosing the basis after seeing the data is choosing the result.

**Which basis, and why it is not a free choice.** If the manipulation changes cell
size, granularity or adherence — and metabolic manipulations routinely do —
**cell-count normalisation is biased**, because equal counts no longer mean equal
biomass. **DNA-based normalisation (e.g. CyQUANT/Hoechst) is more robust to size
change**; protein sits in between and is itself perturbed if the treatment changes
translation.

**Protocol sentence.**
> Seahorse OCR/ECAR are normalised to [chosen basis], fixed prior to the first
> plate read. Post-run normalisation is performed on every plate, including
> failed plates, and the raw per-well values are retained.

**Compliance check.** Normalisation data exist for every plate, including ones
later excluded. Retain raw values so the basis can be re-derived if challenged.

**Also fix now**: seeding density, the injection strategy and concentrations, and
the exclusion rule for wells (e.g. negative OCR, failed injections) — all of which
are equally post-hoc-tunable.

---

## 7. Human-NK arm only — record NK education status

**Rule.** If any arm uses human NK cells, record **KIR / HLA-C genotype per
donor**, or at minimum surface **NKG2A and pan-KIR** staining.

**Failure prevented.** NK education/licensing sets baseline responsiveness:
self-KIR⁺ educated NK respond differently from uneducated NK. In a study of
metabolic polarisation capacity, education is a live, donor-varying, **currently
uncontrolled** variable — so between-donor variance that is really education gets
attributed to the manipulation, or masks it.

**Minimum viable version.** Surface NKG2A + pan-KIR staining on the same panel
costs two channels and no extra donors. Full HLA-C C1/C2 genotyping is better and
can be done retrospectively **if you bank DNA now** — which is the part that must
happen before the run.

**Protocol sentence.**
> For each human donor, NKG2A and pan-KIR surface expression are recorded on the
> phenotyping panel, and a DNA aliquot is banked for retrospective HLA-C (C1/C2)
> and KIR genotyping.

**Compliance check.** A per-donor education table exists before analysis. Donor
identity is carried through to the analysis stage so the readout can be
stratified.

---

## Pre-flight sign-off

Everything below is committed **before run 1**, in one file, timestamped:

- [ ] Assay-start-within-X-hours rule, and hold time recorded per sample
- [ ] List of floor/ceiling-limited readouts, declared unmeasurable
- [ ] Effect scale + primary summary statistic, per readout
- [ ] Statement that n = animals/donors, with technical replicates averaged within
- [ ] Loading-control policy (total protein; GAPDH/ACTB barred from metabolic arms)
- [ ] Seahorse normalisation basis, seeding density, injections, well-exclusion rule
- [ ] Human arm: NKG2A/pan-KIR on panel + DNA banked per donor

**The point of the sign-off is auditability.** Each line has a compliance check
that a third party can run on the finished figures. A constraint that cannot be
checked after the fact does not constrain anything.
