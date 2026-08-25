# Transferable failure modes

Appended only when a run produces a new way for an inference to be decided by
a nuisance parameter instead of by the effect.  Each entry states the failure,
the evidence that caught it, how to detect it elsewhere, and the repair.

---

## T1 — A purity gate whose stringency is set by the nuisance parameter it is meant to exclude

*Found: BM-PB-NK preflight, task A, HCA Census of Immune Cells.*

**The failure.** The inherited protocol-2.1 rule — *an NK call requires zero
CD3D / CD3E / CD3G / TRAC / TRBC2 counts* — is a filter on sequencing depth
times ambient composition, not on identity. A cell needs only one stray TCR
transcript to be discarded, and the chance of one stray transcript rises with
both the depth of the cell and the abundance of TCR transcripts in that
library's soup.

**The evidence.** Applied inside clusters already called NK:

| compartment | cells in NK clusters | kept by the zero rule | dropped |
|---|---|---|---|
| bone marrow (MantonBM1) | 1,315 | 703 | **46.5%** |
| blood (MantonBL1) | 5,125 | 1,347 | **73.7%** |

Blood soup is T-dominated and blood NK are sequenced deeper, so the same rule
is roughly 1.6x more aggressive in one arm of the contrast than the other. The
resulting NK fraction was 1.6% of blood mononuclear cells, against the 5–15%
that flow cytometry gives — the rule was not tightening the gate, it was
deleting the population.

**Why it matters beyond this run.** Task A exists to measure how much of a
BM-vs-PB NK difference is compartment ambient RNA. A gate whose stringency is
itself a function of compartment ambient RNA launders that effect into the cell
selection, where it can no longer be seen. The inference then reports a number
whose value is set by the nuisance parameter, which is the definition of the
class.

**A second, independent instance in the same preflight.** Task C's first gate
definition required zero counts across a long exclusion list that included
LYZ, HBB and MS4A1 — three of the most abundant transcripts in marrow soup. It
returned 1–9 mature NK cells in CD45+ marrow libraries that clustering scores
in the hundreds. Same root cause, different task, different data.

**How to detect it.** For any rule of the form "requires zero counts of gene
set G", report the rule's drop rate **separately in each arm of the contrast**
before using it. If the drop rates differ, the rule is part of the contrast.
Also check G against the ambient profile: any gene in the top of the soup
profile makes the rule a depth filter.

**The repair.** Replace the absolute zero with a depth-normalised comparison:
keep a cell if its own lineage panel outscores the competing panel on
log1p(CP10K). A ratio of two quantities measured in the same cell cancels
depth; a zero-count test does not. Here that moved blood NK from 1.6% to 11.3%
of mononuclear cells and marrow NK to 2.8%, both in the range flow cytometry
reports, and it made the drop rate almost equal in the two compartments.

**Cost of keeping the old rule as a labelled arm.** None; it is cheap to
compute both. It is reported as a sensitivity arm rather than deleted, because
the size of the gap between the two arms is itself the diagnostic.

---

## T2 — A per-gene ambient statistic on a low-count gene is a doublet detector, not an ambient estimate

*Found: BM-PB-NK preflight, task A, acceptance test.*

**The failure.** The acceptance test for the soup model was "genes NK cannot
express must come back with an ambient share near 1". It failed in bone marrow
(MS4A1 median share 0.24, CD79A 0.27) and in blood (HBB 0.14). The natural
reading — *the ambient model is wrong, discard the run* — was wrong.

**The evidence.** On MantonBM1 lane 1 the NK pseudobulk carried 8 MS4A1 counts,
of which **5 sat in a single cell**, and 42 LYZ counts, of which **35 sat in a
single cell**. Those are doublets. Meanwhile the per-gene ambient ratios
r_g = observed_g / (total x soup_g) across a 44-gene foreign panel concentrated
tightly around rho: median 0.0123 in marrow and 0.0145 in blood, against a rho
of 0.0130 and 0.0119. The soup model was fine; the test statistic was not.

**Why it matters beyond this run.** A gene with a dozen counts in a pseudobulk
has an ambient share whose variance is dominated by whether one contaminating
cell landed in the gate. Reading that share as an ambient measurement inverts
the meaning of the test: it fails hardest exactly where doublets are most
likely, which is the tissue with the most diverse cell types — bone marrow —
and thereby produces a spurious compartment difference.

**How to detect it.** Before reading a per-gene ambient share, require the
gene's *predicted* ambient counts to clear a floor (10 here). Then check the
per-cell maximum: if one cell holds most of the gene's counts in the gate, the
number is a doublet, not a measurement.

**The repair.** Two changes, both applied identically to every arm:
(1) remove predicted doublets before forming any gate;
(2) judge acceptance on the *spread* of r_g around rho over foreign genes that
clear the predicted-count floor, instead of on a handful of hand-picked genes.

**A related trap the same check caught.** SPINK2 was in the first foreign panel
as a progenitor marker. Its r_g came back at 1.11 in marrow and 1.69 in blood —
NK carry more SPINK2 than the entire pseudobulk could hold as soup, so it is
not NK-foreign and would have inflated rho. A foreign panel needs the same
r_g > 1 sanity check applied to every member.

---

## T3 — A detection-rate admission criterion is uninterpretable without its ambient floor

*Found: BM-PB-NK preflight, admission criterion 2, GSE120221.*

**The failure.** Criterion 2 admits a gene when its detection rate in the
target population falls inside 5%–85%. A detection rate counts cells with at
least one transcript, and in a library with a heavy soup a large share of those
transcripts are ambient. The criterion therefore has a floor it never
measures, and in a soupy dataset the floor alone can lift a gene over the 5%
threshold.

**The evidence.** GSE120221 (Oetjen, 20 healthy marrow donors, no unfiltered
matrix deposited, median 6,175 UMI per cell):

| gene | detected in ALL cells, every lineage included |
|---|---|
| HBB | **100%** |
| GNLY | **52%** |
| NKG7 | **51%** |
| LYZ | 69% |
| MPO | 27% |

Erythroid and myeloid cells do not transcribe GNLY or NKG7. Half of every
lineage carrying them is the soup, not expression. Every cell in every
lymphoid sub-cluster of donor A carried TCR transcripts at 2.6–3.1 per 1,000
UMI — the cytotoxic sub-cluster included — so the NK gate could not be formed
there either.

**Why it matters beyond this run.** The floor is a property of the library, not
of the gene, so it moves between datasets and between compartments of the same
dataset. Two datasets can return the same detection rate for the same gene and
mean entirely different things by it. Any admission, filtering or "expressed
in X% of cells" statement built on a raw detection rate inherits that.

**The repair, stated as a rule.**

> Criterion 2 (detection rate inside 5%–85%) is not a criterion on its own.
> It holds only when reported together with an **ambient floor**: the same
> gene's detection rate in a population that certainly does not express it,
> measured in the same libraries. The admissible quantity is the gap between
> the two, not the rate. Where no unfiltered matrix exists, the floor cannot
> be estimated and the criterion cannot be applied at all.

Past conclusions that used the bare criterion need re-examination on this
point. This preflight reports the floor for every detection rate it quotes
(`results/ltnk_detection_gate.tsv`, columns `floor_B_*` and `floor_Erythroid_*`).
