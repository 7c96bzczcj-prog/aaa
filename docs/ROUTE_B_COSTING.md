# Route B re-costed on measured numbers, not on a guess

**Why re-costed.** The earlier "TB-scale, days of compute" figure was an estimate
made without measuring anything. These are measured.

**Caveat on "your platform".** This environment is a **4-core container with 30 GB
free disk** — the same class of machine the original estimate came from, and it
cannot host 1.12 TB. So the download figures below are *measured throughput*
extrapolated to full volume; the compute figures are *derived from published
cellranger-atac scaling*, not run here. Anything to be run on your cluster is
labelled as such.

## Measured volume (SRA, BioProject PRJNA1288771)

| quantity | value |
|---|---|
| runs | **38** (all ATAC-seq, NovaSeq 6000) |
| SRA compressed | **1.12 TB** (1,121,977 MB) |
| total bases | **3.52 Tbp** |
| total spots | **11.13 G** |
| spots per run | median **278 M**, range 151–516 M |
| median run size | 28.3 GB |
| FASTQ (gz) after `fasterq-dump`, est. 2–3× | **2.2–3.4 TB** |

**Disk is the first constraint, and it is not the compute.** You need ~1.1 TB for
SRA plus 2.2–3.4 TB for FASTQ plus BAM/outs (cellranger-atac writes a
position-sorted BAM comparable to the FASTQ). Budget **6–8 TB of scratch** if
converting all 38 up front, or stage per-library and delete as you go, which caps
it at ~500 GB.

## Measured download throughput

Single stream from `sra-pub-run-odp.s3.amazonaws.com`, measured here:
**53.5 MB/s** (524 MB in 9.8 s).

| bandwidth | full 1.12 TB |
|---|---|
| 53.5 MB/s (measured, 1 stream) | **5.8 h** |
| ~150 MB/s (4–8 parallel streams) | **2.1 h** |
| ~400 MB/s (well-provisioned 10GbE) | **0.8 h** |

**Download is not the bottleneck.** Under a day on any reasonable link, and
parallelisable across runs with no coordination.

## Compute — derived, not measured

`cellranger-atac count` at ~278 M read pairs per library. Published guidance is
8+ cores / 64 GB, with runtime roughly linear in reads; ~250–300 M read pairs
lands around **6–12 h on 16 cores** per library.

| configuration | 38 libraries |
|---|---|
| 1 × 16-core node, serial | **12–19 days** |
| 4 concurrent nodes | **3–5 days** |
| 8 concurrent nodes | **1.5–2.5 days** |

`mgatk` afterwards is comparatively cheap — it reads the chrM slice of each BAM;
hours, not days, and trivially parallel. Mitotrek on top is minutes.

## Verdict against the stated threshold

The threshold was: **if it is two or three days, do not wait for the email.**

- **On a cluster (≥4 concurrent 16-core nodes): 3–5 days end to end.** That is at
  or just past the threshold — close enough that it should start **in parallel
  with** the email, not instead of it. The email costs nothing to have sent, and
  if the union matrices arrive they supersede the rebuild.
- **On a single node: 2–3 weeks.** Past the threshold; wait for the email, or
  rebuild only the 4 admissible NSCLC donors.

**Cheapest useful variant.** Only 4 donors are admissible for X1 (SU-L-001, -002,
-004, -005), which is **20 of the 38 libraries** — roughly half the download and
half the compute. On 4 nodes that is **~1.5–2.5 days**, which is inside the
threshold. If Route B is started, start it on those 20 libraries only.

**One thing the rebuild also buys.** It produces the mgatk objects for these
libraries, i.e. the un-lossy form the archive lacks — so the ascertainment
confound disappears and Tier 1 as originally preregistered becomes runnable, not
just Route E.
