# Phase 0 — prior-art scan

Date: 2026-08-10. Searches run against the open web. Recorded whether
each returned a hit that would trigger stop rule S1 ("someone has
already systematically reported a Q4-type result").

## Result: S1 does NOT fire

No publication was found that defines, enumerates, or systematically
reports the class "genes that change across multiple immune lineages in
tumour versus adjacent normal, while NK cells do not change". The
negative space is still negative.

Nor was any work found that formalises cross-cell-type effect
sharing/specificity using **equivalence testing**. This appears to be
genuinely unoccupied methodological ground, and is the strongest
novelty claim available to this project — stronger than the biology.

## Query battery and outcomes

| Query theme | Outcome |
|---|---|
| `"cell-type-specific" / "lineage-specific"` + `resistant/refractory/unaffected` + TME + single-cell | No Q4-type result. Returns NK-intrinsic state papers (pan-cancer NK atlases, NK dysfunction) — all Q1-framed. |
| `"shared response"` + cell type + tumour + scRNA-seq + NK | Returns pan-cancer **CD8 T** programme papers, and NK↔CD8 interaction papers. Shared programmes are described *within* a lineage across patients, never *across* lineages within a patient. |
| `"does not change"/"remains unchanged"/"is preserved"` + NK + tumour-infiltrating | Nothing. Absence is not reported as a finding, which is the whole premise. |
| `"differential response"/"divergent response"` + lymphocyte + `"same microenvironment"` | Returns within-lineage heterogeneity (TCR-T subsets responding differently), not between-lineage. |
| `mashr` / condition-specific effects + cell type + tumour | mashr and muscat handle sharing/specificity by empirical-Bayes shrinkage and celltype×condition interaction terms. Neither *accepts* a null; neither is applied to the NK-versus-other-lineage question. |
| `equivalence testing` / `TOST` / `interval null` + differential expression + cell type | **No hits combining these.** Nearest neighbour is a 2025 *Briefings in Bioinformatics* paper on nested-design scRNA DE, which is about pseudoreplication, not about accepting a null. |

## Adjacent occupants — must be engaged in any write-up

- **OSCA / standard practice.** Genes significant across most cell-type
  comparisons are treated as contamination and discarded. The field's
  default is that Q3 is an artefact to be removed, so Q3 has never been
  examined as biology and Q4 has no name because its complement was
  thrown away.
- **Sturm et al., *Bioinformatics* 2019.** Formalised signature spillover;
  explicitly tested NK-versus-CD8. Deconvolution benchmarks routinely
  *exclude* NK when deriving T-cell signatures because the two are
  transcriptionally confusable — the same confusability this protocol
  must survive.
- **Bruni et al., *Front Immunol* 2022.** Argues qualitatively that
  published NK signatures are not NK-specific.
- **Marsh et al., *Nat Neurosci* 2022.** A conserved immediate-early /
  heat-shock dissociation geneset shared across cell types and species.
  This is the empirical basis for the Q3 positive control.
- **mashr** (Urbut/Stephens, *Nat Genet* 2019) and **muscat**
  (*Nat Commun* 2020). Sharing/specificity decomposition and
  celltype×condition interactions are standard; the gap is that neither
  supports an evidenced claim of *no effect*.

## A live example of why this matters

The pan-cancer single-cell NK atlas (*Cell* 2023, 716 patients, 24
cancer types) reports tumour-infiltrating NK cells as characterised by
high **DNAJB1, HSPA1A, FOS, JUN**. That gene set is, essentially
exactly, the Marsh dissociation-stress signature — a textbook Q3
candidate being read as NK tumour biology in a top-tier paper.

This is the clearest available argument for the protocol: nobody is
running the cross-lineage control that would separate the two
interpretations, so a Q3 artefact and a Q1 finding are currently
indistinguishable in practice. It also sharpens the deliverable — even
if Q4 turns up empty, demonstrating that a widely-cited NK signature is
Q3 rather than Q1 is a publishable negative result in its own right.

## Sources

- [Pan-cancer single-cell panorama of human NK cells, *Cell* 2023](https://www.cell.com/cell/fulltext/S0092-8674(23)00849-8)
- [Pan-cancer NK profiling by reference mapping, *Nat Immunol* 2024](https://www.nature.com/articles/s41590-024-01884-z)
- [Marsh et al., dissociation artefacts, *Nat Neurosci* 2022](https://www.nature.com/articles/s41593-022-01022-8)
- [Benchmarking deconvolution pipelines, *Nat Commun* 2020](https://www.nature.com/articles/s41467-020-19015-1)
- [Flexible statistical methods for multi-condition genomics (mashr), *Nat Genet* 2019](https://www.nature.com/articles/s41588-018-0268-8)
- [Single-cell DE between conditions in nested settings, *Brief Bioinform* 2025](https://academic.oup.com/bib/article/doi/10.1093/bib/bbaf397/8232550)
- [Equivalence testing and interval hypotheses (Lakens)](https://lakens.github.io/statistical_inferences/09-equivalencetest.html)
