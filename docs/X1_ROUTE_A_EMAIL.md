# Route A — draft data request to the GSE302113 authors

**Status: DRAFT, not sent.** Sending is the user's call.

**To**: Ansuman T. Satpathy (satpathy@stanford.edu), Caleb A. Lareau
(lareauc@mskcc.org); cc Vincent V. Liu (liuv@stanford.edu)
**Re**: Liu VV et al., *Cancer Cell* 2026;44(7):1509–1521.e4 (PMID 42242233) /
GSE302113

---

## Two asks

1. **The per-donor union heteroplasmy matrices**, or equivalently the Mitotrek
   clone assignments (cell barcode → clone ID) used in the paper.
2. **The numeric source data behind Figures 2E, 2F and S6C** — the clone-sharing
   heatmaps.

## Why each

**(1)** GSE302113 deposits one heteroplasmy matrix per library, each restricted to
that library's own mgatk-selected variants. The Methods describe building "a cell
by variant heteroplasmy matrix combining all samples for each donor", and that
union matrix is what a cross-compartment analysis needs: a variant selected in
tumour and non-involved lung but not in blood has no column in the blood library,
so absence from blood cannot be distinguished from absence of ascertainment.
Restricting to variants selected in all three compartments does not fix this — it
selects variants that are blood-present by construction (0.0% of SU-L-005's 507
tri-compartment variants have zero blood carriers).

**(2)** The heatmap panels include NK/ILC as an annotated cell type, but the
NK entries are not quantified in the text and the only deposited supplementary
file is the patient table. We would like the underlying values rather than
re-deriving numbers that already exist.

## What it is for

We are asking whether intratumoural NK are clonally closer to NK of the matched
non-involved lung than to circulating NK — the NK counterpart of the myeloid
result, which the paper reports for monocytes/macrophages/DC3 but does not break
out for NK.

## Note on the ask

`mitotrek` is already public on GitHub, and mgatk objects are routinely deposited
in this subfield — GSE197037 deposits `*_mgatk.rds.gz` per library, which is
exactly the un-lossy form that would resolve this. So the request is for the
intermediate that the published pipeline already produces.

Happy to work from whatever form is convenient, including a barcode → clone table.

---

## Internal note — do not put this in the email

There is an **open question the archive cannot settle**: whether the paper's own
NK-relevant statement — that innate cells "including … NK cells" show high clone
sharing "suggesting … recent hematopoietic output" — was computed on the union
matrix or on per-library matrices.

- On the **union** matrix, that statement is unbiased and is genuine evidence on
  the blood-origin side.
- On **per-library** matrices, it is inflated by the same ascertainment effect
  that stopped our Tier 1.

From the deposit alone the two are indistinguishable. Their Methods describe the
union procedure, so the former is the more likely reading — but it is not
checkable from what is public, and receiving ask (1) resolves it either way.

**This is not an allegation and must not be phrased as one in correspondence.**
It is a reason the request is worth making, not a claim about their result.

A second, weaker caution, also internal: that sentence aggregates NK with
monocytes, macrophages and DCs. Monocytes are short-lived and blood-derived, so an
aggregate over that group is liable to be carried by the myeloid members, and the
NK cell of the heatmap was never reported on its own. That is a reason to want the
NK number specifically — which is ask (2).
