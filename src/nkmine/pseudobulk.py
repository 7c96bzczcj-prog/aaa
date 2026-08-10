"""Pseudobulk construction with power balancing (Protocol Phase 2.3-2.6).

The single most important function here is `balanced_pseudobulk`.  NK is
a rare lineage; if pseudobulks are built from whatever cells happen to
be present, NK carries a larger standard error by construction, and
"NK did not change" becomes an unfalsifiable default -- the protocol's
stated failure mode.  Balancing subsamples every lineage down to the
rarest one within each (individual x condition) before summing, so the
four lineages enter Phase 4 on comparable footing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PseudobulkSet:
    """Genes x samples counts plus the sample annotation columns."""

    counts: np.ndarray          # genes x samples
    genes: np.ndarray
    patient: np.ndarray
    condition: np.ndarray
    lineage: np.ndarray
    n_cells: np.ndarray         # cells contributing to each pseudobulk

    def subset_lineage(self, lineage: str) -> "PseudobulkSet":
        m = self.lineage == lineage
        if not m.any():
            raise ValueError(f"no pseudobulk samples for lineage {lineage!r}")
        return PseudobulkSet(
            counts=self.counts[:, m],
            genes=self.genes,
            patient=self.patient[m],
            condition=self.condition[m],
            lineage=self.lineage[m],
            n_cells=self.n_cells[m],
        )

    def complete_pairs(self) -> "PseudobulkSet":
        """Keep only individuals observed in both conditions, per lineage.

        A paired design cannot use a half pair; dropping them here keeps
        `~ patient + condition` full rank.
        """
        keep = np.zeros(len(self.patient), bool)
        for lin in np.unique(self.lineage):
            for pat in np.unique(self.patient):
                m = (self.lineage == lin) & (self.patient == pat)
                if len(np.unique(self.condition[m])) == 2:
                    keep |= m
        return PseudobulkSet(
            counts=self.counts[:, keep],
            genes=self.genes,
            patient=self.patient[keep],
            condition=self.condition[keep],
            lineage=self.lineage[keep],
            n_cells=self.n_cells[keep],
        )


def balanced_pseudobulk(
    cell_counts,
    genes,
    cell_patient: np.ndarray,
    cell_condition: np.ndarray,
    cell_lineage: np.ndarray,
    lineages=("NK", "CD8T", "CD4T", "B", "Myeloid"),
    n_reps: int = 20,
    min_cells: int = 30,
    rng: np.random.Generator | None = None,
) -> PseudobulkSet:
    """Protocol 2.4: equalise cell numbers across lineages, then sum.

    For each (individual x condition), take n_min = the smallest lineage
    cell count, subsample every lineage to n_min, sum to pseudobulk, and
    repeat `n_reps` times taking the elementwise median.  Groups where
    the rarest lineage falls below `min_cells` are dropped entirely --
    partial balancing would reintroduce the very asymmetry this exists
    to remove.

    `cell_counts` may be a dense array or a scipy sparse matrix
    (genes x cells).
    """
    import scipy.sparse as sp

    rng = rng or np.random.default_rng(0)
    genes = np.asarray(genes)
    cell_patient = np.asarray(cell_patient)
    cell_condition = np.asarray(cell_condition)
    cell_lineage = np.asarray(cell_lineage)
    sparse = sp.issparse(cell_counts)
    if sparse:
        cell_counts = cell_counts.tocsc()

    out_counts, out_pat, out_cond, out_lin, out_n = [], [], [], [], []
    dropped = []

    for pat in np.unique(cell_patient):
        for cond in np.unique(cell_condition):
            base = (cell_patient == pat) & (cell_condition == cond)
            if not base.any():
                continue
            idx_by_lin = {l: np.flatnonzero(base & (cell_lineage == l)) for l in lineages}
            sizes = {l: len(v) for l, v in idx_by_lin.items()}
            n_min = min(sizes.values())
            if n_min < min_cells:
                limiting = min(sizes, key=sizes.get)
                dropped.append((pat, cond, dict(sizes), limiting))
                continue

            for l in lineages:
                idx = idx_by_lin[l]
                reps = np.empty((n_reps, len(genes)))
                for r in range(n_reps):
                    pick = rng.choice(idx, size=n_min, replace=False)
                    if sparse:
                        reps[r] = np.asarray(cell_counts[:, pick].sum(axis=1)).ravel()
                    else:
                        reps[r] = cell_counts[:, pick].sum(axis=1)
                out_counts.append(np.median(reps, axis=0))
                out_pat.append(pat)
                out_cond.append(cond)
                out_lin.append(l)
                out_n.append(n_min)

    if not out_counts:
        raise ValueError(
            "no (individual x condition) group met the balancing threshold; "
            f"min_cells={min_cells}. Dropped groups: {dropped[:5]}"
        )

    ps = PseudobulkSet(
        counts=np.column_stack(out_counts),
        genes=genes,
        patient=np.array(out_pat),
        condition=np.array(out_cond),
        lineage=np.array(out_lin),
        n_cells=np.array(out_n),
    )
    ps.dropped_groups = dropped  # attached for the QC report
    return ps


def simple_pseudobulk(
    cell_counts,
    genes,
    cell_patient: np.ndarray,
    cell_condition: np.ndarray,
    cell_lineage: np.ndarray,
    lineages=("NK", "CD8T", "CD4T", "B", "Myeloid"),
) -> PseudobulkSet:
    """Unbalanced aggregation: sum every available cell.

    This is the naive alternative to `balanced_pseudobulk`, kept so the
    cost of skipping Phase 2.4 can be measured rather than asserted.
    A rare lineage keeps its larger standard error here, which is the
    mechanism that turns "underpowered" into a spurious "unchanged".
    """
    import scipy.sparse as sp

    genes = np.asarray(genes)
    sparse = sp.issparse(cell_counts)
    if sparse:
        cell_counts = cell_counts.tocsc()

    cols, pats, conds, lins, ns = [], [], [], [], []
    for pat in np.unique(cell_patient):
        for cond in np.unique(cell_condition):
            for l in lineages:
                idx = np.flatnonzero(
                    (cell_patient == pat) & (cell_condition == cond) & (cell_lineage == l)
                )
                if len(idx) == 0:
                    continue
                sub = cell_counts[:, idx]
                col = np.asarray(sub.sum(axis=1)).ravel() if sparse else sub.sum(axis=1)
                cols.append(col)
                pats.append(pat)
                conds.append(cond)
                lins.append(l)
                ns.append(len(idx))

    return PseudobulkSet(
        counts=np.column_stack(cols),
        genes=genes,
        patient=np.array(pats),
        condition=np.array(conds),
        lineage=np.array(lins),
        n_cells=np.array(ns),
    )


def detection_filter(ps: PseudobulkSet, min_count: int = 10, min_frac: float = 0.5,
                     lineages=("NK", "CD8T", "CD4T", "B", "Myeloid")) -> np.ndarray:
    """Protocol 2.5: keep genes detected in *every* lineage.

    Without this, "NK did not change" is frequently just "NK does not
    express this gene", which is a different statement.
    """
    keep = np.ones(ps.counts.shape[0], bool)
    for l in lineages:
        sub = ps.subset_lineage(l)
        frac = (sub.counts >= min_count).mean(axis=1)
        keep &= frac >= min_frac
    return keep


def saturation_flags(
    cell_counts,
    genes,
    cell_patient: np.ndarray,
    cell_condition: np.ndarray,
    cell_lineage: np.ndarray,
    lineage: str = "NK",
    hi: float = 0.90,
    lo: float = 0.05,
) -> np.ndarray:
    """Protocol 2.6: flag genes with no room to move in `lineage`.

    If detection rate in both conditions is >90% (ceiling) or <5%
    (floor), an apparent "no change" is uninformative and the gene is
    barred from Q4 calls.
    """
    import scipy.sparse as sp

    genes = np.asarray(genes)
    sparse = sp.issparse(cell_counts)
    if sparse:
        cell_counts = cell_counts.tocsc()

    rates = {}
    for cond in np.unique(cell_condition):
        m = np.flatnonzero((cell_lineage == lineage) & (cell_condition == cond))
        if len(m) == 0:
            return np.zeros(len(genes), bool)
        sub = cell_counts[:, m]
        if sparse:
            det = np.asarray((sub > 0).sum(axis=1)).ravel() / len(m)
        else:
            det = (sub > 0).sum(axis=1) / len(m)
        rates[cond] = det

    vals = list(rates.values())
    both_hi = np.all([v > hi for v in vals], axis=0)
    both_lo = np.all([v < lo for v in vals], axis=0)
    return both_hi | both_lo


def se_balance_report(de_by_lineage: dict, ratio_stop: float = 2.0) -> dict:
    """Protocol 2.4 QC / stop rule S3.

    Compares the median moderated SE of NK against the other lineages.
    If NK's is more than `ratio_stop` times the others', Q4 cannot be
    assessed on this dataset and the protocol says to stop.
    """
    med = {l: float(np.median(d.se)) for l, d in de_by_lineage.items()}
    others = [v for l, v in med.items() if l != "NK"]
    ratio = med["NK"] / float(np.median(others))
    return {
        "median_SE": med,
        "NK_vs_others_ratio": ratio,
        "passes": bool(ratio <= ratio_stop),
        "stop_rule": "S3",
        "verdict": (
            "PASS - NK standard errors comparable to other lineages"
            if ratio <= ratio_stop
            else f"STOP (S3) - NK SE is {ratio:.2f}x the other lineages; "
            "'NK unchanged' would be unfalsifiable on this dataset"
        ),
    }
