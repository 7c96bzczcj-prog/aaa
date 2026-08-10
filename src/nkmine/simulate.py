"""Synthetic multi-lineage data with planted quadrant ground truth.

Used to answer the question the protocol says must be answered before
any real result is looked at (Phase 5.3): does this pipeline actually
have the power to find a Q4 gene, and -- more importantly -- does it
refuse to invent one when a lineage is merely underpowered?
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

LINEAGES = ("NK", "CD8T", "CD4T", "B", "Myeloid")


@dataclass
class SimData:
    counts: np.ndarray          # genes x cells
    genes: np.ndarray
    cell_patient: np.ndarray
    cell_condition: np.ndarray
    cell_lineage: np.ndarray
    truth: dict = field(default_factory=dict)   # gene -> planted quadrant
    planted_lfc: dict = field(default_factory=dict)  # gene -> {lineage: lfc}


def simulate(
    n_patients: int = 12,
    n_genes_per_class: int = 30,
    n_null_genes: int = 200,
    cells_per_lineage: dict | None = None,
    effect: float = 1.0,
    dispersion: float = 0.3,
    patient_sd: float = 0.35,
    baseline_mean: float = 3.0,
    rng: np.random.Generator | None = None,
) -> SimData:
    """Generate a cell-level count matrix with known quadrant membership.

    `cells_per_lineage` controls the power imbalance being studied; the
    default gives NK far fewer cells than the others, as in real CD45+
    tumour data.
    """
    rng = rng or np.random.default_rng(0)
    if cells_per_lineage is None:
        cells_per_lineage = {"NK": 60, "CD8T": 300, "CD4T": 300, "B": 200, "Myeloid": 300}

    classes = {
        "Q1": {l: (effect if l == "NK" else 0.0) for l in LINEAGES},
        "Q2": {l: (effect if l in ("NK", "CD8T") else 0.0) for l in LINEAGES},
        "Q3": {l: effect for l in LINEAGES},
        "Q4": {l: (0.0 if l == "NK" else effect) for l in LINEAGES},
    }

    genes, truth, planted = [], {}, {}
    for cls, spec in classes.items():
        for i in range(n_genes_per_class):
            g = f"{cls}_{i:03d}"
            genes.append(g)
            truth[g] = cls
            # Alternate the sign of the planted effect.  If every planted
            # effect pointed the same way, the induced compositional shift
            # would show up as a uniform negative offset on the true-null
            # genes after library-size normalisation -- an artefact of the
            # simulation, not of the method.
            sign = 1.0 if i % 2 == 0 else -1.0
            planted[g] = {l: sign * v for l, v in spec.items()}
    for i in range(n_null_genes):
        g = f"NULL_{i:03d}"
        genes.append(g)
        truth[g] = "unclassified"
        planted[g] = {l: 0.0 for l in LINEAGES}

    genes = np.array(genes)
    n_genes = len(genes)

    # gene baseline abundance, on a log2 scale, shared across lineages
    base = rng.normal(baseline_mean, 1.2, n_genes)
    # each lineage has its own expression offset per gene
    lin_offset = {l: rng.normal(0, 0.5, n_genes) for l in LINEAGES}
    # individual-level random effects
    pat_effect = {
        (p, l): rng.normal(0, patient_sd, n_genes)
        for p in range(n_patients)
        for l in LINEAGES
    }

    cols, c_pat, c_cond, c_lin = [], [], [], []
    for p in range(n_patients):
        for cond in ("Normal", "Tumor"):
            for l in LINEAGES:
                n_cells = cells_per_lineage[l]
                lfc = np.array([planted[g][l] for g in genes])
                mu_log = base + lin_offset[l] + pat_effect[(p, l)]
                if cond == "Tumor":
                    mu_log = mu_log + lfc
                mu = np.exp2(mu_log)
                # negative binomial per cell
                size = 1.0 / dispersion
                prob = size / (size + mu[:, None])
                block = rng.negative_binomial(size, prob, size=(n_genes, n_cells))
                cols.append(block)
                c_pat.extend([f"P{p:02d}"] * n_cells)
                c_cond.extend([cond] * n_cells)
                c_lin.extend([l] * n_cells)

    return SimData(
        counts=np.concatenate(cols, axis=1),
        genes=genes,
        cell_patient=np.array(c_pat),
        cell_condition=np.array(c_cond),
        cell_lineage=np.array(c_lin),
        truth=truth,
        planted_lfc=planted,
    )
