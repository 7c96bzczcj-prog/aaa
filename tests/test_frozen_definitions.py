"""Guard for the frozen round-1 definitions file.

The file must not change after it is committed. The hash check catches any
edit, including whitespace. The value checks restate the task
description independently, so the file and this test would both have to
be edited to change a gene list or threshold. A deliberate change means a
new versioned file, not an edit to this one.
"""

from __future__ import annotations

import hashlib
import tomllib
from pathlib import Path

DEF = Path(__file__).resolve().parents[1] / "definitions" / "nk_atlas_round1_frozen.toml"
SHA256 = "0c544303b5590c9586cad5313052447a614067a6567ec7df7c9d598f476d1d36"


def _load():
    with DEF.open("rb") as fh:
        return tomllib.load(fh)


def test_file_is_byte_identical_to_the_frozen_commit():
    assert hashlib.sha256(DEF.read_bytes()).hexdigest() == SHA256


def test_nk_gate_matches_task_spec():
    g = _load()["nk_gate"]
    assert g["keep_require_any"] == ["NCAM1", "KLRD1"]
    assert g["keep_require_all"] == ["NCR1"]
    assert g["drop_if_any"] == ["CD3D", "CD3E", "IL7R", "RORA", "GATA3"]
    # The ILC1 exclusion genes must be present and mandatory.
    assert set(g["ilc1_exclusion_genes"]) == {"IL7R", "RORA", "GATA3"}
    assert set(g["ilc1_exclusion_genes"]) <= set(g["drop_if_any"])
    assert g["ilc1_exclusion_mandatory"] is True


def test_signatures_match_task_spec_and_are_disjoint():
    d = _load()
    s, single = d["signatures"], d["single_genes"]
    assert s["loss"] == ["PRF1", "GZMB", "GZMH", "GNLY", "FCGR3A",
                         "KLRK1", "CD226", "NCR3", "CX3CR1", "TBX21"]
    assert s["gain"] == ["CD9", "ITGA1", "VEGFA", "CXCR4", "PAEP",
                         "ITGAE", "CD69", "EOMES"]
    assert single["recent_degranulation"] == ["SELL", "LAMP1"]
    assert single["tissue_origin"] == ["CXCR6", "S1PR5"]
    singles = set(single["recent_degranulation"]) | set(single["tissue_origin"])
    assert not singles & (set(s["loss"]) | set(s["gain"]))
    assert not set(s["loss"]) & set(s["gain"])


def test_call_rule_and_references_match_task_spec():
    d = _load()
    r = d["call_rule"]
    assert (r["signature"], r["min_genes"], r["padj_threshold"]) == ("gain", 3, 0.05)
    assert r["padj_strict_less_than"] is True
    assert r["require_consistent_fold_change_direction"] is True
    assert d["reference_groups"] == {
        "primary": "healthy_donor_peripheral_blood_NK_all",
        "independent": "healthy_donor_CD56bright_NK",
    }
    assert all(v is False for v in d["prohibitions"].values())
