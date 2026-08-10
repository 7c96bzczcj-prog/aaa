"""Shared helpers for the phase driver scripts."""
import pandas as pd


def hto_map(annots: pd.DataFrame, batch: int) -> dict:
    """HTO tag name -> tissue, for one library."""
    out = {}
    for _, r in annots[annots.amp_batch_ID == batch].iterrows():
        if pd.isna(r.HTO):
            continue
        for tag in str(r.HTO).split(","):
            tag = tag.strip()
            if tag:
                out[f"HTO_{tag}"] = r.tissue
    return out
