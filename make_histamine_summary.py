"""
Generate confirmed_histaminergic_summary.csv from cross-dataset results.

Combines FAFB entropy with MCNS histamine confirmation for the four
canonical photoreceptor / histaminergic types including R1-6.
"""
from __future__ import annotations

import pandas as pd

from mcns_matching import resolve_fafb_to_mcns
from paths import (
    CONFIRMED_HISTAMINERGIC,
    CONFIRMED_HISTAMINERGIC_SUMMARY,
    ENTROPY_RAW,
    FAFB_MERGED,
    MCNS_MERGED,
    ensure_output_dirs,
)

ensure_output_dirs()

TARGET_TYPES = ["R7", "R8", "Lai", "R1-6"]

entropy = pd.read_csv(ENTROPY_RAW)
fafb = pd.read_csv(FAFB_MERGED)
mcns = pd.read_csv(MCNS_MERGED).dropna(subset=["primary_type", "nt_type"])
mcns_names = mcns["primary_type"].unique().tolist()

rows = []
for cell_type in TARGET_TYPES:
    ent_row = entropy[entropy["cell_type"] == cell_type]
    if ent_row.empty:
        continue
    ent_row = ent_row.iloc[0]

    n_labeled = int((fafb["primary_type"] == cell_type).sum())
    n_total = n_labeled  # same merged table; extend if unlabeled counts needed later

    matched, _ = resolve_fafb_to_mcns(cell_type, mcns_names)
    mcns_hist_frac = None
    if matched:
        sub = mcns[mcns["primary_type"].isin(matched)]
        if len(sub):
            mcns_hist_frac = (sub["nt_type"] == "HIST").mean()

    rows.append({
        "cell_type": cell_type,
        "n_labeled": n_labeled,
        "n_total": n_total,
        "entropy": round(ent_row["entropy"], 3),
        "dominant_nt": ent_row["dominant_nt"],
        "dominant_frac": round(ent_row["dominant_frac"], 3),
        "mcns_match": ",".join(matched) if matched else None,
        "mcns_hist_frac": round(mcns_hist_frac, 3) if mcns_hist_frac is not None else None,
    })

summary = pd.DataFrame(rows)
summary.to_csv(CONFIRMED_HISTAMINERGIC_SUMMARY, index=False)
print(f"Saved {CONFIRMED_HISTAMINERGIC_SUMMARY}")
print(summary.to_string(index=False))

# Also refresh the detailed file if we have all four types
if len(summary) >= 3:
    detail = entropy[entropy["cell_type"].isin(summary["cell_type"])]
    detail.to_csv(CONFIRMED_HISTAMINERGIC, index=False)
