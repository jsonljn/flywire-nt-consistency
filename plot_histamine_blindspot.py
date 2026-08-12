"""
Regenerate the histamine blindspot figure from saved entropy results.

Shows the full FAFB entropy distribution with confirmed-histaminergic cell types
(R7, R8, Lai, R1-6) highlighted — the visual anchor for the headline finding.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paths import (
    CONFIRMED_HISTAMINERGIC,
    CONFIRMED_HISTAMINERGIC_SUMMARY,
    ENTROPY_RAW,
    FIGURES,
    ensure_output_dirs,
)

ensure_output_dirs()

# Primary n>=20 screen
entropy = pd.read_csv(ENTROPY_RAW)

# Confirmed histaminergic types — prefer summary table, else dedicated results file
if CONFIRMED_HISTAMINERGIC_SUMMARY.exists():
    confirmed = pd.read_csv(CONFIRMED_HISTAMINERGIC_SUMMARY)
    highlight_types = confirmed["cell_type"].tolist()
    highlight_data = confirmed.set_index("cell_type")["entropy"].to_dict()
elif CONFIRMED_HISTAMINERGIC.exists():
    confirmed = pd.read_csv(CONFIRMED_HISTAMINERGIC)
    highlight_types = confirmed["cell_type"].tolist()
    highlight_data = confirmed.set_index("cell_type")["entropy"].to_dict()
else:
    highlight_types = ["R7", "R8", "Lai", "R1-6"]
    highlight_data = entropy.set_index("cell_type")["entropy"].reindex(highlight_types).dropna().to_dict()

all_entropies = entropy["entropy"].values
mean_all = all_entropies.mean()
mean_confirmed = np.mean(list(highlight_data.values()))

fig, ax = plt.subplots(figsize=(10, 6))

# Background distribution
ax.hist(
    all_entropies,
    bins=40,
    color="lightgray",
    edgecolor="white",
    alpha=0.9,
    label=f"All FAFB types (n={len(entropy)}, mean={mean_all:.2f})",
)

# Highlight confirmed histaminergic types
colors = {"R7": "#d62728", "R8": "#ff7f0e", "Lai": "#9467bd", "R1-6": "#2ca02c"}
for ct, ent in highlight_data.items():
    pct = (entropy["entropy"] < ent).mean() * 100
    ax.axvline(
        ent,
        color=colors.get(ct, "black"),
        linewidth=2,
        alpha=0.85,
        label=f"{ct} (H={ent:.2f}, ~{pct:.0f}th pct)",
    )

ax.set_xlabel("Shannon entropy of NT predictions (bits)")
ax.set_ylabel("Number of cell types")
ax.set_title(
    "FAFB histamine blind spot\n"
    "Confirmed histaminergic types sit in the extreme right tail; FAFB predicts zero HIST"
)
ax.legend(fontsize=8, loc="upper right")
plt.tight_layout()

out = FIGURES / "histamine_blindspot.png"
plt.savefig(out, dpi=150)
plt.close()
print(f"Saved {out}")
print(f"Mean entropy — all types: {mean_all:.3f}, confirmed histaminergic: {mean_confirmed:.3f}")
