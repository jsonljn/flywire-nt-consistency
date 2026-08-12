"""
Summary visualization of the three systematic confusion patterns
(histamine blind spot, ORN serotonin confusion, Dm glutamate confusion).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paths import FIGURES, THREE_PATTERNS, ensure_output_dirs

ensure_output_dirs()

patterns = pd.read_csv(THREE_PATTERNS)

PATTERN_LABELS = {
    "categorical_blindspot_HIST": "Histamine blind spot\n(categorical gap)",
    "ORN_SER_confusion": "ORN serotonin confusion\n(classifier error)",
    "Dm_GLUT_confusion": "Dm glutamate confusion\n(classifier error)",
}
PATTERN_COLORS = {
    "categorical_blindspot_HIST": "#d62728",
    "ORN_SER_confusion": "#1f77b4",
    "Dm_GLUT_confusion": "#2ca02c",
}

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

for ax, (pattern_key, label) in zip(axes, PATTERN_LABELS.items()):
    subset = patterns[patterns["pattern"] == pattern_key].sort_values(
        "fafb_entropy", ascending=True
    )
    y = np.arange(len(subset))
    colors = PATTERN_COLORS[pattern_key]

    ax.barh(y, subset["fafb_entropy"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(subset["fafb_cell_type"], fontsize=8)
    ax.set_xlabel("FAFB entropy (bits)")
    ax.set_title(label, fontsize=10)
    ax.axvline(0.3, color="gray", linestyle="--", alpha=0.5, linewidth=1)

    for i, (_, row) in enumerate(subset.iterrows()):
        ax.text(
            row["fafb_entropy"] + 0.02,
            i,
            f"{row['fafb_dominant_nt']}→{row['mcns_confirmed_nt']}",
            va="center",
            fontsize=7,
            color="dimgray",
        )

axes[0].set_ylabel("Cell type")
fig.suptitle(
    "Three systematic NT confusion patterns in FAFB\n"
    "(MCNS-confirmed consistent types with elevated FAFB entropy)",
    fontsize=12,
    y=1.02,
)
plt.tight_layout()

out = FIGURES / "three_confusion_patterns.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved {out}")

# Donut chart: pattern counts
fig, ax = plt.subplots(figsize=(6, 6))
counts = patterns.groupby("pattern").size()
labels = [PATTERN_LABELS.get(k, k).replace("\n", " ") for k in counts.index]
colors = [PATTERN_COLORS.get(k, "gray") for k in counts.index]
ax.pie(
    counts.values,
    labels=[f"{lab}\n(n={n})" for lab, n in zip(labels, counts.values)],
    colors=colors,
    autopct="%1.0f%%",
    startangle=90,
    wedgeprops={"width": 0.45, "edgecolor": "white"},
)
ax.set_title("Flagged cell types by pattern category")
plt.tight_layout()
out2 = FIGURES / "pattern_counts_donut.png"
plt.savefig(out2, dpi=150)
plt.close()
print(f"Saved {out2}")
