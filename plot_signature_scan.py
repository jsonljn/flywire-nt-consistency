"""Visualize confusion-signature scan: simplex neighborhoods vs entropy."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paths import FIGURES, RESULTS, ensure_output_dirs

ensure_output_dirs()

scored = pd.read_csv(RESULTS / "signature_scan.csv")

PATTERN_COLORS = {
    "histamine_blindspot": "#d62728",
    "ORN_SER_confusion": "#1f77b4",
    "Dm_GLUT_confusion": "#2ca02c",
}
PATTERN_LABELS = {
    "histamine_blindspot": "Histamine blind spot",
    "ORN_SER_confusion": "ORN serotonin confusion",
    "Dm_GLUT_confusion": "Dm glutamate confusion",
}


# ── Simplex projection: SER vs GLUT, size = n, alpha by neighborhood ──
fig, ax = plt.subplots(figsize=(10, 7))
background = scored[~scored["in_neighborhood"]]
ax.scatter(
    background["p_GLUT"],
    background["p_SER"],
    s=np.clip(background["n_neurons"] / 8, 8, 80),
    c="#d0d0d0",
    alpha=0.45,
    linewidths=0,
    label="Other FAFB types",
    zorder=1,
)

for pattern, color in PATTERN_COLORS.items():
    subset = scored[(scored["best_pattern"] == pattern) & scored["in_neighborhood"]]
    seeds = subset[subset["is_seed"]]
    novel = subset[subset["is_novel_candidate"]]
    known = subset[~subset["is_seed"] & ~subset["is_novel_candidate"]]
    ax.scatter(
        seeds["p_GLUT"],
        seeds["p_SER"],
        s=np.clip(seeds["n_neurons"] / 6, 28, 140),
        c=color,
        alpha=0.95,
        edgecolors="black",
        linewidths=0.6,
        label=f"{PATTERN_LABELS[pattern]} (seed)",
        zorder=3,
    )
    if len(known):
        ax.scatter(
            known["p_GLUT"],
            known["p_SER"],
            s=np.clip(known["n_neurons"] / 6, 20, 100),
            c=color,
            alpha=0.7,
            marker="s",
            linewidths=0,
            label=f"{PATTERN_LABELS[pattern]} (already flagged)",
            zorder=2,
        )
    if len(novel):
        ax.scatter(
            novel["p_GLUT"],
            novel["p_SER"],
            s=np.clip(novel["n_neurons"] / 6, 24, 110),
            c=color,
            alpha=0.9,
            marker="D",
            edgecolors="black",
            linewidths=0.4,
            label=f"{PATTERN_LABELS[pattern]} (novel)",
            zorder=4,
        )

for _, row in scored[scored["is_seed"] | scored["is_novel_candidate"]].iterrows():
    if row["n_neurons"] < 30 and not row["is_seed"]:
        continue
    ax.annotate(
        row["cell_type"],
        (row["p_GLUT"], row["p_SER"]),
        fontsize=7,
        xytext=(4, 4),
        textcoords="offset points",
    )

ax.set_xlabel("FAFB P(GLUT)")
ax.set_ylabel("FAFB P(SER)")
ax.set_title(
    "Confusion fingerprints on the FAFB NT simplex\n"
    "Seeds define neighborhoods; diamonds are types the name-matched scan never saw"
)
ax.legend(fontsize=7, loc="upper right")
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
plt.tight_layout()
out = FIGURES / "signature_simplex.png"
plt.savefig(out, dpi=150)
plt.close()
print(f"Saved {out}")

# ── Why entropy fails for R1-6 ──
fig, ax = plt.subplots(figsize=(9, 6))
ax.scatter(
    scored["entropy"],
    scored["js_histamine_blindspot"],
    s=18,
    c="#bbbbbb",
    alpha=0.5,
    linewidths=0,
    label="All types",
)
hist = scored[scored["cell_type"].isin(["R7", "R8", "R1-6", "Lai"])]
for _, row in hist.iterrows():
    ax.scatter(
        row["entropy"],
        row["js_histamine_blindspot"],
        s=70,
        c=PATTERN_COLORS["histamine_blindspot"],
        zorder=3,
        edgecolors="black",
        linewidths=0.5,
    )
    ax.annotate(
        row["cell_type"],
        (row["entropy"], row["js_histamine_blindspot"]),
        fontsize=9,
        xytext=(6, 4),
        textcoords="offset points",
    )
ax.set_xlabel("Shannon entropy (bits)")
ax.set_ylabel("JS divergence to histamine prototype (bits)")
ax.set_title(
    "Entropy misses R1-6; simplex distance does not\n"
    "R1-6 is large and 82% ACH, so Dale's-law z-score is strongly negative"
)
ax.legend(fontsize=8)
plt.tight_layout()
out2 = FIGURES / "signature_vs_entropy.png"
plt.savefig(out2, dpi=150)
plt.close()
print(f"Saved {out2}")

# ── Novel candidates ──
novel = scored[scored["is_novel_candidate"]].sort_values("best_js").head(20)
if len(novel):
    fig, ax = plt.subplots(figsize=(10, max(4, 0.38 * len(novel) + 1.5)))
    colors = [PATTERN_COLORS[p] for p in novel["best_pattern"]]
    y = np.arange(len(novel))
    ax.barh(y, novel["best_js"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(novel["cell_type"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Leave-one-out JS divergence to nearest confusion prototype (bits)")
    ax.set_title("Novel simplex-neighborhood candidates (not in the MCNS name-matched list)")
    plt.tight_layout()
    out3 = FIGURES / "signature_novel_candidates.png"
    plt.savefig(out3, dpi=150)
    plt.close()
    print(f"Saved {out3}")
