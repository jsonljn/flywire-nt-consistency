"""
Correction summary visualization: neurons flagged per cell type and pattern.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from paths import CORRECTIONS, FIGURES, LITERATURE_VALIDATED, ensure_output_dirs

ensure_output_dirs()

fafb = pd.read_csv(CORRECTIONS / "corrections_fafb.csv")
lit = pd.read_csv(LITERATURE_VALIDATED)

pattern_map = lit.set_index("fafb_cell_type")["pattern"].to_dict()
fafb["pattern"] = fafb["cell_type"].map(pattern_map)

# Bar chart: corrections per cell type, colored by pattern
PATTERN_COLORS = {
    "categorical_blindspot_HIST": "#d62728",
    "ORN_SER_confusion": "#1f77b4",
    "Dm_GLUT_confusion": "#2ca02c",
}

counts = fafb.groupby(["cell_type", "pattern"]).size().reset_index(name="n")
counts = counts.sort_values("n", ascending=True)

fig, ax = plt.subplots(figsize=(10, 8))
colors = [PATTERN_COLORS.get(p, "gray") for p in counts["pattern"]]
ax.barh(counts["cell_type"], counts["n"], color=colors, alpha=0.85)
ax.set_xlabel("Neurons proposed for correction")
ax.set_ylabel("Cell type")
ax.set_title(f"FAFB correction list: {len(fafb)} neurons across {fafb['cell_type'].nunique()} cell types")

# Legend
from matplotlib.patches import Patch
legend_items = [
    Patch(color=c, label=p.replace("_", " "))
    for p, c in PATTERN_COLORS.items()
    if p in counts["pattern"].values
]
ax.legend(handles=legend_items, loc="lower right", fontsize=8)
plt.tight_layout()

out = FIGURES / "corrections_summary.png"
plt.savefig(out, dpi=150)
plt.close()
print(f"Saved {out}")

# Stacked bar: wrong NT categories per pattern
fig, ax = plt.subplots(figsize=(8, 5))
cross = pd.crosstab(fafb["pattern"], fafb["current_predicted_nt"])
cross.plot(kind="bar", stacked=True, ax=ax, colormap="tab10", alpha=0.85)
ax.set_xlabel("Pattern")
ax.set_ylabel("Neurons")
ax.set_title("Wrong NT predictions in correction list, by pattern")
ax.legend(title="Predicted NT", bbox_to_anchor=(1.02, 1), fontsize=8)
plt.xticks(rotation=15, ha="right")
plt.tight_layout()
out2 = FIGURES / "corrections_by_wrong_nt.png"
plt.savefig(out2, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved {out2}")
