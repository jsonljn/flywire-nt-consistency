import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from paths import FIGURES, RESULTS, ensure_output_dirs

ensure_output_dirs()

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, cell_type in zip(axes, ["R7", "R8"]):
    df = pd.read_csv(RESULTS / f"connectivity_profiles_{cell_type}.csv")
    label_col = df["nt_label"]
    X = df.drop(columns=["root_id", "nt_label"]).values

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)

    colors = {
        "GLUT": "#1f77b4", "GABA": "#ff7f0e", "ACH": "#2ca02c",
        "OCT": "#d62728", "DA": "#9467bd", "SER": "#8c564b",
    }

    for label in label_col.unique():
        mask = label_col == label
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            label=f"{label} (n={mask.sum()})",
            alpha=0.6, s=25, color=colors.get(label, "gray"),
        )

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.set_title(f"{cell_type} output connectivity profile\ncolored by (wrong) NT prediction")
    ax.legend(fontsize=8)

plt.tight_layout()
out_path = FIGURES / "connectivity_pca.png"
plt.savefig(out_path, dpi=150)
plt.close()
print(f"Saved {out_path}")
