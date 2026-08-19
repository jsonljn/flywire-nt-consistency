"""
Neurotransmitter Consistency Within Cell Types Across FlyWire Connectomes
=========================================================================

Goal: Find which cell types violate Dale's law — where neurons of the same
type disagree on their predicted neurotransmitter — and determine whether
that inconsistency is genuine biology or annotation noise.
"""

from __future__ import annotations

import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

from nt_simplex import entropy_from_count_matrix
from paths import FIGURES, ensure_output_dirs, entropy_paths

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────
# STEP 1: Load data
# ─────────────────────────────────────────────

def load_annotations(path):
    """Load the FlyWire annotations CSV."""
    print(f"Loading {path}...")
    df = pd.read_csv(path)
    print(f"  Shape: {df.shape}")
    print(f"  Columns: {list(df.columns)}")
    return df


def detect_columns(df):
    """Auto-detect relevant columns regardless of exact naming."""
    cols = {c.lower(): c for c in df.columns}
    mapping = {}

    for candidate in ["nt_type", "neurotransmitter", "nt", "predicted_nt", "top_nt"]:
        if candidate in cols:
            mapping["nt"] = cols[candidate]
            break

    for candidate in ["primary_type", "cell_type", "type", "cell type", "celltype", "hemibrain_type"]:
        if candidate in cols:
            mapping["cell_type"] = cols[candidate]
            break

    for candidate in ["hemilineage", "lineage", "birth_time"]:
        if candidate in cols:
            mapping["hemilineage"] = cols[candidate]
            break

    for candidate in ["neuropil", "top_neuropil", "region", "sub_class", "class", "side"]:
        if candidate in cols:
            mapping["neuropil"] = cols[candidate]
            break

    print(f"  Detected columns: {mapping}")
    return mapping


# ─────────────────────────────────────────────
# STEP 2: Shannon entropy per cell type
# ─────────────────────────────────────────────

def shannon_entropy(counts) -> float:
    """Shannon entropy in bits (log base 2)."""
    counts = np.asarray(counts)
    counts = counts[counts > 0]
    probs = counts / counts.sum()
    return float(-np.sum(probs * np.log2(probs)))


def _counts_to_json(counts: dict) -> str:
    """Serialize NT count dict without numpy types."""
    return json.dumps({str(k): int(v) for k, v in counts.items()})


def compute_entropy_per_type(df, cell_type_col, nt_col, min_members=20):
    """Compute Shannon entropy of NT distribution per cell type."""
    results = []

    for cell_type, group in df.groupby(cell_type_col):
        if len(group) < min_members:
            continue

        nt_counts = group[nt_col].value_counts()
        entropy = shannon_entropy(nt_counts.values)

        results.append({
            "cell_type": cell_type,
            "n_neurons": len(group),
            "n_nt_types": int((nt_counts > 0).sum()),
            "entropy": entropy,
            "dominant_nt": nt_counts.index[0],
            "dominant_frac": nt_counts.iloc[0] / len(group),
            "nt_distribution": _counts_to_json(nt_counts.to_dict()),
        })

    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values("entropy", ascending=False).reset_index(drop=True)
    print(f"\nCell types with >={min_members} members: {len(df_results)}")
    print(f"Cell types with entropy > 0 (any inconsistency): {(df_results['entropy'] > 0).sum()}")
    return df_results


# ─────────────────────────────────────────────
# STEP 3: Stratified permutation null
# ─────────────────────────────────────────────

def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg FDR q-values for a vector of p-values."""
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p

    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        prev = min(prev, value)
        q[i] = prev
    q = np.clip(q, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def stratified_permutation_null(
    df,
    cell_type_col,
    nt_col,
    min_members=20,
    n_permutations=1000,
    seed=42,
):
    """
    Stratified permutation null preserving marginal NT counts and group sizes.

    Returns raw entropy, null moments, z-score, and empirical one-sided p-value
    P(null >= observed), with Benjamini-Hochberg q-values across all types.
    """
    rng = np.random.default_rng(seed)

    eligible = df.groupby(cell_type_col).filter(lambda x: len(x) >= min_members)
    type_codes, type_names = pd.factorize(eligible[cell_type_col], sort=False)
    nt_codes, _nt_names = pd.factorize(eligible[nt_col].astype(str), sort=False)
    n_types = len(type_names)
    n_nt = len(_nt_names)
    sizes = np.bincount(type_codes, minlength=n_types)

    observed_counts = np.zeros((n_types, n_nt), dtype=np.int64)
    np.add.at(observed_counts, (type_codes, nt_codes), 1)
    observed = entropy_from_count_matrix(observed_counts)

    null_sum = np.zeros(n_types, dtype=np.float64)
    null_sumsq = np.zeros(n_types, dtype=np.float64)
    null_ge = np.zeros(n_types, dtype=np.int64)

    for _ in range(n_permutations):
        shuffled = rng.permutation(nt_codes)
        counts = np.zeros((n_types, n_nt), dtype=np.int64)
        np.add.at(counts, (type_codes, shuffled), 1)
        ent = entropy_from_count_matrix(counts)
        null_sum += ent
        null_sumsq += ent * ent
        null_ge += (ent >= observed - 1e-12).astype(np.int64)

    mean_null = null_sum / n_permutations
    var_null = np.maximum(null_sumsq / n_permutations - mean_null ** 2, 0.0)
    std_null = np.sqrt(var_null)
    z = np.zeros(n_types, dtype=np.float64)
    positive_std = std_null > 0
    z[positive_std] = (observed[positive_std] - mean_null[positive_std]) / std_null[positive_std]
    # (k+1)/(n+1) avoids zero p-values under permutation
    p_values = (null_ge + 1) / (n_permutations + 1)

    results = []
    for i, ct in enumerate(type_names):
        results.append({
            "cell_type": ct,
            "n_neurons": int(sizes[i]),
            "entropy": float(observed[i]),
            "mean_null_entropy": float(mean_null[i]),
            "std_null_entropy": float(std_null[i]),
            "z_score": float(z[i]),
            "p_value": float(p_values[i]),
        })

    df_null = pd.DataFrame(results)
    df_null["q_value"] = benjamini_hochberg(df_null["p_value"].values)
    df_null = df_null.sort_values("z_score", ascending=False).reset_index(drop=True)
    return df_null


# ─────────────────────────────────────────────
# STEP 4: Biological cross-reference
# ─────────────────────────────────────────────

def crossref_outliers(
    df_annotations,
    df_entropy,
    cell_type_col,
    nt_col,
    hemilineage_col=None,
    neuropil_col=None,
    top_n=10,
):
    """Print hemilineage/neuropil breakdown for top outlier cell types."""
    top_types = df_entropy.head(top_n)["cell_type"].tolist()

    print(f"\nTop {top_n} outlier cell types:")
    print("=" * 60)

    for ct in top_types:
        group = df_annotations[df_annotations[cell_type_col] == ct]
        nt_dist = group[nt_col].value_counts()
        entropy_row = df_entropy[df_entropy["cell_type"] == ct].iloc[0]

        print(f"\nCell type: {ct}")
        print(f"  N neurons: {len(group)}")
        print(f"  Entropy: {entropy_row['entropy']:.3f}  |  Z-score: {entropy_row.get('z_score', 'N/A')}")
        print(f"  NT distribution: {dict(nt_dist)}")

        if hemilineage_col and hemilineage_col in group.columns:
            print(f"  Top hemilineages: {dict(group[hemilineage_col].value_counts().head(3))}")

        if neuropil_col and neuropil_col in group.columns:
            print(f"  Top neuropils: {dict(group[neuropil_col].value_counts().head(3))}")


# ─────────────────────────────────────────────
# STEP 5: Visualizations
# ─────────────────────────────────────────────

def plot_entropy_distribution(df_entropy, output_path):
    top30 = df_entropy.head(30)
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(range(len(top30)), top30["entropy"], color="steelblue", alpha=0.8)
    ax.set_xticks(range(len(top30)))
    ax.set_xticklabels(top30["cell_type"], rotation=90, fontsize=7)
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.set_title("Top 30 Cell Types by Neurotransmitter Inconsistency (Raw Entropy)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved: {output_path}")


def plot_zscore_distribution(df_null, output_path, min_members=20):
    top30 = df_null.head(30)
    sig = df_null[df_null["q_value"] < 0.05]

    fig, ax = plt.subplots(figsize=(14, 6))
    colors = ["firebrick" if ct in set(sig["cell_type"]) else "indianred" for ct in top30["cell_type"]]
    ax.bar(range(len(top30)), top30["z_score"], color=colors, alpha=0.85)
    ax.set_xticks(range(len(top30)))
    ax.set_xticklabels(top30["cell_type"], rotation=90, fontsize=7)
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Z-score (size-corrected)")
    ax.set_title(
        f"Top 30 Cell Types by Size-Corrected Inconsistency (n>={min_members})\n"
        "dark red = FDR q < 0.05"
    )
    ax.axhline(y=2, color="black", linestyle="--", alpha=0.5, label="z=2 (reference)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved: {output_path}")


def plot_entropy_vs_size(df_null, output_path):
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(
        df_null["n_neurons"],
        df_null["entropy"],
        c=df_null["z_score"],
        cmap="RdYlGn_r",
        alpha=0.6,
        s=30,
    )
    plt.colorbar(sc, ax=ax, label="Z-score")
    ax.set_xlabel("Number of neurons in cell type")
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.set_title("Entropy vs Cell Type Size\n(color = size-corrected z-score)")
    ax.set_xscale("log")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved: {output_path}")


def plot_fdr_volcano(df_null, output_path, min_members=20):
    """Volcano-style plot: z-score vs -log10(FDR q-value)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    x = df_null["z_score"].values
    q = np.clip(df_null["q_value"].values, 1e-300, 1.0)
    y = -np.log10(q)
    sig = df_null["q_value"] < 0.05

    ax.scatter(x[~sig], y[~sig], alpha=0.45, s=28, color="gray", label="not significant")
    ax.scatter(x[sig], y[sig], alpha=0.85, s=40, color="firebrick", label="FDR q < 0.05")

    for _, row in df_null[sig].head(8).iterrows():
        ax.annotate(
            row["cell_type"],
            (row["z_score"], -np.log10(max(row["q_value"], 1e-300))),
            fontsize=7,
            xytext=(4, 4),
            textcoords="offset points",
        )

    ax.axhline(-np.log10(0.05), color="black", linestyle="--", alpha=0.4, label="q = 0.05")
    ax.set_xlabel("Z-score (size-corrected)")
    ax.set_ylabel("-log10(FDR q-value)")
    ax.set_title(f"Significance vs effect size (n>={min_members})")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved: {output_path}")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def run_analysis(annotations_path, min_members=20, n_permutations=1000):
    ensure_output_dirs()
    raw_path, corrected_path = entropy_paths(min_members)

    df = load_annotations(annotations_path)
    col_map = detect_columns(df)

    cell_type_col = col_map.get("cell_type")
    nt_col = col_map.get("nt")
    hemilineage_col = col_map.get("hemilineage")
    neuropil_col = col_map.get("neuropil")

    if not cell_type_col or not nt_col:
        print("\nERROR: Could not detect cell_type or nt column.")
        return None

    df = df.dropna(subset=[cell_type_col, nt_col])
    print(f"\nRows after dropping NA in key columns: {len(df)}")

    print("\n--- Computing raw entropy per cell type ---")
    df_entropy = compute_entropy_per_type(df, cell_type_col, nt_col, min_members)
    df_entropy.to_csv(raw_path, index=False)
    print(f"Saved: {raw_path}")

    print(f"\n--- Running stratified permutation null ({n_permutations} permutations) ---")
    df_null = stratified_permutation_null(
        df, cell_type_col, nt_col, min_members, n_permutations
    )

    df_null = df_null.merge(
        df_entropy[["cell_type", "dominant_nt", "dominant_frac", "n_nt_types", "nt_distribution"]],
        on="cell_type",
        how="left",
    )
    df_null.to_csv(corrected_path, index=False)
    print(f"Saved: {corrected_path}")

    crossref_outliers(
        df, df_null, cell_type_col, nt_col, hemilineage_col, neuropil_col, top_n=10
    )

    suffix = "" if min_members == 20 else f"_n{min_members}"
    print("\n--- Generating plots ---")
    plot_entropy_distribution(df_entropy, FIGURES / f"entropy_distribution{suffix}.png")
    plot_zscore_distribution(df_null, FIGURES / f"zscore_distribution{suffix}.png", min_members)
    plot_entropy_vs_size(df_null, FIGURES / f"entropy_vs_size{suffix}.png")
    plot_fdr_volcano(df_null, FIGURES / f"fdr_volcano{suffix}.png", min_members)

    sig_z = df_null[df_null["z_score"] > 2]
    sig_fdr = df_null[df_null["q_value"] < 0.05]
    print("\n=== SUMMARY ===")
    print(f"Total cell types analyzed: {len(df_null)}")
    print(f"Outliers (z > 2): {len(sig_z)}")
    print(f"Outliers (FDR q < 0.05): {len(sig_fdr)}")
    print("\nTop 10 by z-score:")
    print(
        df_null[
            ["cell_type", "n_neurons", "entropy", "z_score", "p_value", "q_value", "dominant_nt"]
        ].head(10).to_string(index=False)
    )

    return df, df_entropy, df_null


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Neurotransmitter consistency analysis")
    parser.add_argument("path", nargs="?", default="data/merged_annotations.csv")
    parser.add_argument("--min-members", type=int, default=20)
    parser.add_argument("--n-permutations", type=int, default=1000)
    args = parser.parse_args()
    run_analysis(args.path, min_members=args.min_members, n_permutations=args.n_permutations)
