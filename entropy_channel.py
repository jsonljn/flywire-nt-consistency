"""Estimate entropy significance from aggregated neurotransmitter counts.

Input tables provide each cell type's nt_distribution and n_neurons.
Under a label-shuffle null, a group's category counts follow a multivariate
hypergeometric distribution parameterized by pooled category counts and
group size. Sampling this distribution reconstructs per-type null entropies
without a raw per-neuron annotation table.

The reference population is the pooled counts in the supplied table, so its
coverage determines the null. Monte Carlo estimates provide one-sided
p-values, within-table BH q-values, and reconstructed z-scores.
validate_reconstruction_against_real_zscores compares these estimates with
the saved analysis z-scores.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from nt_simplex import NT_ORDER, counts_to_vector, entropy_from_count_matrix


def pooled_counts(df: pd.DataFrame, counts_col: str = "counts", order=NT_ORDER) -> np.ndarray:
    """Sum every row's NT count dict into one dataset-wide count vector."""
    total = np.zeros(len(order), dtype=np.int64)
    for counts in df[counts_col]:
        total += counts_to_vector(counts, order).astype(np.int64)
    return total


def reconstruct_entropy_significance(
    df: pd.DataFrame,
    n_permutations: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    """Estimate entropy p-values, BH q-values, and z-scores from aggregated counts.

    Requires cell_type, n_neurons, entropy, and counts columns. The reference
    population consists of the pooled category counts in the supplied table.
    """
    from analysis import benjamini_hochberg

    rng = np.random.default_rng(seed)
    pooled = pooled_counts(df)
    total_pool = int(pooled.sum())

    null_cache: dict[int, np.ndarray] = {}
    records = []
    for _, row in df.iterrows():
        n = min(int(row["n_neurons"]), total_pool)
        if n not in null_cache:
            draws = rng.multivariate_hypergeometric(pooled, n, size=n_permutations)
            null_cache[n] = entropy_from_count_matrix(draws.astype(np.float64))
        null_ent = null_cache[n]

        observed = float(row["entropy"])
        p_value = (np.sum(null_ent >= observed - 1e-12) + 1) / (n_permutations + 1)
        mean_null, std_null = float(null_ent.mean()), float(null_ent.std())
        z = (observed - mean_null) / std_null if std_null > 0 else 0.0

        records.append({
            "cell_type": row["cell_type"],
            "entropy_p_value": float(p_value),
            "entropy_z_reconstructed": float(z),
        })

    out = pd.DataFrame(records)
    out["entropy_q_value"] = benjamini_hochberg(out["entropy_p_value"].to_numpy())
    return out


def validate_reconstruction_against_real_zscores(df_with_real_z: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Compare reconstructed z-scores with saved values, allowing for Monte Carlo variation."""
    recon = reconstruct_entropy_significance(df_with_real_z, **kwargs)
    merged = df_with_real_z[["cell_type", "z_score"]].merge(recon, on="cell_type")
    merged["abs_diff"] = (merged["z_score"] - merged["entropy_z_reconstructed"]).abs()
    return merged
