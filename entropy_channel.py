"""
Reconstructs the entropy channel's permutation-based significance (p-value,
q-value) from the aggregated entropy table alone -- no raw per-neuron table
needed.

WHY THIS EXISTS
----------------
analysis.py's `stratified_permutation_null` computes an exact one-sided
permutation p-value for every cell type's entropy from the *raw* per-neuron
annotation table, and does save that p-value (and its BH q-value) into its
output. But the currently-committed `results/entropy_corrected.csv` snapshot
in this bundle predates that column being written (it only has `z_score`),
and `data/merged_annotations.csv` (the raw table) is gitignored and not part
of this bundle -- so the p-value/q-value can't simply be read back or re-run
from source.

This reconstructs it anyway, *exactly*, from what the committed table already
has: each type's `nt_distribution` counts and `n_neurons`. The permutation
argument still applies without the raw table because of a standard fact
about random partitions: under `stratified_permutation_null`'s full-dataset
label shuffle, a single type's resulting category counts are marginally
distributed as an exact multivariate hypergeometric draw -- population =
the pooled counts across every labeled neuron in the dataset, sample size =
that type's own n_neurons. Sampling directly from that exact distribution
reproduces the same per-type null the original function would have built
from the raw table. This is not a normal-theory approximation on z (which is
what an earlier version of the dual-channel check used, out of necessity --
see CHANGELOG.md); it's the same permutation logic, applied to a sufficient
statistic (pooled counts + group size) instead of the full raw table. See
`validate_reconstruction_against_real_zscores` for a direct empirical check
of this claim against the project's own real, already-computed z-scores.
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
    """
    df needs columns: cell_type, n_neurons, entropy, counts (dict, NT->count).
    Returns cell_type, entropy_p_value, entropy_q_value (BH within this
    table), entropy_z_reconstructed -- statistically equivalent to what
    analysis.py's stratified_permutation_null would report, reconstructed
    without the raw per-neuron table (see module docstring).
    """
    from analysis import benjamini_hochberg  # reuse the project's one BH implementation

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
    """
    Sanity check: reconstruct significance while *ignoring* the real z_score
    column, then compare reconstructed vs. real z side by side. If the
    reconstruction argument above is correct, these should agree closely
    (up to Monte Carlo noise) with no systematic bias -- this is a direct,
    checkable claim, not an assertion.
    """
    recon = reconstruct_entropy_significance(df_with_real_z, **kwargs)
    merged = df_with_real_z[["cell_type", "z_score"]].merge(recon, on="cell_type")
    merged["abs_diff"] = (merged["z_score"] - merged["entropy_z_reconstructed"]).abs()
    return merged
