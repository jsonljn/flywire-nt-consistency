"""Calibrate confusion-signature distances against non-seed cell types.

For a candidate and k available seeds, compare the nearest-seed
Jensen-Shannon divergence with distances to M reference types. If i_obs
reference types are at least as close, use:

    q_close = (i_obs + 0.5) / (M + 1)
    p = 1 - (1 - q_close) ** k

Seed candidates use leave-one-out comparisons. Reference profiles represent
cell-type diversity; sampling around pooled transmitter frequencies would
concentrate on the dataset average. Small reference pools limit significance
resolution, including after continuity correction.

Report within-pattern Benjamini-Hochberg q-values. Candidate selection uses
p < 0.05 and requires subsequent literature validation. Rank patterns by
p-value, then q-value and observed divergence. Entropy significance provides
complementary evidence for types with atypical seed geometry.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from nt_simplex import batch_js_divergence

ALPHA = 0.05


def build_reference_pool(df: pd.DataFrame, all_seed_names: set[str]) -> pd.DataFrame:
    """Select reference cell types that are not seeds of any pattern."""
    return df[~df["cell_type"].astype(str).isin(all_seed_names)].reset_index(drop=True)


def exact_p_value(observed: float, pool_distances: np.ndarray, k: int) -> tuple[float, int, int]:
    """Compute the probability of an equally close reference match in k independent draws.

    Use q_close = (i_obs + 0.5) / (M + 1) and p = 1 - (1 - q_close) ** k.
    Continuity correction avoids zero probabilities when no reference type
    matches the observed distance. Small M limits significance resolution.
    Returns (p_value, i_obs, M).
    """
    M = len(pool_distances)
    if M == 0:
        return 1.0, 0, 0
    i_obs = int(np.sum(pool_distances <= observed))
    q_close = (i_obs + 0.5) / (M + 1)
    p_value = 1.0 - (1.0 - q_close) ** k
    return float(p_value), i_obs, M


def calibrate(
    df: pd.DataFrame,
    pattern_seeds: dict[str, tuple[str, ...]],
    lookup: dict[str, np.ndarray],
) -> pd.DataFrame:
    """For every (cell_type, pattern) pair, compute the observed nearest-seed JS
    distance (leave-one-out when the type is itself a seed of that pattern)
    and the exact empirical-reference-pool p-value described above. Returns a
    long dataframe with one row per (cell_type, pattern); q_value is
    BH-corrected within each pattern's family of tests.
    """
    from analysis import benjamini_hochberg

    all_seed_names = {s for seeds in pattern_seeds.values() for s in seeds}
    pool = build_reference_pool(df, all_seed_names)
    pool_names = pool["cell_type"].astype(str).to_numpy()
    pool_matrix = np.stack(pool["p_vec"].to_numpy())

    records = []
    for _, row in df.iterrows():
        name = str(row["cell_type"])
        p_vec = row["p_vec"]

        # Reuse candidate-to-reference distances across patterns.
        if name in pool_names:
            keep = pool_names != name
            pool_dists = batch_js_divergence(pool_matrix[keep], p_vec)
        else:
            pool_dists = batch_js_divergence(pool_matrix, p_vec)

        for pattern, seeds in pattern_seeds.items():
            is_own_seed = name in seeds
            active_names = [s for s in seeds if s in lookup and s != (name if is_own_seed else None)]
            active_seeds = [lookup[s] for s in active_names]
            if not active_seeds:
                continue

            dists_to_active_seeds = np.array([
                float(batch_js_divergence(p_vec[None, :], s)[0]) for s in active_seeds
            ])
            observed = float(dists_to_active_seeds.min())
            k = len(active_seeds)
            p_value, i_obs, M = exact_p_value(observed, pool_dists, k)

            records.append({
                "cell_type": name,
                "pattern": pattern,
                "n_neurons": int(row["n_neurons"]),
                "observed_js": observed,
                "is_pattern_seed": is_own_seed,
                "n_active_seeds": k,
                "reference_pool_size": M,
                "reference_pool_at_least_as_close": i_obs,
                "p_value": p_value,
            })

    long_df = pd.DataFrame(records)
    long_df["q_value"] = long_df.groupby("pattern")["p_value"].transform(
        lambda s: benjamini_hochberg(s.to_numpy())
    )
    return long_df


def summarize_best_pattern(long_df: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Select one pattern per cell type by calibrated p-value.

    Use q-value and observed divergence as tie breakers. Report raw p-values
    and within-pattern BH q-values; significant denotes p < alpha.
    Raw divergence alone is not comparable across patterns with different seed
    spreads. Geometric matches require independent transmitter validation.
    """
    wide_q = long_df.pivot(index="cell_type", columns="pattern", values="q_value")
    wide_q.columns = [f"q_{c}" for c in wide_q.columns]
    wide_p = long_df.pivot(index="cell_type", columns="pattern", values="p_value")
    wide_p.columns = [f"p_{c}" for c in wide_p.columns]

    best = (
        long_df.sort_values(["cell_type", "p_value", "q_value", "observed_js"])
        .groupby("cell_type", as_index=False)
        .first()[["cell_type", "pattern", "observed_js", "p_value", "q_value"]]
        .rename(columns={
            "pattern": "best_pattern_calibrated",
            "observed_js": "best_js_calibrated",
            "p_value": "best_p_calibrated",
            "q_value": "best_q_calibrated",
        })
    )
    best["significant"] = best["best_p_calibrated"] < alpha

    out = best.merge(wide_q.reset_index(), on="cell_type").merge(wide_p.reset_index(), on="cell_type")
    return out
