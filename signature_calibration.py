"""
Statistically calibrated significance test for confusion-signature matching.

For candidate X and a pattern with k available seeds (excluding X under
leave-one-out), measure X's Jensen-Shannon divergence to its nearest seed.
Compare that distance with X's distances to other non-seed cell types.
If i_obs of M reference types are at least as close, use the continuity-
corrected fraction q_close = (i_obs + 0.5) / (M + 1). For k independent
random reference draws, the probability of at least one equally close match is:

    p = 1 - (1 - q_close) ** k

The reference pool represents the diversity of real cell-type profiles.
Sampling around pooled transmitter frequencies instead concentrates near
the dataset average, which is not representative of the many real types
with a strongly dominant transmitter. Likewise, a shared distance threshold
can include unrelated acetylcholine-dominant types when seeds are dispersed.
R7 and R1-6 illustrate this: both have verified histamine signaling, but their
FAFB prediction profiles differ substantially.

Reference-pool size limits significance resolution. Without continuity
correction, no close matches in a small pool would force p=0; with correction,
a small pool can lack power even for close seed matches. Tests use reference
pools appropriate to their seed counts. These limits apply to interpretation
of the reported evidence as well as to synthetic test fixtures.

Benjamini-Hochberg q-values are reported within each pattern's family of
tests. Candidate selection uses raw p < 0.05, followed by independent
literature validation; q-values provide a stricter evidence measure.
The best pattern is selected by smallest p-value, with q-value and observed
distance as tie breakers. Raw distances alone are not comparable across
patterns with different seed spreads. Entropy significance provides a
complementary channel for types such as R7 whose geometry does not reliably
identify their literature-supported pattern.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from nt_simplex import batch_js_divergence

ALPHA = 0.05


def build_reference_pool(df: pd.DataFrame, all_seed_names: set[str]) -> pd.DataFrame:
    """Cell types that are not a seed of *any* pattern -- the empirical background."""
    return df[~df["cell_type"].astype(str).isin(all_seed_names)].reset_index(drop=True)


def exact_p_value(observed: float, pool_distances: np.ndarray, k: int) -> tuple[float, int, int]:
    """
    P(at least one of k independently, randomly chosen ordinary cell types is
    at least as close to the candidate as `observed`), computed exactly from
    the empirical fraction of the reference pool that is that close.

    Uses add-one-half (Laplace/Krichevsky-Trofimov) continuity correction on
    the underlying proportion, q_close = (i_obs + 0.5) / (M + 1), rather than
    the raw i_obs/M. At the reference-pool sizes this project actually scores
    against (M ~ 385-400), the correction changes nothing meaningful (e.g.
    i_obs=0, M=388: raw 0.0 vs corrected 0.0013). It matters at small M --
    without it, i_obs=0 forces p_value=0.0 exactly regardless of how small M
    is, which is overconfident when M is small (see
    tests/test_core.py::test_small_reference_pool_is_not_overconfident: 3
    reference points is nowhere near enough to call anything "impossible by
    chance," and the uncorrected version does exactly that).

    Returns (p_value, i_obs, M) where i_obs/M is the raw (uncorrected)
    fraction, for inspection of the reference-pool calibration.
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
    """
    For every (cell_type, pattern) pair, compute the observed nearest-seed JS
    distance (leave-one-out when the type is itself a seed of that pattern)
    and the exact empirical-reference-pool p-value described above. Returns a
    long dataframe with one row per (cell_type, pattern); q_value is
    BH-corrected within each pattern's family of tests.
    """
    from analysis import benjamini_hochberg  # reuse the project's one BH implementation

    all_seed_names = {s for seeds in pattern_seeds.values() for s in seeds}
    pool = build_reference_pool(df, all_seed_names)
    pool_names = pool["cell_type"].astype(str).to_numpy()
    pool_matrix = np.stack(pool["p_vec"].to_numpy())

    records = []
    for _, row in df.iterrows():
        name = str(row["cell_type"])
        p_vec = row["p_vec"]

        # This candidate's distance to every ordinary (non-seed) type, computed
        # once and reused across all three patterns.
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
    """
    Collapse the long (cell_type, pattern) table to one row per cell type.

    Two tiers are reported, deliberately mirroring how analysis.py's entropy
    screen already reports both a raw effect-size bar (z > 2) *and* a
    stricter FDR bar (q < 0.05) side by side rather than picking one:

      - `best_p_calibrated` / `significant` (p < alpha): the raw exact
        p-value described in this module's docstring. This is the primary,
        reported tier -- appropriate for a candidate-generating screen with
        very low prevalence (a handful of true members among ~400 tests per
        pattern), the same regime the name-matched three-pattern
        scan already operates in before its own literature cross-check.
      - `best_q_calibrated` (BH-corrected within each pattern's family of
        ~400 tests): kept for full transparency, but at this prevalence and
        these effect sizes BH-FDR has essentially no power -- even the
        seeds themselves mostly do not survive it (see module docstring) --
        so it is reported, not used to gate anything.

    best_pattern is chosen by smallest p-value because raw JS distances
    are not comparable across pattern families with different seed spreads.
    R7 can still match Dm_GLUT_confusion under leave-one-out; entropy provides
    complementary evidence for its histamine classification (see README).
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
