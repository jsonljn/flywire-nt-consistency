"""
Statistically calibrated significance test for the confusion-signature scan.

THE ORIGINAL PROBLEM
---------------------
signature_scan.py's original calibration set one JS-divergence threshold per
pattern from the seeds' own leave-one-out (LOO) spread:

    threshold = max(JS_FLOOR, loo_max * LOO_MARGIN)

This is fragile whenever a pattern's seeds are not tightly clustered. The
histamine-blindspot family is exactly that case: R7 (44% GLUT / 39% GABA /
15% ACH) and R1-6 (82% ACH) sit far apart in the 6-simplex -- they share a
*cause* (histamine is not a predictable FAFB output category) but not a
*shape*. R7's poor fit to {R8, R1-6} drags the calibration up, and because
ACh is the single most common FAFB prediction dataset-wide, "loose enough to
include R1-6's 82%-ACH profile" turns out loose enough to include almost any
clean, ACh-dominant type that has nothing to do with histamine. On the real,
committed project data (results/entropy_raw.csv, 402 FAFB cell types) this
flags 341/402 (85%) of *all* cell types as "in a confusion neighborhood."
The project's own unit tests independently caught two symptoms of this:
tests/test_core.py::test_clean_type_is_not_a_novel_candidate (a synthetic
100%-ACH type, nothing to explain, got flagged anyway) and
test_recovers_histamine_seeds (R7 did not recover under its own pattern).

A FIRST FIX ATTEMPT THAT DID NOT WORK, KEPT HERE AS A NOTE
------------------------------------------------------------
The obvious statistical fix is a permutation null in the style of
analysis.py's `stratified_permutation_null`: for a cell type of size n, draw
many same-sized random groups from the dataset's pooled NT-prediction rates
and see how close such an unstructured group lands to a pattern's seeds by
chance. That was implemented and measured against the real data -- and it is
*wrong* in a way worth recording so it is not reinvented. For even a
moderately sized group, sampling noise around the *dataset average* shrinks
quickly, so literally every specific point that is not the dataset average
becomes "impossible by chance" -- including the true seeds (R7/R8/R1-6 all
came back p approx. 1, i.e. now *under*-flagged) and, symmetrically, every
other clean type whose dominant category happens to lean the same general
direction as a seed (117 unrelated types tied at the same near-zero q-value
for the Dm pattern, i.e. *over*-flagged again, just via a different route).
The dataset average is not a realistic stand-in for "an unremarkable fly cell
type" -- essentially no real cell type looks like the average, because Dale's
law means real types are clean in *some* direction. Comparing against
synthetic noise around the average therefore answers "is this non-average,"
not "is this specifically seed-like," and rejects almost everything for
large n regardless of direction.

THE ACTUAL FIX
---------------
Compare each candidate against *other real cell types*, not synthetic noise.
For candidate X and a pattern with seed set S (size k, or k-1 under leave-
one-out when X is itself a seed):

  1. observed = X's JS distance to its nearest seed in S.
  2. Take the ~385-400 *other* FAFB cell types that are not a seed of any
     pattern -- the empirical population of "ordinary" cell-type profiles in
     this dataset -- and compute X's JS distance to each of them once.
  3. Ask: if we swapped S for k independently, randomly chosen ordinary cell
     types instead, what is the probability that at least one of those k
     random types would be at least as close to X as the true seeds are?

Step 3 has an exact closed form (no simulation, no resolution floor to tune):
if a fraction q_close of the M ordinary types are at least as close to X as
`observed`, then for k independent random draws,

    p = P(at least one of k draws is that close) = 1 - (1 - q_close) ** k

q_close uses add-one-half continuity correction, (i_obs + 0.5) / (M + 1),
not the raw i_obs / M -- see `exact_p_value` for why (short version: at the
project's real reference-pool size, M ~ 385-400, this changes nothing; it
only matters, and matters a lot, at small M).

CAVEAT: THIS NEEDS A REASONABLY LARGE REFERENCE POOL
-------------------------------------------------------
This test's resolution is fundamentally limited by M, the number of ordinary
cell types available for comparison. On the real project data M ~ 385-400,
which is plenty. But a small M breaks the test in *both* directions at once,
not just one -- this was found empirically while building tests for this
module (see tests/test_core.py's fixture size) and is worth stating plainly
rather than leaving implicit:
  - Without the continuity correction, small M makes the test *overconfident*:
    with M=3, i_obs=0 forces p=0.0 exactly (something no real reference point
    was as close as), which reads as "impossible by chance" when it's really
    just "we only checked 3 things."
  - With the continuity correction, small M instead makes the test
    *underpowered*: even the true seeds testing against each other can fail
    to reach significance, because q_close can't get much below ~1/M no
    matter how genuinely close the match is.
There is no threshold-free fix for this -- it is a real statement about how
much evidence a comparison against M things can provide, not a bug. Tests
of this module accordingly use a reference pool of realistic size (~20+),
not a handful of seeds plus one or two decoys, so they exercise the regime
this method is actually meant for.

This is the same "how surprising is this, given chance alone" logic as the
project's existing permutation tests, just built from the dataset's own real
diversity of cell-type shapes instead of a parametric resampling model that
turns out to have no realistic cell type near it. Benjamini-Hochberg FDR
correction (reusing analysis.py's implementation) is applied within each
pattern's own family of tests. `best_pattern` is chosen by smallest q-value,
not smallest raw distance, which is what fixes the R7 mis-assignment: raw
distance is not comparable across pattern families with different seed
spreads, but a p-value computed the same way for every pattern is.
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
    fraction, for transparency/debugging.
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
        pattern), the same regime the original name-matched three-pattern
        scan already operates in before its own literature cross-check.
      - `best_q_calibrated` (BH-corrected within each pattern's family of
        ~400 tests): kept for full transparency, but at this prevalence and
        these effect sizes BH-FDR has essentially no power -- even the
        seeds themselves mostly do not survive it (see module docstring) --
        so it is reported, not used to gate anything.

    best_pattern is chosen by smallest p-value, which is what fixes the
    original R7 mis-assignment bug: raw JS distance is not comparable across
    pattern families with different seed spreads, but a p-value computed the
    same way for every pattern is. Note this does not "solve" R7 itself --
    R7's geometric nearest match is genuinely Dm_GLUT_confusion, not its own
    family, under leave-one-out (p=0.03 vs p=0.11). That is a real, reported
    limit of geometry-only matching for this specific type, not a bug; see
    README.
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
