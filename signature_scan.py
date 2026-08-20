"""
Confusion-signature scan.

Entropy is a scalar: it detects *that* a cell type is inconsistent, but not
*why*, and it systematically misses large types that are consistently assigned
the *wrong* transmitter (R1-6: 82% ACH, z = -16, yet canonically
histaminergic).

This module treats each cell type as a point on the 6-simplex of FAFB
classifier outputs and scores it against literature-confirmed confusion
fingerprints. That recovers the original three patterns *without
name-matching to MCNS*, and surfaces additional types that sit in the same
simplex neighborhoods — candidates the name-matched scan could never see.

CALIBRATION -- read this before trusting a number out of this file
--------------------------------------------------------------------
Membership ("in_neighborhood" / "is_novel_candidate" / "best_pattern") is
decided by `signature_calibration.py`'s exact permutation p-value against an
empirical reference-pool null, not by a fixed JS-divergence threshold. See
that module's docstring for why and how: the earlier fixed-threshold
heuristic flagged 80% of all FAFB cell types (322/402) as "novel candidates"
on the real project data, which the project's own unit tests independently
caught as wrong (a synthetic 100%-ACH type has nothing to do with histamine,
and got flagged anyway). The old heuristic's numbers are still computed and
kept in this file's output under an explicit `_heuristic` suffix, purely so
the fix is auditable side-by-side with what shipped before.

DUAL-CHANNEL RECOVERY
----------------------
R1-6 (and, on real data, Dm12/Dm1) are only reachable through simplex
geometry: their entropy z-scores are negative (unusually *consistent*, not
inconsistent), which the entropy screen structurally cannot flag -- that
blind spot is the entire reason this module exists. R7 is close to the
opposite case: its predicted-NT profile is a genuine outlier *within its own
seed family* (nearer, under leave-one-out, to the unrelated Dm_GLUT_confusion
seeds than to R8/R1-6 -- a real geometric fact, not a bug), so simplex
matching alone is not reliable for it. But R7 was never a hard case in the
first place: it is one of the largest entropy z-score outliers in the entire
dataset, already caught by the project's existing, established channel.
`recovery_report` checks both channels and reports which one(s) actually
caught each seed (`recovered_via`) -- nothing recovers silently. The entropy
channel's q-value here is reconstructed exactly (not approximated) from the
committed entropy table via `entropy_channel.py`; see that module for why an
exact reconstruction is possible without the raw per-neuron table, and its
validation against this project's own real z-scores.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import signature_calibration as calib
from entropy_channel import reconstruct_entropy_significance
from nt_simplex import (
    NT_ORDER,
    counts_to_vector,
    js_divergence,
    normalize,
    parse_nt_distribution,
)
from paths import (
    ENTROPY_CORRECTED,
    ENTROPY_CORRECTED_N10,
    ENTROPY_RAW,
    ENTROPY_RAW_N10,
    GT_DATA,
    LITERATURE_VALIDATED,
    RESULTS,
    THREE_PATTERNS,
    ensure_output_dirs,
)

# Literature-confirmed seeds (Lai excluded: literature contradicts MCNS).
# R1-6 is the motivating example: entropy z-score misses it; simplex match does not.
HIST_SEEDS = ("R7", "R8", "R1-6")
ORN_SEEDS = (
    "ORN_V", "ORN_VM3", "ORN_VA2", "ORN_DA3", "ORN_DA4m",
    "ORN_DA4l", "ORN_DM2", "ORN_DM3", "ORN_DL4", "ORN_DL3",
)
DM_SEEDS = ("Dm12", "Dm19", "Dm1")  # literature-confirmed GLUT, FAFB-confused

PATTERN_SEEDS = {
    "histamine_blindspot": HIST_SEEDS,
    "ORN_SER_confusion": ORN_SEEDS,
    "Dm_GLUT_confusion": DM_SEEDS,
}

# Legacy heuristic, kept only for side-by-side comparison -- see module docstring.
JS_FLOOR = 0.12
LOO_MARGIN = 1.35


def load_entropy_table() -> pd.DataFrame:
    """Prefer n>=10 coverage; prefer the z-scored table when present so the
    dual-channel recovery check has a real z_score to reconstruct against."""
    candidates = [ENTROPY_CORRECTED_N10, ENTROPY_CORRECTED, ENTROPY_RAW_N10, ENTROPY_RAW]
    path = next((p for p in candidates if p.exists()), ENTROPY_RAW)
    df = pd.read_csv(path)
    df["counts"] = df["nt_distribution"].map(parse_nt_distribution)
    vectors = df["counts"].map(lambda c: normalize(counts_to_vector(c)))
    df["p_vec"] = list(vectors)
    if "z_score" not in df.columns:
        df["z_score"] = np.nan
    return df


def _vector_lookup(df: pd.DataFrame) -> dict[str, np.ndarray]:
    return {row["cell_type"]: row["p_vec"] for _, row in df.iterrows()}


def nearest_seed_js(
    p: np.ndarray,
    pattern: str,
    lookup: dict[str, np.ndarray],
    exclude: str | None = None,
) -> float:
    """Minimum JS divergence to any seed of `pattern` (leave-one-out)."""
    dists = [
        js_divergence(p, lookup[seed])
        for seed in PATTERN_SEEDS[pattern]
        if seed != exclude and seed in lookup
    ]
    return min(dists) if dists else float("nan")


def _legacy_heuristic_columns(df: pd.DataFrame, lookup: dict[str, np.ndarray]) -> pd.DataFrame:
    """As-shipped fixed-threshold heuristic. Kept for comparison only -- see module docstring."""
    loo_max: dict[str, float] = {}
    for pattern, seeds in PATTERN_SEEDS.items():
        dists = []
        for seed in seeds:
            if seed not in lookup:
                continue
            d = nearest_seed_js(lookup[seed], pattern, lookup, exclude=seed)
            if np.isfinite(d):
                dists.append(d)
        loo_max[pattern] = max(dists) if dists else JS_FLOOR
    thresholds = {p: max(JS_FLOOR, loo_max[p] * LOO_MARGIN) for p in PATTERN_SEEDS}

    rows = []
    for _, row in df.iterrows():
        name = str(row["cell_type"])
        p = row["p_vec"]
        js_by_pattern = {pat: nearest_seed_js(p, pat, lookup, exclude=name) for pat in PATTERN_SEEDS}
        best_pattern = min(js_by_pattern, key=js_by_pattern.get)
        best_js = js_by_pattern[best_pattern]
        rows.append({
            "cell_type": name,
            "js_histamine_blindspot": js_by_pattern["histamine_blindspot"],
            "js_ORN_SER_confusion": js_by_pattern["ORN_SER_confusion"],
            "js_Dm_GLUT_confusion": js_by_pattern["Dm_GLUT_confusion"],
            "best_pattern_heuristic": best_pattern,
            "best_js_heuristic": best_js,
            "pattern_threshold_heuristic": thresholds[best_pattern],
            "in_neighborhood_heuristic": bool(best_js <= thresholds[best_pattern]),
        })
    out = pd.DataFrame(rows)
    out.attrs["thresholds"] = thresholds
    out.attrs["loo_max"] = loo_max
    return out


def score_types(df: pd.DataFrame) -> pd.DataFrame:
    lookup = _vector_lookup(df)
    known_flagged: set[str] = set()
    if THREE_PATTERNS.exists():
        known_flagged = set(pd.read_csv(THREE_PATTERNS)["fafb_cell_type"].astype(str))

    lit_confirmed: set[str] = set()
    if LITERATURE_VALIDATED.exists():
        lit = pd.read_csv(LITERATURE_VALIDATED)
        lit_confirmed = set(
            lit.loc[lit["agrees_with_literature"] == True, "fafb_cell_type"].astype(str)  # noqa: E712
        )

    legacy = _legacy_heuristic_columns(df, lookup)

    long_df = calib.calibrate(df, PATTERN_SEEDS, lookup)
    calibrated = calib.summarize_best_pattern(long_df)

    matched = (
        long_df[long_df["p_value"] < calib.ALPHA]
        .groupby("cell_type")["pattern"]
        .apply(lambda s: ",".join(sorted(s)))
        .rename("matched_patterns")
        .reset_index()
    )

    entropy_sig = reconstruct_entropy_significance(df, n_permutations=20_000)
    entropy_sig["entropy_channel_significant"] = (
        (entropy_sig["entropy_q_value"] < calib.ALPHA) & (entropy_sig["entropy_z_reconstructed"] > 0)
    )

    base_rows = []
    for _, row in df.iterrows():
        name = str(row["cell_type"])
        p = row["p_vec"]
        p_fast = float(sum(p[NT_ORDER.index(nt)] for nt in ("ACH", "GABA", "GLUT")))
        base_rows.append({
            "cell_type": name,
            "n_neurons": int(row["n_neurons"]),
            "entropy": float(row["entropy"]),
            "dominant_nt": row["dominant_nt"],
            "dominant_frac": float(row["dominant_frac"]),
            "p_ACH": float(p[0]),
            "p_GABA": float(p[1]),
            "p_GLUT": float(p[2]),
            "p_DA": float(p[3]),
            "p_SER": float(p[4]),
            "p_OCT": float(p[5]),
            "p_fast_transmitters": p_fast,
            "is_seed": any(name in seeds for seeds in PATTERN_SEEDS.values()),
            "already_name_matched": name in known_flagged,
            "literature_confirmed": name in lit_confirmed,
            "entropy_z": float(row.get("z_score", np.nan)),
        })
    base = pd.DataFrame(base_rows)

    result = (
        base.merge(legacy, on="cell_type", how="left")
        .merge(calibrated, on="cell_type", how="left")
        .merge(matched, on="cell_type", how="left")
        .merge(entropy_sig, on="cell_type", how="left")
    )
    result["matched_patterns"] = result["matched_patterns"].fillna("")

    # Authoritative (calibrated) membership columns, named to match what
    # downstream consumers (plot_signature_scan.py, validate_results.py) expect.
    result["best_pattern"] = result["best_pattern_calibrated"]
    result["best_js"] = result["best_js_calibrated"]
    result["in_neighborhood"] = result["significant"]
    result["is_novel_candidate"] = (
        result["significant"] & ~result["is_seed"] & ~result["already_name_matched"]
    )

    result = result.sort_values(
        ["is_novel_candidate", "significant", "best_q_calibrated"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    result.attrs["thresholds"] = legacy.attrs["thresholds"]
    result.attrs["loo_max"] = legacy.attrs["loo_max"]
    annotated = _annotate_literature(result)
    annotated.attrs["thresholds"] = legacy.attrs["thresholds"]
    annotated.attrs["loo_max"] = legacy.attrs["loo_max"]
    return annotated


def _annotate_literature(scored: pd.DataFrame) -> pd.DataFrame:
    """Attach literature NT when gt_data.csv is present — novel candidates get a check."""
    scored = scored.copy()
    scored["gt_verified_nt"] = pd.NA
    scored["gt_agrees_with_pattern"] = pd.NA
    if not GT_DATA.exists():
        return scored

    from name_matching import build_match_index, find_match

    gt = pd.read_csv(GT_DATA)
    if "species" in gt.columns:
        gt = gt[gt["species"] == "adult_drosophila_melanogaster"]
    if gt.empty or "cell_type" not in gt.columns:
        return scored

    nt_map = {
        "acetylcholine": "ACH",
        "glutamate": "GLUT",
        "gaba": "GABA",
        "dopamine": "DA",
        "serotonin": "SER",
        "octopamine": "OCT",
        "histamine": "HIST",
    }
    present_cols = [c for c in nt_map if c in gt.columns]

    def combine(group: pd.DataFrame) -> str:
        found = []
        for col in present_cols:
            if (group[col] == 1).any():
                found.append(nt_map[col])
        return ",".join(sorted(found))

    combined = {}
    for cell_type, group in gt.groupby("cell_type"):
        combined[cell_type] = combine(group)
    combined = pd.Series(combined)
    canon_idx, range_idx = build_match_index(list(combined.index))

    pattern_expected = {
        "histamine_blindspot": {"HIST"},
        "ORN_SER_confusion": {"ACH"},
        "Dm_GLUT_confusion": {"GLUT"},
    }

    verified = []
    agrees = []
    for name, pattern in zip(scored["cell_type"], scored["best_pattern"]):
        matched, _ = find_match(name, canon_idx, range_idx)
        if matched is None:
            verified.append(pd.NA)
            agrees.append(pd.NA)
            continue
        nts = set(part for part in str(combined.loc[matched]).split(",") if part)
        verified.append(",".join(sorted(nts)) if nts else pd.NA)
        expected = pattern_expected.get(pattern, set())
        agrees.append(bool(nts & expected) if nts else pd.NA)
    scored["gt_verified_nt"] = verified
    scored["gt_agrees_with_pattern"] = agrees
    return scored


def recovery_report(scored: pd.DataFrame) -> dict:
    """Did we recover the literature-confirmed seeds -- via simplex
    neighborhood, the reconstructed entropy z-score/FDR channel, or both?

    See module docstring: R1-6 (and, on real data, Dm12/Dm1) are only
    reachable via simplex geometry; R7 is the reverse -- geometrically an
    outlier within its own family, but a large, unambiguous entropy z-score
    outlier already caught by the existing channel. A seed counts as
    recovered if *either* channel independently catches it; `recovered_via`
    reports which one(s), so nothing recovers silently.
    """
    report = {}
    for pattern, seeds in PATTERN_SEEDS.items():
        present = [s for s in seeds if s in set(scored["cell_type"])]
        rows = scored[scored["cell_type"].isin(present)]

        recovered_via: dict[str, str] = {}
        for _, r in rows.iterrows():
            via = []
            if pattern in str(r["matched_patterns"]).split(","):
                via.append("simplex")
            if bool(r["entropy_channel_significant"]):
                via.append("entropy")
            if via:
                recovered_via[r["cell_type"]] = "+".join(via)

        report[pattern] = {
            "n_seeds": len(present),
            "n_recovered": len(recovered_via),
            "seeds": present,
            "recovered": list(recovered_via.keys()),
            "recovered_via": recovered_via,
            "missed": sorted(set(present) - set(recovered_via)),
        }
    return report


def run_scan() -> pd.DataFrame:
    ensure_output_dirs()
    df = load_entropy_table()
    print(f"Loaded {len(df)} cell types from entropy table")
    scored = score_types(df)

    thresholds = scored.attrs.get("thresholds", {})
    loo_max = scored.attrs.get("loo_max", {})
    print("\n[legacy heuristic, kept for comparison only] LOO max and thresholds:")
    for pattern in PATTERN_SEEDS:
        print(
            f"  {pattern}: LOO max = {loo_max.get(pattern, float('nan')):.4f}  "
            f"threshold = {thresholds.get(pattern, float('nan')):.4f}"
        )
    n_heuristic = int(scored["in_neighborhood_heuristic"].sum())
    n_calibrated = int(scored["in_neighborhood"].sum())
    print(
        f"\n'In neighborhood' under legacy heuristic: {n_heuristic}/{len(scored)} "
        f"({n_heuristic / len(scored):.0%})"
    )
    print(
        f"'In neighborhood' under calibrated permutation test (q<{calib.ALPHA}): "
        f"{n_calibrated}/{len(scored)} ({n_calibrated / len(scored):.0%})"
    )

    report = recovery_report(scored)
    print("\nSeed recovery (simplex neighborhood and/or reconstructed entropy z/FDR channel):")
    for pattern, info in report.items():
        print(f"  {pattern}: {info['n_recovered']}/{info['n_seeds']}")
        for seed, via in info["recovered_via"].items():
            print(f"    {seed}: recovered via {via}")
        if info["missed"]:
            print(f"    NOT recovered by either channel: {info['missed']}")

    novel = scored[scored["is_novel_candidate"]].copy()
    print(f"\nNovel simplex-neighborhood candidates (calibrated, not in name-matched list): {len(novel)}")
    if len(novel):
        cols = [
            "cell_type", "n_neurons", "entropy", "dominant_nt",
            "best_pattern", "best_js", "best_q_calibrated",
        ]
        if "gt_verified_nt" in novel.columns:
            cols.append("gt_verified_nt")
        print(novel[cols].head(30).to_string(index=False))

    out = RESULTS / "signature_scan.csv"
    scored.to_csv(out, index=False)
    print(f"\nSaved {out}")

    novel_out = RESULTS / "signature_scan_novel.csv"
    novel.to_csv(novel_out, index=False)
    print(f"Saved {novel_out}")

    # Compact summary for README / canvas
    summary_rows = [
        {"metric": "cell_types_scored", "value": len(scored)},
        {"metric": "seeds_recovered", "value": sum(info["n_recovered"] for info in report.values())},
        {"metric": "seeds_total", "value": sum(info["n_seeds"] for info in report.values())},
        {"metric": "novel_candidates_calibrated", "value": int(scored["is_novel_candidate"].sum())},
        {"metric": "novel_candidates_legacy_heuristic", "value": int(
            (scored["in_neighborhood_heuristic"] & ~scored["is_seed"] & ~scored["already_name_matched"]).sum()
        )},
        {
            "metric": "r16_js_histamine",
            "value": float(scored.loc[scored["cell_type"] == "R1-6", "js_histamine_blindspot"].iloc[0])
            if "R1-6" in set(scored["cell_type"]) else np.nan,
        },
        {
            "metric": "r16_q_histamine_calibrated",
            "value": float(scored.loc[scored["cell_type"] == "R1-6", "q_histamine_blindspot"].iloc[0])
            if "R1-6" in set(scored["cell_type"]) else np.nan,
        },
    ]
    pd.DataFrame(summary_rows).to_csv(RESULTS / "signature_scan_summary.csv", index=False)
    return scored


def main() -> None:
    parser = argparse.ArgumentParser(description="NT confusion signature scan")
    parser.parse_args()
    run_scan()


if __name__ == "__main__":
    main()
