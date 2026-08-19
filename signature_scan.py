"""
Confusion-signature scan.

Entropy is a scalar: it detects *that* a cell type is inconsistent, but not
*why*, and it systematically misses large types that are consistently assigned
the *wrong* transmitter (R1-6: 82% ACH, z = -16, yet canonically
histaminergic).

This module treats each cell type as a point on the 6-simplex of FAFB
classifier outputs and scores it against literature-confirmed confusion
fingerprints with leave-one-out Jensen-Shannon distance.

That recovers the original three patterns *without name-matching to MCNS*,
and surfaces additional types that sit in the same simplex neighborhoods —
candidates the name-matched scan could never see.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from nt_simplex import (
    NT_ORDER,
    counts_to_vector,
    js_divergence,
    normalize,
    parse_nt_distribution,
)
from paths import (
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

# JS divergence below this is "same neighborhood" after LOO calibration.
# Overridden at runtime by max LOO seed distance * margin if larger.
JS_FLOOR = 0.12
LOO_MARGIN = 1.35


def load_entropy_table() -> pd.DataFrame:
    """Prefer n>=10 coverage; fall back to the primary n>=20 screen."""
    path = ENTROPY_RAW_N10 if ENTROPY_RAW_N10.exists() else ENTROPY_RAW
    df = pd.read_csv(path)
    df["counts"] = df["nt_distribution"].map(parse_nt_distribution)
    vectors = df["counts"].map(lambda c: normalize(counts_to_vector(c)))
    df["p_vec"] = list(vectors)
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


def score_types(df: pd.DataFrame) -> pd.DataFrame:
    lookup = _vector_lookup(df)
    known_flagged: set[str] = set()
    if THREE_PATTERNS.exists():
        known_flagged = set(pd.read_csv(THREE_PATTERNS)["fafb_cell_type"].astype(str))

    lit_confirmed: set[str] = set()
    if LITERATURE_VALIDATED.exists():
        lit = pd.read_csv(LITERATURE_VALIDATED)
        lit_confirmed = set(
            lit.loc[lit["agrees_with_literature"] == True, "fafb_cell_type"].astype(str)
        )

    # Calibrate per-pattern threshold from leave-one-out nearest-seed distances.
    # Histamine fingerprints are multi-modal (R7 ≠ R1-6), so the mean prototype
    # would smear them; nearest-seed distance keeps the family intact.
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

    thresholds = {
        pattern: max(JS_FLOOR, loo_max[pattern] * LOO_MARGIN)
        for pattern in PATTERN_SEEDS
    }

    rows = []
    for _, row in df.iterrows():
        name = str(row["cell_type"])
        p = row["p_vec"]
        js_by_pattern = {
            pattern: nearest_seed_js(p, pattern, lookup, exclude=name)
            for pattern in PATTERN_SEEDS
        }

        best_pattern = min(js_by_pattern, key=js_by_pattern.get)
        best_js = js_by_pattern[best_pattern]
        threshold = thresholds[best_pattern]
        is_seed = any(name in seeds for seeds in PATTERN_SEEDS.values())
        in_neighborhood = bool(best_js <= threshold)

        p_fast = float(sum(p[NT_ORDER.index(nt)] for nt in ("ACH", "GABA", "GLUT")))
        p_ser = float(p[NT_ORDER.index("SER")])

        rows.append({
            "cell_type": name,
            "n_neurons": int(row["n_neurons"]),
            "entropy": float(row["entropy"]),
            "dominant_nt": row["dominant_nt"],
            "dominant_frac": float(row["dominant_frac"]),
            "p_ACH": float(p[0]),
            "p_GABA": float(p[1]),
            "p_GLUT": float(p[2]),
            "p_DA": float(p[3]),
            "p_SER": p_ser,
            "p_OCT": float(p[5]),
            "p_fast_transmitters": p_fast,
            "js_histamine_blindspot": js_by_pattern["histamine_blindspot"],
            "js_ORN_SER_confusion": js_by_pattern["ORN_SER_confusion"],
            "js_Dm_GLUT_confusion": js_by_pattern["Dm_GLUT_confusion"],
            "best_pattern": best_pattern,
            "best_js": best_js,
            "pattern_threshold": threshold,
            "in_neighborhood": in_neighborhood,
            "is_seed": is_seed,
            "already_name_matched": name in known_flagged,
            "literature_confirmed": name in lit_confirmed,
            "is_novel_candidate": in_neighborhood and not is_seed and name not in known_flagged,
        })

    result = pd.DataFrame(rows)
    result = result.sort_values(
        ["is_novel_candidate", "in_neighborhood", "best_js"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    result.attrs["thresholds"] = thresholds
    result.attrs["loo_max"] = loo_max
    annotated = _annotate_literature(result)
    annotated.attrs["thresholds"] = thresholds
    annotated.attrs["loo_max"] = loo_max
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
    """Did simplex geometry recover the literature-confirmed seeds?"""
    report = {}
    for pattern, seeds in PATTERN_SEEDS.items():
        present = [s for s in seeds if s in set(scored["cell_type"])]
        recovered = scored[
            scored["cell_type"].isin(present)
            & (scored["best_pattern"] == pattern)
            & scored["in_neighborhood"]
        ]
        report[pattern] = {
            "n_seeds": len(present),
            "n_recovered": int(len(recovered)),
            "seeds": present,
            "recovered": recovered["cell_type"].tolist(),
        }
    return report


def run_scan() -> pd.DataFrame:
    ensure_output_dirs()
    df = load_entropy_table()
    print(f"Loaded {len(df)} cell types from entropy table")
    scored = score_types(df)

    thresholds = scored.attrs.get("thresholds", {})
    loo_max = scored.attrs.get("loo_max", {})
    print("\nLeave-one-out seed JS (max) and match thresholds:")
    for pattern in PATTERN_SEEDS:
        print(
            f"  {pattern}: LOO max = {loo_max.get(pattern, float('nan')):.4f}  "
            f"threshold = {thresholds.get(pattern, float('nan')):.4f}"
        )

    report = recovery_report(scored)
    print("\nSeed recovery (assigned to own pattern, inside threshold):")
    for pattern, info in report.items():
        print(f"  {pattern}: {info['n_recovered']}/{info['n_seeds']}")
        missed = set(info["seeds"]) - set(info["recovered"])
        if missed:
            print(f"    missed: {sorted(missed)}")

    novel = scored[scored["is_novel_candidate"]].copy()
    print(f"\nNovel simplex-neighborhood candidates (not in name-matched list): {len(novel)}")
    if len(novel):
        cols = [
            "cell_type", "n_neurons", "entropy", "dominant_nt",
            "best_pattern", "best_js",
        ]
        print(novel[cols].head(25).to_string(index=False))

    out = RESULTS / "signature_scan.csv"
    scored.drop(columns=[], errors="ignore").to_csv(out, index=False)
    print(f"\nSaved {out}")

    novel_out = RESULTS / "signature_scan_novel.csv"
    novel.to_csv(novel_out, index=False)
    print(f"Saved {novel_out}")

    # Compact summary for README / canvas
    summary_rows = [
        {
            "metric": "cell_types_scored",
            "value": len(scored),
        },
        {
            "metric": "seeds_recovered",
            "value": sum(info["n_recovered"] for info in report.values()),
        },
        {
            "metric": "seeds_total",
            "value": sum(info["n_seeds"] for info in report.values()),
        },
        {
            "metric": "novel_candidates",
            "value": int(scored["is_novel_candidate"].sum()),
        },
        {
            "metric": "r16_js_histamine",
            "value": float(
                scored.loc[scored["cell_type"] == "R1-6", "js_histamine_blindspot"].iloc[0]
            )
            if "R1-6" in set(scored["cell_type"])
            else np.nan,
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
