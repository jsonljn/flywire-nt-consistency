"""
Shared MCNS cell-type matching for cross-dataset comparisons.

Uses name_matching.py for safe exact/range matches. MCNS subtype splits
(R7y, R7p, ...) are handled via explicit groups — not a generic prefix rule,
which could conflate unrelated types (e.g. Dm1 vs Dm12).
"""
import re
from typing import Optional

from name_matching import build_match_index, find_match

# FAFB coarse type -> MCNS finer subtypes to aggregate
EXPLICIT_SUBTYPE_GROUPS: dict[str, list[str]] = {
    "R7": ["R7y", "R7p", "R7d", "R7_unclear", "ExR7"],
    "R8": ["R8y", "R8p", "R8d", "R8_unclear", "ExR8"],
}

# Safe MCNS subtype suffix: letters/underscore only, no digits (Dm1 != Dm12)
_SUBTYPE_SUFFIX = re.compile(r"^[a-z_]+$")


def build_mcns_nt_lookup(
    mcns_df,
    type_col: str = "primary_type",
    nt_col: str = "nt_type",
) -> dict[str, dict]:
    """Per MCNS cell type: n, majority_nt, majority_frac."""
    lookup: dict[str, dict] = {}
    for cell_type, group in mcns_df.groupby(type_col):
        counts = group[nt_col].value_counts()
        lookup[cell_type] = {
            "n": len(group),
            "majority_nt": counts.index[0],
            "majority_frac": counts.iloc[0] / len(group),
        }
    return lookup


def _safe_subtype_variants(fafb_type: str, mcns_type_names: set[str]) -> list[str]:
    """
    MCNS names that are fafb_type plus a lowercase suffix (e.g. Dm3 -> Dm3a).
    Rejects suffixes containing digits so Dm1 does not match Dm12.
    """
    prefix = fafb_type
    variants = []
    for name in mcns_type_names:
        if not name.startswith(prefix) or name == prefix:
            continue
        suffix = name[len(prefix) :]
        if suffix and _SUBTYPE_SUFFIX.match(suffix):
            variants.append(name)
    return sorted(variants)


def resolve_fafb_to_mcns(
    fafb_type: str,
    mcns_type_names: list[str],
) -> tuple[Optional[list[str]], Optional[str]]:
    """
    Map a FAFB cell type to one or more MCNS type names.

    Returns (mcns_names, method) or (None, None) if no safe match.
    """
    mcns_set = set(mcns_type_names)

    if fafb_type in EXPLICIT_SUBTYPE_GROUPS:
        names = [n for n in EXPLICIT_SUBTYPE_GROUPS[fafb_type] if n in mcns_set]
        if names:
            return names, "explicit_subtype_group"
        # Fall through — MCNS may use the coarse name (plain "R7") instead of subtypes

    if fafb_type in mcns_set:
        return [fafb_type], "exact"

    canon_idx, range_idx = build_match_index(mcns_type_names)
    matched, method = find_match(fafb_type, canon_idx, range_idx)
    if matched:
        return [matched], method

    variants = _safe_subtype_variants(fafb_type, mcns_set)
    if variants:
        return variants, "subtype_group"

    return None, None


def aggregate_mcns_stats(
    mcns_names: list[str],
    mcns_lookup: dict[str, dict],
) -> dict:
    """Weighted aggregate across one or more MCNS cell types."""
    rows = [mcns_lookup[n] for n in mcns_names if n in mcns_lookup]
    if not rows:
        return {"n": 0, "majority_nt": None, "majority_frac": 0.0}

    total_n = sum(r["n"] for r in rows)
    if len(rows) == 1:
        return rows[0].copy()

    # Weighted fraction per NT label across subtypes
    nt_weight: dict[str, float] = {}
    for row in rows:
        nt = row["majority_nt"]
        nt_weight[nt] = nt_weight.get(nt, 0.0) + row["majority_frac"] * row["n"]

    dominant_nt = max(nt_weight, key=nt_weight.get)
    dominant_frac = nt_weight[dominant_nt] / total_n
    all_same_nt = len(nt_weight) == 1
    all_high_frac = all(r["majority_frac"] >= 0.9 for r in rows)

    return {
        "n": total_n,
        "majority_nt": dominant_nt,
        "majority_frac": dominant_frac,
        "all_subtypes_consistent": all_same_nt and all_high_frac,
        "subtype_nts": set(r["majority_nt"] for r in rows),
    }
