"""
Build the deliverable Arie requested: a list of corrections, one per dataset.

For every FAFB cell type that is:
  (a) flagged as inconsistent by the entropy method, AND
  (b) confirmed by literature ground truth (not contradicted, not just MCNS-only)

...list every individual neuron whose current predicted NT does not match the
literature-verified transmitter(s), with the evidence attached.

Excludes:
  - Lai (contradicted by literature -- MCNS said HIST, literature says GLUT)
  - Dm16, Dm20, Dm6, Dm9 (no literature match found -- unconfirmed, listed separately)

Also checks MCNS itself for any individual neurons within these same cell types
whose prediction disagrees with the literature-verified transmitter, since a
correction list should cover the dataset it applies to, not just FAFB.
"""
import pandas as pd

from mcns_matching import EXPLICIT_SUBTYPE_GROUPS
from name_matching import build_match_index, find_match
from nt_utils import parse_verified_nts, prediction_needs_correction
from paths import (
    CORRECTIONS,
    FAFB_MERGED,
    LITERATURE_VALIDATED,
    MCNS_MERGED,
    ensure_output_dirs,
)

ensure_output_dirs()

CONFIRMED_CANDIDATES = pd.read_csv(LITERATURE_VALIDATED)
CONFIRMED_CANDIDATES = CONFIRMED_CANDIDATES[
    CONFIRMED_CANDIDATES["agrees_with_literature"] == True  # noqa: E712
].copy()

print(f"Confirmed cell types eligible for corrections: {len(CONFIRMED_CANDIDATES)}")
print(CONFIRMED_CANDIDATES[["fafb_cell_type", "literature_verified_nt", "gt_confidence"]].to_string(index=False))

verified_lookup = {
    row["fafb_cell_type"]: parse_verified_nts(row["literature_verified_nt"])
    for _, row in CONFIRMED_CANDIDATES.iterrows()
}

confirmed_types = list(verified_lookup.keys())

# ─────────────────────────────────────────────
# FAFB corrections
# ─────────────────────────────────────────────

print("\nLoading FAFB annotations...")
fafb = pd.read_csv(FAFB_MERGED)
fafb_subset = fafb[fafb["primary_type"].isin(confirmed_types) & fafb["nt_type"].notna()].copy()

fafb_corrections = []
for _, row in fafb_subset.iterrows():
    ct = row["primary_type"]
    verified_nts = verified_lookup[ct]
    current = row["nt_type"]
    if prediction_needs_correction(current, verified_nts):
        gt_row = CONFIRMED_CANDIDATES[CONFIRMED_CANDIDATES["fafb_cell_type"] == ct].iloc[0]
        fafb_corrections.append({
            "root_id": row["root_id"],
            "cell_type": ct,
            "current_predicted_nt": current,
            "verified_nt": ",".join(sorted(verified_nts)),
            "evidence_source": gt_row["gt_sources"],
            "evidence_confidence": gt_row["gt_confidence"],
            "pattern": gt_row["pattern"],
            "proposed_action": (
                f"Review for relabel: classifier predicted {current}, "
                f"literature verifies {'/'.join(sorted(verified_nts))} "
                f"(source: {gt_row['gt_sources']}, confidence {gt_row['gt_confidence']}/5)"
            ),
        })

fafb_corrections_df = pd.DataFrame(fafb_corrections)
print(f"\nFAFB: {len(fafb_subset)} neurons checked across {len(confirmed_types)} confirmed cell types")
print(f"FAFB: {len(fafb_corrections_df)} individual neuron corrections proposed")
if len(fafb_corrections_df) > 0:
    print(fafb_corrections_df["cell_type"].value_counts())

fafb_out = CORRECTIONS / "corrections_fafb.csv"
fafb_corrections_df.to_csv(fafb_out, index=False)
print(f"Saved {fafb_out}")

# ─────────────────────────────────────────────
# MCNS corrections
# ─────────────────────────────────────────────

print("\nLoading MCNS annotations...")
mcns = pd.read_csv(MCNS_MERGED)

mcns_type_names = mcns["primary_type"].dropna().unique().tolist()
mcns_canon_idx, mcns_range_idx = build_match_index(mcns_type_names)

mcns_corrections = []
mcns_types_checked = []
for fafb_ct in confirmed_types:
    verified_nts = verified_lookup[fafb_ct]
    gt_row = CONFIRMED_CANDIDATES[CONFIRMED_CANDIDATES["fafb_cell_type"] == fafb_ct].iloc[0]

    if fafb_ct in EXPLICIT_SUBTYPE_GROUPS:
        mcns_names_to_check = [n for n in EXPLICIT_SUBTYPE_GROUPS[fafb_ct] if n in mcns_type_names]
        mcns_types_checked.append((fafb_ct, ",".join(mcns_names_to_check)))
    else:
        matched_mcns_name, method = find_match(fafb_ct, mcns_canon_idx, mcns_range_idx)
        if matched_mcns_name is None:
            continue
        mcns_names_to_check = [matched_mcns_name]
        mcns_types_checked.append((fafb_ct, matched_mcns_name))

    subset = mcns[mcns["primary_type"].isin(mcns_names_to_check) & mcns["nt_type"].notna()]

    for _, row in subset.iterrows():
        current = row["nt_type"]
        if prediction_needs_correction(current, verified_nts):
            mcns_corrections.append({
                "root_id": row["root_id"],
                "cell_type": row["primary_type"],
                "matched_fafb_type": fafb_ct,
                "current_predicted_nt": current,
                "verified_nt": ",".join(sorted(verified_nts)),
                "evidence_source": gt_row["gt_sources"],
                "evidence_confidence": gt_row["gt_confidence"],
                "pattern": gt_row["pattern"],
                "proposed_action": (
                    f"Review for relabel: classifier predicted {current}, "
                    f"literature verifies {'/'.join(sorted(verified_nts))} "
                    f"(source: {gt_row['gt_sources']}, confidence {gt_row['gt_confidence']}/5)"
                ),
            })

print(f"\nMCNS: matched {len(mcns_types_checked)}/{len(confirmed_types)} confirmed FAFB types to an MCNS name")
for fafb_ct, mcns_name in mcns_types_checked:
    if fafb_ct != mcns_name:
        print(f"  {fafb_ct} (FAFB) -> {mcns_name} (MCNS)")

mcns_corrections_df = pd.DataFrame(mcns_corrections)
print(f"MCNS: {len(mcns_corrections_df)} individual neuron corrections proposed")
if len(mcns_corrections_df) > 0:
    print(mcns_corrections_df["cell_type"].value_counts())

mcns_out = CORRECTIONS / "corrections_mcns.csv"
mcns_corrections_df.to_csv(mcns_out, index=False)
print(f"Saved {mcns_out}")

# ─────────────────────────────────────────────
# Excluded / unconfirmed, for transparency
# ─────────────────────────────────────────────

all_flagged = pd.read_csv(LITERATURE_VALIDATED)
excluded = all_flagged[all_flagged["agrees_with_literature"] != True]
excluded_out = CORRECTIONS / "excluded_unconfirmed_candidates.csv"
excluded.to_csv(excluded_out, index=False)
print(f"\nExcluded/unconfirmed candidates (NOT included in corrections): {len(excluded)}")
print(excluded[["fafb_cell_type", "literature_status", "mcns_says", "literature_verified_nt"]].to_string(index=False))
print(f"Saved {excluded_out}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(
    f"FAFB corrections proposed: {len(fafb_corrections_df)} neurons across "
    f"{fafb_corrections_df['cell_type'].nunique() if len(fafb_corrections_df) else 0} cell types"
)
print(
    f"MCNS corrections proposed: {len(mcns_corrections_df)} neurons across "
    f"{mcns_corrections_df['cell_type'].nunique() if len(mcns_corrections_df) else 0} cell types"
)
print(f"Excluded (contradicted or unconfirmed by literature): {len(excluded)} cell types")
