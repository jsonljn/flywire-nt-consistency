"""
Test whether the R7/R8 pattern (spurious entropy from a missing HIST category)
generalizes to other known-histaminergic cell types in FAFB.

Strategy:
1. From MCNS (whose classifier includes histamine), find cell types that are
   confidently histaminergic (majority NT = HIST).
2. Match those cell type names to FAFB's cell types via mcns_matching.
3. Check FAFB's raw entropy for those matched types -- if the hypothesis is
   right, they should show unusually high entropy despite being a single,
   real, consistent (histaminergic) population.
"""
import pandas as pd

from mcns_matching import EXPLICIT_SUBTYPE_GROUPS
from name_matching import build_match_index, find_match
from paths import (
    CONFIRMED_HISTAMINERGIC,
    ENTROPY_RAW,
    FAFB_MERGED,
    MCNS_MERGED,
    ensure_output_dirs,
)

ensure_output_dirs()

mcns = pd.read_csv(MCNS_MERGED)
mcns = mcns.dropna(subset=["primary_type", "nt_type"])
fafb = pd.read_csv(FAFB_MERGED)
fafb = fafb.dropna(subset=["primary_type", "nt_type"])
fafb_entropy = pd.read_csv(ENTROPY_RAW)
fafb_all_types = set(fafb["primary_type"].unique())
fafb_entropy_types = set(fafb_entropy["cell_type"].unique())

print("=" * 70)
print("STEP 1: Find MCNS cell types that are confidently histaminergic")
print("=" * 70)

rows = []
for cell_type, group in mcns.groupby("primary_type"):
    counts = group["nt_type"].value_counts()
    rows.append({
        "mcns_cell_type": cell_type,
        "n_neurons": len(group),
        "majority_nt": counts.index[0],
        "majority_frac": counts.iloc[0] / len(group),
    })
mcns_types = pd.DataFrame(rows)

hist_types = mcns_types[
    (mcns_types["majority_nt"] == "HIST")
    & (mcns_types["majority_frac"] >= 0.9)
    & (mcns_types["n_neurons"] >= 10)
].sort_values("n_neurons", ascending=False)

print(f"\nMCNS cell types that are >=90% histaminergic (n>=10): {len(hist_types)}")
print(hist_types.to_string(index=False))

print("\n" + "=" * 70)
print("STEP 2: Match these to FAFB cell types and check entropy")
print("=" * 70)

# Reverse map: MCNS subtype -> FAFB coarse type for explicit groups
mcns_to_fafb_explicit = {}
for fafb_type, subtypes in EXPLICIT_SUBTYPE_GROUPS.items():
    for subtype in subtypes:
        mcns_to_fafb_explicit[subtype] = fafb_type

fafb_canon_idx, fafb_range_idx = build_match_index(list(fafb_all_types))

matches = []
for _, row in hist_types.iterrows():
    mcns_ct = row["mcns_cell_type"]
    matched_fafb_type = None
    match_method = None

    if mcns_ct in mcns_to_fafb_explicit:
        matched_fafb_type = mcns_to_fafb_explicit[mcns_ct]
        match_method = "explicit_subtype_group"
    elif mcns_ct in fafb_all_types:
        matched_fafb_type = mcns_ct
        match_method = "exact"
    else:
        matched, method = find_match(mcns_ct, fafb_canon_idx, fafb_range_idx)
        if matched:
            matched_fafb_type = matched
            match_method = method

    if matched_fafb_type and matched_fafb_type in fafb_entropy_types:
        matches.append({
            "mcns_cell_type": mcns_ct,
            "mcns_n": row["n_neurons"],
            "mcns_hist_frac": row["majority_frac"],
            "fafb_cell_type": matched_fafb_type,
            "match_method": match_method,
        })

matches_df = pd.DataFrame(matches)
print(f"\nMatched to FAFB cell types: {len(matches_df)} MCNS types matched")
if len(matches_df) > 0:
    print(matches_df.to_string(index=False))

if len(matches_df) > 0:
    unique_fafb_types = matches_df["fafb_cell_type"].unique()
    print(f"\nUnique FAFB cell types confirmed histaminergic via MCNS: {len(unique_fafb_types)}")
    print(sorted(unique_fafb_types))

    print("\n" + "=" * 70)
    print("STEP 3: Entropy of these confirmed-histaminergic types in FAFB")
    print("=" * 70)

    result = fafb_entropy[fafb_entropy["cell_type"].isin(unique_fafb_types)].copy()
    result = result.sort_values("entropy", ascending=False)
    print(result[["cell_type", "n_neurons", "entropy", "dominant_nt", "dominant_frac"]].to_string(index=False))

    print(f"\nMean entropy of these types: {result['entropy'].mean():.3f}")
    print(f"Mean entropy of ALL FAFB cell types (for comparison): {fafb_entropy['entropy'].mean():.3f}")
    print(f"Fraction of these types with entropy > 0.5: {(result['entropy'] > 0.5).mean():.1%}")
    print(f"Fraction of ALL FAFB types with entropy > 0.5: {(fafb_entropy['entropy'] > 0.5).mean():.1%}")

    result.to_csv(CONFIRMED_HISTAMINERGIC, index=False)
    print(f"\nSaved {CONFIRMED_HISTAMINERGIC}")
else:
    print("\nNo matches found.")
