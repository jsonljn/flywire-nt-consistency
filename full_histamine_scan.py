"""
Full scan: for every FAFB cell type (not just the 4 already confirmed),
check if a matching MCNS cell type exists and what MCNS's dominant NT is.

This flags any FAFB cell type where:
- FAFB entropy is elevated (potential inconsistency)
- MCNS confidently identifies the matching type as histaminergic

...as an additional histamine-blind-spot candidate beyond R7/R8/Lai/R1-6.
"""
import pandas as pd

from mcns_matching import (
    aggregate_mcns_stats,
    build_mcns_nt_lookup,
    resolve_fafb_to_mcns,
)
from paths import ENTROPY_RAW, FULL_CROSS_DATASET, MCNS_MERGED, ensure_output_dirs

ensure_output_dirs()

fafb_entropy = pd.read_csv(ENTROPY_RAW)
mcns = pd.read_csv(MCNS_MERGED)
mcns = mcns.dropna(subset=["primary_type", "nt_type"])

print(f"FAFB cell types to check: {len(fafb_entropy)}")

mcns_lookup = build_mcns_nt_lookup(mcns)
mcns_type_names = mcns["primary_type"].dropna().unique().tolist()
print(f"MCNS cell types available: {len(mcns_lookup)}")

# See general_scan_n10.py's WATCH_TYPES for why this matters: these are the
# project's own confirmed histaminergic seeds and must never drop silently.
WATCH_TYPES = {"R7", "R8", "R1-6", "Lai"}

results = []
dropped = []
for _, row in fafb_entropy.iterrows():
    fafb_type = row["cell_type"]
    watched = fafb_type in WATCH_TYPES

    mcns_names, method = resolve_fafb_to_mcns(fafb_type, mcns_type_names)
    if mcns_names is None:
        if watched:
            dropped.append((fafb_type, "no MCNS name match at all"))
        continue

    stats = aggregate_mcns_stats(mcns_names, mcns_lookup)
    if stats["n"] == 0:
        if watched:
            dropped.append((fafb_type, f"matched {mcns_names} but 0 rows in mcns_lookup"))
        continue

    if len(mcns_names) > 1:
        if not stats.get("all_subtypes_consistent", False):
            if watched:
                dropped.append((fafb_type, f"multi-subtype match {mcns_names} not all_subtypes_consistent"))
            continue
        mcns_nt = "HIST" if stats["majority_nt"] == "HIST" else stats["majority_nt"]
        mcns_frac = stats["majority_frac"]
        mcns_n = stats["n"]
    else:
        mcns_nt = stats["majority_nt"]
        mcns_frac = stats["majority_frac"]
        mcns_n = stats["n"]

    results.append({
        "fafb_cell_type": fafb_type,
        "fafb_n": row["n_neurons"],
        "fafb_entropy": row["entropy"],
        "fafb_dominant_nt": row["dominant_nt"],
        "mcns_match": ",".join(mcns_names),
        "match_method": method,
        "mcns_n": mcns_n,
        "mcns_majority_nt": mcns_nt,
        "mcns_majority_frac": mcns_frac,
    })

df = pd.DataFrame(results)
print(f"\nTotal FAFB types with any MCNS match: {len(df)}")

surviving_watched = set(df["fafb_cell_type"]) if len(df) else set()
missing_watched = WATCH_TYPES - surviving_watched
if missing_watched:
    print(f"\nWARNING: {len(missing_watched)} watched seed type(s) did not reach results:")
    for name, reason in dropped:
        print(f"  - {name}: {reason}")

hist_candidates = df[(df["mcns_majority_nt"] == "HIST") & (df["mcns_majority_frac"] >= 0.9)]
print(f"\nMCNS-confirmed histaminergic types with an FAFB match: {len(hist_candidates)}")
print(hist_candidates[[
    "fafb_cell_type", "fafb_n", "fafb_entropy", "fafb_dominant_nt",
    "mcns_n", "mcns_majority_frac",
]].sort_values("fafb_entropy", ascending=False).to_string(index=False))

df.to_csv(FULL_CROSS_DATASET, index=False)
print(f"\nSaved {FULL_CROSS_DATASET}")
