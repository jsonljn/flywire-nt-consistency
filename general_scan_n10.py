"""
General scan at n>=10: for every FAFB cell type, check if MCNS confirms it as
highly NT-consistent (any single dominant transmitter, not just histamine),
and flag cases where FAFB shows meaningfully elevated entropy despite that
external confirmation of consistency.

This generalizes the histamine-blindspot check to catch:
1. Categorical blind spots (true NT not in FAFB's 6-category output, e.g. HIST)
2. Genuine classifier confusion on a predictable category (e.g. potentially PFGs/SER)
"""
import pandas as pd

from mcns_matching import (
    aggregate_mcns_stats,
    build_mcns_nt_lookup,
    resolve_fafb_to_mcns,
)
from paths import (
    ENTROPY_RAW_N10,
    GENERAL_SCAN_FLAGGED,
    GENERAL_SCAN_FULL,
    MCNS_MERGED,
    ensure_output_dirs,
)

ensure_output_dirs()

# Prefer the n>=10 entropy file; fall back to the primary n>=20 table.
if ENTROPY_RAW_N10.exists():
    entropy_path = ENTROPY_RAW_N10
else:
    from paths import ENTROPY_RAW
    entropy_path = ENTROPY_RAW
    print(f"NOTE: {ENTROPY_RAW_N10.name} not found, using {ENTROPY_RAW.name}")

fafb_entropy = pd.read_csv(entropy_path)
mcns = pd.read_csv(MCNS_MERGED)
mcns = mcns.dropna(subset=["primary_type", "nt_type"])

print(f"FAFB cell types to check (n>=10): {len(fafb_entropy)}")

mcns_lookup = build_mcns_nt_lookup(mcns)
mcns_type_names = mcns["primary_type"].dropna().unique().tolist()

# Report drop reasons for the four MCNS histamine comparison types.
# Literature validation separately determines which support corrections.
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
        mcns_nt = stats["majority_nt"]
        mcns_frac = stats["majority_frac"]
        mcns_n = stats["n"]
    else:
        mcns_nt = stats["majority_nt"]
        mcns_frac = stats["majority_frac"]
        mcns_n = stats["n"]

    if mcns_frac < 0.9 or mcns_n < 10:
        if watched:
            dropped.append((fafb_type, f"matched {mcns_names}, n={mcns_n}, frac={mcns_frac:.3f} (below 0.9/10 threshold)"))
        continue

    results.append({
        "fafb_cell_type": fafb_type,
        "fafb_n": row["n_neurons"],
        "fafb_entropy": row["entropy"],
        "fafb_dominant_nt": row["dominant_nt"],
        "fafb_dominant_frac": row["dominant_frac"],
        "mcns_match": ",".join(mcns_names),
        "mcns_match_method": method,
        "mcns_n": mcns_n,
        "mcns_confirmed_nt": mcns_nt,
        "mcns_confirmed_frac": mcns_frac,
        "agrees_with_mcns": row["dominant_nt"] == mcns_nt,
        "is_categorical_blindspot": mcns_nt == "HIST",
    })

df = pd.DataFrame(results)
print(f"\nFAFB types with an MCNS match confirming >=90% single-NT consistency: {len(df)}")

surviving_watched = set(df["fafb_cell_type"]) if len(df) else set()
missing_watched = WATCH_TYPES - surviving_watched
if missing_watched:
    print(f"\nWARNING: {len(missing_watched)} watched seed type(s) did not reach results:")
    for name, reason in dropped:
        print(f"  - {name}: {reason}")
    unexplained = missing_watched - {name for name, _ in dropped}
    for name in unexplained:
        print(f"  - {name}: not in fafb_entropy input at all (check {entropy_path.name})")

# Flag cases with meaningful FAFB entropy despite MCNS-confirmed consistency.
# Use both an absolute floor (0.3 bits) and the 90th percentile among matched types.
entropy_90th = df["fafb_entropy"].quantile(0.90) if len(df) else 0.3
entropy_threshold = max(0.3, entropy_90th)
print(f"\nEntropy flag threshold: {entropy_threshold:.3f} (max of 0.3 and 90th percentile)")

flagged = df[df["fafb_entropy"] > entropy_threshold].sort_values("fafb_entropy", ascending=False)
print(f"\nFlagged (FAFB entropy > {entropy_threshold:.3f} despite MCNS confirming consistency): {len(flagged)}")
print(flagged[[
    "fafb_cell_type", "fafb_n", "fafb_entropy", "fafb_dominant_nt",
    "mcns_confirmed_nt", "mcns_n", "is_categorical_blindspot",
]].to_string(index=False))

df.to_csv(GENERAL_SCAN_FULL, index=False)
flagged.to_csv(GENERAL_SCAN_FLAGGED, index=False)
print(f"\nSaved {GENERAL_SCAN_FULL} and {GENERAL_SCAN_FLAGGED}")
