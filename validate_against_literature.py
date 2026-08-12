"""
Validate flagged cell types against the literature ground truth
(flyconnectome/drosophila_neurotransmitters), and produce the deliverable
Arie asked for: a list of proposed corrections, one per dataset.

Ground truth semantics (gt_data.csv):
  1  = neurotransmitter verified present
 -1  = neurotransmitter verified absent
  0  = not assessed / no data

A cell type's "verified NT set" is every column marked 1. Most cell types
have exactly one; a genuine co-transmitter can have more than one (this is
represented in the data, e.g. some types show both acetylcholine=1 and
another transmitter=1).
"""
import pandas as pd
import numpy as np
from name_matching import build_match_index, find_match
from paths import GT_DATA, LITERATURE_VALIDATED, THREE_PATTERNS, ensure_output_dirs

ensure_output_dirs()

NT_COLUMNS = ['acetylcholine', 'glutamate', 'gaba', 'glycine', 'dopamine',
              'serotonin', 'octopamine', 'tyramine', 'histamine', 'nitric_oxide']

# Map ground truth column names -> FlyWire classifier NT codes
GT_TO_CODE = {
    'acetylcholine': 'ACH', 'glutamate': 'GLUT', 'gaba': 'GABA',
    'dopamine': 'DA', 'serotonin': 'SER', 'octopamine': 'OCT',
    'histamine': 'HIST', 'glycine': 'GLY', 'tyramine': 'TYR',
    'nitric_oxide': 'NO',
}

print("Loading ground truth...")
gt = pd.read_csv(GT_DATA)
gt = gt[gt['species'] == 'adult_drosophila_melanogaster'].copy()
print(f"  {len(gt)} adult D. melanogaster rows, {gt['cell_type'].nunique()} unique cell types")

# Some cell types have multiple ground-truth rows (different sources). Combine:
# a transmitter counts as verified-present for the type if ANY source marks it 1,
# unless another equally-or-more confident source marks it -1 with no support for 1.
def combine_gt_rows(group):
    verified_present = set()
    verified_absent = set()
    max_confidence = group['neurotransmitter_verified_confidence'].max()
    sources = '; '.join(group['neurotransmitter_verified_source'].unique())
    for nt_col in NT_COLUMNS:
        if (group[nt_col] == 1).any():
            verified_present.add(GT_TO_CODE[nt_col])
        elif (group[nt_col] == -1).all() and len(group) > 0:
            verified_absent.add(GT_TO_CODE[nt_col])
    return pd.Series({
        'verified_present': verified_present,
        'verified_absent': verified_absent,
        'max_confidence': max_confidence,
        'sources': sources,
        'n_gt_rows': len(group),
    })

gt_combined = gt.groupby('cell_type').apply(combine_gt_rows, include_groups=False).reset_index()
print(f"  Combined to {len(gt_combined)} unique cell types with verified NT calls")

gt_names = gt_combined['cell_type'].tolist()
gt_canon_idx, gt_range_idx = build_match_index(gt_names)
gt_lookup = gt_combined.set_index('cell_type').to_dict('index')


def match_against_gt(cell_type_name):
    matched, method = find_match(cell_type_name, gt_canon_idx, gt_range_idx)
    if matched is None:
        return None
    info = gt_lookup[matched]
    return {
        'gt_matched_name': matched,
        'gt_match_method': method,
        'gt_verified_present': info['verified_present'],
        'gt_verified_absent': info['verified_absent'],
        'gt_confidence': info['max_confidence'],
        'gt_sources': info['sources'],
    }


# ─────────────────────────────────────────────
# Load all flagged candidates from prior analysis (both n>=20 and n>=10 runs)
# ─────────────────────────────────────────────

print("\nLoading previously flagged candidates...")
pattern_flagged = pd.read_csv(THREE_PATTERNS)
print(f"  {len(pattern_flagged)} previously flagged FAFB cell types (histamine blindspot + ORN + Dm patterns)")

# ─────────────────────────────────────────────
# Cross-check every flagged type against literature ground truth
# ─────────────────────────────────────────────

print("\nCross-checking each flagged type against literature ground truth...")
records = []
for _, row in pattern_flagged.iterrows():
    ct = row['fafb_cell_type']
    gt_match = match_against_gt(ct)

    record = {
        'fafb_cell_type': ct,
        'pattern': row['pattern'],
        'fafb_dominant_nt_wrong': row['fafb_dominant_nt'],
        'mcns_says': row.get('mcns_confirmed_nt', row.get('mcns_majority_nt')),
    }

    if gt_match is None:
        record['literature_status'] = 'no_match_found'
        record['literature_verified_nt'] = None
        record['agrees_with_literature'] = None
    else:
        record["gt_matched_name"] = gt_match["gt_matched_name"]
        record["gt_match_method"] = gt_match["gt_match_method"]
        record["gt_confidence"] = gt_match["gt_confidence"]
        record["gt_sources"] = gt_match["gt_sources"]
        record["gt_verified_present"] = ",".join(sorted(gt_match["gt_verified_present"]))
        record["gt_verified_absent"] = ",".join(sorted(gt_match["gt_verified_absent"]))
        present = gt_match["gt_verified_present"]
        if len(present) == 0:
            record["literature_status"] = "matched_no_verified_nt"
            record["literature_verified_nt"] = None
            record["agrees_with_literature"] = None
        else:
            record["literature_status"] = "matched_with_verified_nt"
            record["literature_verified_nt"] = ",".join(sorted(present))
            mcns_or_pattern_claim = record["mcns_says"]
            record["agrees_with_literature"] = mcns_or_pattern_claim in present

    records.append(record)

results = pd.DataFrame(records)
print("\n" + "=" * 100)
print(results[['fafb_cell_type', 'pattern', 'fafb_dominant_nt_wrong', 'mcns_says',
                'literature_status', 'literature_verified_nt', 'agrees_with_literature',
                'gt_confidence']].to_string(index=False))

results.to_csv(LITERATURE_VALIDATED, index=False)
print(f"\nSaved {LITERATURE_VALIDATED}")

# Flag any disagreements for manual review
disagreements = results[results['agrees_with_literature'] == False]
if len(disagreements) > 0:
    print(f"\n*** WARNING: {len(disagreements)} candidate(s) where MCNS/pattern claim "
          f"DISAGREES with literature ground truth -- exclude from corrections list ***")
    print(disagreements[['fafb_cell_type', 'mcns_says', 'literature_verified_nt', 'gt_sources']].to_string(index=False))

no_match = results[results['literature_status'] == 'no_match_found']
print(f"\n{len(no_match)} candidate(s) with no literature match found (not confirmed, not contradicted):")
print(no_match['fafb_cell_type'].tolist())
