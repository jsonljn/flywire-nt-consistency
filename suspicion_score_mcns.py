"""Rank MCNS corrections by confidence in the flagged prediction.

The MCNS export lacks per-category probabilities, so scores use E1 only
and record score_type as E1_only.
Requires data_mcns/merged_annotations.csv and corrections/corrections_mcns.csv.
Writes corrections/corrections_mcns_scored.csv.
"""
import pandas as pd


def main():
    print("Loading MCNS annotations and correction candidates...")
    mcns = pd.read_csv('data_mcns/merged_annotations.csv')
    corrections = pd.read_csv('corrections/corrections_mcns.csv')
    print(f"  {len(corrections)} MCNS correction candidates to score")

    mcns_scores = mcns[['root_id', 'nt_type_score']]
    scored = corrections.merge(mcns_scores, on='root_id', how='left')

    missing = scored['nt_type_score'].isna().sum()
    if missing:
        print(f"  WARNING: {missing} neurons missing nt_type_score after merge")

    scored['E1_confident_wrongness'] = scored['nt_type_score']
    scored['E2_evidence_for_correct'] = pd.NA
    scored['suspicion_score'] = scored['E1_confident_wrongness']
    scored['score_type'] = 'E1_only (MCNS has no per-category probability breakdown)'

    scored = scored.sort_values('suspicion_score', ascending=False).reset_index(drop=True)
    scored['rank'] = scored.index + 1

    cols = ['rank', 'root_id', 'cell_type', 'pattern', 'current_predicted_nt', 'verified_nt',
            'E1_confident_wrongness', 'E2_evidence_for_correct', 'suspicion_score', 'score_type',
            'evidence_source', 'evidence_confidence', 'proposed_action']
    scored = scored[cols]

    scored.to_csv('corrections/corrections_mcns_scored.csv', index=False)
    print(f"\nSaved corrections/corrections_mcns_scored.csv ({len(scored)} rows)")
    print(scored[['rank', 'cell_type', 'current_predicted_nt', 'verified_nt',
                  'suspicion_score']].to_string(index=False))


if __name__ == "__main__":
    main()
