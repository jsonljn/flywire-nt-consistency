"""
Suspicion score for the MCNS correction candidates, extending the same
approach used for FAFB (suspicion_score.py).

MCNS's Neuron Attributes export exposes a single "Predicted NT confidence"
(nt_type_score) per neuron, not a per-category probability breakdown like
FAFB's ach_avg/glut_avg/etc. So only E1 (confidence in the wrong prediction)
is computable here -- E2 is not available for any MCNS correction, not just
the categorical blind-spot ones. This is a data-availability limit, not the
same structural reasoning as FAFB's histamine gap, and is documented as such
(score_type = "E1_only") rather than silently reusing the FAFB label.

Usage:
    python suspicion_score_mcns.py
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
