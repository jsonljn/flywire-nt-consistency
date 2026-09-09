"""Rank FAFB correction candidates by classifier confidence and transmitter support.

E1 is nt_type_score, the confidence in the flagged prediction.
E2 is the probability assigned to the literature-verified transmitter set.
E2 is undefined when the verified transmitter is outside FAFB's categories.

    suspicion = E1 * (1 - E2)  when E2 is defined
    suspicion = E1             otherwise

Scores prioritize review and are not calibrated probabilities of error.
Requires data/merged_annotations.csv and corrections/corrections_fafb.csv.
Writes corrections/corrections_fafb_scored.csv.
"""
import pandas as pd
import numpy as np

AVG_COL_FOR_NT = {
    'ACH': 'ach_avg', 'GLUT': 'glut_avg', 'GABA': 'gaba_avg',
    'DA': 'da_avg', 'SER': 'ser_avg', 'OCT': 'oct_avg',
}


def compute_e2(row):
    verified = str(row['verified_nt']).split(',')
    real_cats = [nt for nt in verified if nt in AVG_COL_FOR_NT]
    if not real_cats:
        return np.nan
    return max(row[AVG_COL_FOR_NT[nt]] for nt in real_cats)


def main():
    print("Loading FAFB annotations and correction candidates...")
    fafb = pd.read_csv('data/merged_annotations.csv')
    corrections = pd.read_csv('corrections/corrections_fafb.csv')
    print(f"  {len(corrections)} correction candidates to score")

    fafb_scores = fafb[['root_id', 'nt_type_score', 'ach_avg', 'glut_avg', 'gaba_avg',
                         'da_avg', 'ser_avg', 'oct_avg']]
    scored = corrections.merge(fafb_scores, on='root_id', how='left')

    missing = scored['nt_type_score'].isna().sum()
    if missing:
        print(f"  WARNING: {missing} neurons missing nt_type_score after merge")

    scored['E1_confident_wrongness'] = scored['nt_type_score']
    if len(scored):
        scored['E2_evidence_for_correct'] = scored.apply(compute_e2, axis=1)
        scored['suspicion_score'] = scored.apply(
            lambda r: r['E1_confident_wrongness'] if pd.isna(r['E2_evidence_for_correct'])
            else r['E1_confident_wrongness'] * (1 - r['E2_evidence_for_correct']),
            axis=1,
        )
    else:
        # Preserve the output schema when no corrections exist.
        scored['E2_evidence_for_correct'] = pd.Series(dtype=float)
        scored['suspicion_score'] = pd.Series(dtype=float)
    scored['score_type'] = np.where(
        scored['E2_evidence_for_correct'].isna(),
        'categorical_blindspot (E1 only)', 'full (E1 x E2)',
    )

    scored = scored.sort_values('suspicion_score', ascending=False).reset_index(drop=True)
    scored['rank'] = scored.index + 1

    cols = ['rank', 'root_id', 'cell_type', 'pattern', 'current_predicted_nt', 'verified_nt',
            'E1_confident_wrongness', 'E2_evidence_for_correct', 'suspicion_score', 'score_type',
            'evidence_source', 'evidence_confidence', 'proposed_action']
    scored = scored[cols]

    scored.to_csv('corrections/corrections_fafb_scored.csv', index=False)
    print(f"\nSaved corrections/corrections_fafb_scored.csv ({len(scored)} rows)")
    print("\nTop 10 highest-priority corrections:")
    print(scored.head(10)[['rank', 'cell_type', 'current_predicted_nt', 'verified_nt',
                            'suspicion_score']].to_string(index=False))
    print("\nScore type breakdown:")
    print(scored['score_type'].value_counts())


if __name__ == "__main__":
    main()
