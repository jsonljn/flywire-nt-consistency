"""
Literature cross-check for signature_scan.py's novel candidates.

build_corrections.py already produces corrections/corrections_fafb.csv from
the *name-matched* three-pattern scan (three_confusion_patterns.csv ->
validate_against_literature.py -> literature_validated_candidates.csv). This
script does the equivalent last step for the *geometry-matched* candidates
signature_scan.py finds independently of any MCNS name match -- kept as a
separate, clearly-attributed output rather than merged into corrections_fafb.csv,
because these candidates were found by a different method (simplex distance +
an exact permutation test, not cross-dataset name matching) and deserve
independent provenance, in the same spirit as this project's existing
name_matching.py philosophy: don't quietly merge two different kinds of
evidence into one list.

A candidate only becomes an actual correction here if BOTH:
  1. gt_data.csv (literature) confirms a transmitter for it, AND
  2. that transmitter actually DISAGREES with FAFB's own dominant prediction
     (nt_utils.prediction_needs_correction) -- literature agreeing with the
     *pattern's expected transmitter* is not sufficient on its own, because a
     type can sit geometrically near a confusion fingerprint while still
     being correctly predicted (e.g. Dm4: literature GLUT, FAFB also predicts
     GLUT -- geometrically Dm-like, but there is nothing to correct).
"""
from __future__ import annotations

import pandas as pd

from nt_utils import parse_verified_nts, prediction_needs_correction
from paths import CORRECTIONS, RESULTS, ensure_output_dirs


def build_signature_scan_corrections() -> pd.DataFrame:
    ensure_output_dirs()
    scored = pd.read_csv(RESULTS / "signature_scan.csv")
    novel = scored[scored["is_novel_candidate"]].copy()

    novel["verified_set"] = novel["gt_verified_nt"].map(parse_verified_nts)
    has_lit = novel["verified_set"].map(bool)

    confirmed = novel[has_lit].copy()
    confirmed["needs_correction"] = confirmed.apply(
        lambda r: prediction_needs_correction(r["dominant_nt"], r["verified_set"]), axis=1
    )

    corrections = confirmed[confirmed["needs_correction"]].copy()
    corrections["corrected_nt"] = corrections["verified_set"].map(lambda s: ",".join(sorted(s)))
    cols = [
        "cell_type", "n_neurons", "dominant_nt", "corrected_nt", "best_pattern",
        "best_js", "best_p_calibrated", "gt_agrees_with_pattern",
    ]
    corrections = corrections[cols].sort_values("best_p_calibrated")

    already_correct = confirmed[~confirmed["needs_correction"]][
        ["cell_type", "n_neurons", "dominant_nt", "gt_verified_nt", "best_pattern", "best_p_calibrated"]
    ].sort_values("best_p_calibrated")

    unconfirmed = novel[~has_lit][
        ["cell_type", "n_neurons", "dominant_nt", "best_pattern", "best_p_calibrated"]
    ].sort_values("best_p_calibrated")

    out_dir = CORRECTIONS
    corrections.to_csv(out_dir / "corrections_signature_scan_novel.csv", index=False)
    already_correct.to_csv(out_dir / "signature_scan_novel_already_correct.csv", index=False)
    unconfirmed.to_csv(out_dir / "signature_scan_novel_unconfirmed.csv", index=False)

    print(f"Novel candidates from calibrated signature scan: {len(novel)}")
    print(f"  Literature-confirmed, genuinely mispredicted (new corrections): {len(corrections)}")
    if len(corrections):
        print(corrections.to_string(index=False))
    print(f"  Literature-confirmed but already correctly predicted (geometric near-miss): {len(already_correct)}")
    print(f"  No literature entry found (unconfirmed candidates): {len(unconfirmed)}")
    return corrections


if __name__ == "__main__":
    build_signature_scan_corrections()
