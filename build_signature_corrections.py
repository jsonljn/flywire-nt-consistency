"""Validate signature-scan candidates against literature evidence.

A correction requires a verified transmitter set that excludes the dominant
FAFB prediction. Agreement with a pattern label alone is insufficient.
Results are stored separately from name-matched corrections to retain
method provenance.
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
    if len(confirmed):
        confirmed["needs_correction"] = confirmed.apply(
            lambda r: prediction_needs_correction(r["dominant_nt"], r["verified_set"]), axis=1
        )
    else:
        # Preserve the output schema when no rows have literature support.
        confirmed["needs_correction"] = pd.Series(dtype=bool)

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
