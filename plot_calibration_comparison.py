"""Before/after summary figure for the signature_scan calibration fix."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paths import FIGURES, RESULTS, ensure_output_dirs

ensure_output_dirs()
scored = pd.read_csv(RESULTS / "signature_scan.csv")

n_total = len(scored)
n_heuristic = int(
    (scored["in_neighborhood_heuristic"] & ~scored["is_seed"] & ~scored["already_name_matched"]).sum()
)
n_calibrated = int(scored["is_novel_candidate"].sum())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1, 1.3]})

# ── Left: flagged fraction, before vs after ──
bars = ax1.bar(
    ["Original heuristic\n(LOO max x 1.35)", "Calibrated exact test\n(this fix)"],
    [n_heuristic / n_total * 100, n_calibrated / n_total * 100],
    color=["#d62728", "#2ca02c"],
    width=0.55,
)
for bar, n in zip(bars, [n_heuristic, n_calibrated]):
    ax1.annotate(
        f"{n}/{n_total}\n({n/n_total:.0%})",
        (bar.get_x() + bar.get_width() / 2, bar.get_height()),
        ha="center", va="bottom", fontsize=11, fontweight="bold",
    )
ax1.set_ylabel("FAFB cell types flagged as\n'novel confusion candidate' (%)")
ax1.set_title("The bug: 80% of all cell types\nwere being flagged")
ax1.set_ylim(0, 95)
ax1.spines[["top", "right"]].set_visible(False)

# ── Right: what survives the fix, literature-checked ──
# Uses the same precise logic as build_signature_corrections.py (agreement
# with the *actual* verified transmitter set, not just the pattern's
# expected direction) so this figure and that script's CSVs always tell the
# same story -- see that script's docstring for why "literature confirms the
# pattern's expected direction" alone is not sufficient (a type can sit near
# a confusion fingerprint while still being correctly predicted).
from nt_utils import parse_verified_nts, prediction_needs_correction  # noqa: E402

novel = scored[scored["is_novel_candidate"]].copy()
novel["verified_set"] = novel["gt_verified_nt"].map(parse_verified_nts)
has_lit = novel["verified_set"].map(bool)
confirmed = novel[has_lit].copy()
confirmed["needs_correction"] = confirmed.apply(
    lambda r: prediction_needs_correction(r["dominant_nt"], r["verified_set"]), axis=1
)

counts = pd.Series({
    "Genuine correction\n(literature disagrees\nwith FAFB)": int(confirmed["needs_correction"].sum()),
    "Already correct\n(literature agrees\nwith FAFB)": int((~confirmed["needs_correction"]).sum()),
    "No literature entry\n(unconfirmed)": int((~has_lit).sum()),
})
colors2 = ["#2ca02c", "#7f7f7f", "#bbbbbb"]
bars2 = ax2.bar(counts.index, counts.values, color=colors2, width=0.6)
for bar, n in zip(bars2, counts.values):
    if n > 0:
        ax2.annotate(str(n), (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                     ha="center", va="bottom", fontsize=11, fontweight="bold")
ax2.set_ylabel(f"Novel candidates (of {n_calibrated} total)")
ax2.set_title("What survives the fix, checked against\nreal literature ground truth")
ax2.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    "signature_scan.py: fixed-threshold heuristic vs. permutation-calibrated exact test",
    fontsize=12, y=1.02,
)
plt.tight_layout()
out = FIGURES / "calibration_before_after.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved {out}")
