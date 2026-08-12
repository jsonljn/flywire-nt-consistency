"""
Validate key README claims against saved result files.

Run after the pipeline completes (or on existing results) to catch regressions.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from paths import (
    CONFIRMED_HISTAMINERGIC_SUMMARY,
    CONNECTIVITY_SUMMARY,
    CORRECTIONS,
    ENTROPY_CORRECTED,
    ENTROPY_RAW,
    LITERATURE_VALIDATED,
    RESULTS,
    THREE_PATTERNS,
    VALIDATION_REPORT,
    ensure_output_dirs,
)

ensure_output_dirs()

checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    checks.append((name, condition, detail))
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""))


print("=" * 70)
print("VALIDATION REPORT")
print("=" * 70)

# ── Entropy screen (n>=20) ──
if ENTROPY_CORRECTED.exists():
    ec = pd.read_csv(ENTROPY_CORRECTED)
    n_types = len(ec)
    sig_z2 = (ec["z_score"] > 2).sum()
    top = ec.sort_values("z_score", ascending=False).head(2)

    check("402 cell types at n>=20", n_types == 402, f"got {n_types}")
    check("Exactly 2 outliers at z>2", sig_z2 == 2, f"got {sig_z2}")
    if len(top) >= 2:
        check(
            "Top outliers are R7 and R8",
            set(top["cell_type"].head(2)) == {"R7", "R8"},
            f"got {list(top['cell_type'].head(2))}",
        )
        r7 = ec[ec["cell_type"] == "R7"].iloc[0]
        check(
            "R7 z-score ~7.8",
            7.5 < r7["z_score"] < 8.1,
            f"z={r7['z_score']:.2f}",
        )
else:
    check("entropy_corrected.csv exists", False, "file missing")

if ENTROPY_RAW.exists():
    er = pd.read_csv(ENTROPY_RAW)
    hist_in_fafb = er["dominant_nt"].eq("HIST").any() or (
        er["nt_distribution"].astype(str).str.contains("HIST").any()
        if "nt_distribution" in er.columns
        else False
    )
    check("FAFB predicts zero HIST at type level", not hist_in_fafb, "HIST found in distributions")

# ── Three patterns ──
if THREE_PATTERNS.exists():
    tp = pd.read_csv(THREE_PATTERNS)
    check("20 flagged cell types in three_patterns", len(tp) == 20, f"got {len(tp)}")
    check(
        "Pattern 1 includes R7, R8, Lai",
        set(tp[tp["pattern"] == "categorical_blindspot_HIST"]["fafb_cell_type"]) >= {"R7", "R8", "Lai"},
    )
    check(
        "Pattern 2 has 10 ORN types",
        (tp["pattern"] == "ORN_SER_confusion").sum() == 10,
    )
    check(
        "Pattern 3 has 7 Dm types",
        (tp["pattern"] == "Dm_GLUT_confusion").sum() == 7,
    )

# ── Literature validation ──
if LITERATURE_VALIDATED.exists():
    lv = pd.read_csv(LITERATURE_VALIDATED)
    confirmed = lv[lv["agrees_with_literature"] == True]
    excluded = lv[lv["agrees_with_literature"] == False]
    unconfirmed = lv[lv["literature_status"] == "no_match_found"]

    check("15 literature-confirmed types", len(confirmed) == 15, f"got {len(confirmed)}")
    check("1 literature-contradicted type (Lai)", len(excluded) == 1 and excluded.iloc[0]["fafb_cell_type"] == "Lai")
    check("4 unconfirmed types", len(unconfirmed) == 4, f"got {len(unconfirmed)}")

# ── Corrections ──
fafb_corr = CORRECTIONS / "corrections_fafb.csv"
mcns_corr = CORRECTIONS / "corrections_mcns.csv"
if fafb_corr.exists():
    fc = pd.read_csv(fafb_corr)
    check("1118 FAFB neuron corrections", len(fc) == 1118, f"got {len(fc)}")
    check("15 cell types in FAFB corrections", fc["cell_type"].nunique() == 15, f"got {fc['cell_type'].nunique()}")
    if "R8" in fc["cell_type"].values:
        r8_wrong = fc[fc["cell_type"] == "R8"]
        check(
            "R8 corrections exclude ACH (co-transmitter)",
            (r8_wrong["current_predicted_nt"] == "ACH").sum() == 0,
            f"{(r8_wrong['current_predicted_nt'] == 'ACH').sum()} ACH flagged",
        )
if mcns_corr.exists():
    mc = pd.read_csv(mcns_corr)
    check("4 MCNS neuron corrections", len(mc) == 4, f"got {len(mc)}")

# ── Connectivity ──
if CONNECTIVITY_SUMMARY.exists():
    cs = pd.read_csv(CONNECTIVITY_SUMMARY)
    r7 = cs[cs["cell_type"] == "R7"].iloc[0]
    r8 = cs[cs["cell_type"] == "R8"].iloc[0]
    check("R7 connectivity p < 0.001", r7["p_value"] < 0.001, f"p={r7['p_value']}")
    check("R8 connectivity p ~ 0.095", 0.05 < r8["p_value"] < 0.15, f"p={r8['p_value']}")

# ── Confirmed histaminergic summary ──
if CONFIRMED_HISTAMINERGIC_SUMMARY.exists():
    ch = pd.read_csv(CONFIRMED_HISTAMINERGIC_SUMMARY)
    check("4 confirmed histaminergic types", len(ch) == 4, f"got {len(ch)}")
    check("Includes R1-6", "R1-6" in ch["cell_type"].values)

# ── Write report ──
n_pass = sum(1 for _, ok, _ in checks if ok)
n_fail = sum(1 for _, ok, _ in checks if not ok)
summary = [
    "FlyWire NT Consistency — Validation Report",
    "=" * 50,
    f"Passed: {n_pass}/{len(checks)}",
    f"Failed: {n_fail}/{len(checks)}",
    "",
]
for name, ok, detail in checks:
    summary.append(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))

report_text = "\n".join(summary)
VALIDATION_REPORT.write_text(report_text, encoding="utf-8")
print("\n" + report_text)
print(f"\nSaved {VALIDATION_REPORT}")

if n_fail:
    raise SystemExit(1)
