"""
End-to-end pipeline runner for the FlyWire NT consistency project.

Usage:
    py -3.13 run_pipeline.py              # full pipeline (requires data/)
    py -3.13 run_pipeline.py --figures    # regenerate figures from existing results
    py -3.13 run_pipeline.py --validate   # validate README claims only
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYTHON = sys.executable


def run_step(label: str, script: str, *args: str) -> None:
    cmd = [PYTHON, str(ROOT / script), *args]
    print("\n" + "=" * 70)
    print(f"STEP: {label}")
    print(" ".join(cmd))
    print("=" * 70)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run FlyWire NT consistency pipeline")
    parser.add_argument("--figures-only", action="store_true", help="Regenerate figures only")
    parser.add_argument("--validate-only", action="store_true", help="Run validation checks only")
    parser.add_argument("--skip-connectivity", action="store_true", help="Skip connectivity analysis")
    args = parser.parse_args()

    if args.validate_only:
        run_step("Validate README claims", "validate_results.py")
        return

    if args.figures_only:
        for script in (
            "plot_histamine_blindspot.py",
            "plot_patterns.py",
            "connectivity_plots.py",
        ):
            run_step(f"Plot: {script}", script)
        return

    # Data prep
    if (ROOT / "data" / "neurons.csv").exists():
        run_step("Merge FAFB annotations", "merge_data.py")
    else:
        print("SKIP merge_data.py — data/neurons.csv not found")

    if (ROOT / "data_mcns" / "neurons.csv").exists():
        run_step("Normalize MCNS annotations", "normalize_mcns.py")
    else:
        print("SKIP normalize_mcns.py — data_mcns/neurons.csv not found")

    # Primary entropy screen (n>=20)
    if (ROOT / "data" / "merged_annotations.csv").exists():
        run_step("FAFB entropy analysis (n>=20)", "analysis.py", "data/merged_annotations.csv")
        run_step("Full histamine cross-dataset scan", "full_histamine_scan.py")
        run_step("Histamine pattern check", "histamine_pattern_check.py")
    else:
        print("SKIP FAFB analysis — data/merged_annotations.csv not found")

    # Sensitivity screen (n>=10)
    if (ROOT / "data" / "merged_annotations.csv").exists():
        run_step(
            "FAFB entropy analysis (n>=10)",
            "analysis.py",
            "data/merged_annotations.csv",
            "--min-members",
            "10",
        )
        run_step("General MCNS scan (n>=10)", "general_scan_n10.py")
        run_step("Three patterns summary", "three_patterns_summary.py")

    # Literature validation + corrections
    if (ROOT / "data" / "gt_data.csv").exists() and (ROOT / "results" / "three_confusion_patterns.csv").exists():
        run_step("Literature validation", "validate_against_literature.py")
        run_step("Build correction lists", "build_corrections.py")
    else:
        print("SKIP literature validation — gt_data.csv or three_confusion_patterns.csv missing")

    # Connectivity
    if not args.skip_connectivity and (ROOT / "data" / "connections.csv").exists():
        run_step("Connectivity comparison", "connectivity_comparison.py")
        run_step("Connectivity PCA plots", "connectivity_plots.py")

    # Figures + validation
    run_step("Histamine blindspot figure", "plot_histamine_blindspot.py")
    if (ROOT / "results" / "three_confusion_patterns.csv").exists():
        run_step("Pattern summary figures", "plot_patterns.py")
    run_step("Validate README claims", "validate_results.py")

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
