"""Project paths — single source of truth for data, results, and figures."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DATA_MCNS = ROOT / "data_mcns"
RESULTS = ROOT / "results"
CORRECTIONS = ROOT / "corrections"
FIGURES = ROOT / "figures"

FAFB_MERGED = DATA / "merged_annotations.csv"
MCNS_MERGED = DATA_MCNS / "merged_annotations.csv"
GT_DATA = DATA / "gt_data.csv"
FAFB_CONNECTIONS = DATA / "connections.csv"

# Primary screen (n >= 20) — keep legacy filenames for backward compatibility
ENTROPY_RAW = RESULTS / "entropy_raw.csv"
ENTROPY_CORRECTED = RESULTS / "entropy_corrected.csv"

# Sensitivity screen (n >= 10) — separate files so runs do not overwrite each other
ENTROPY_RAW_N10 = RESULTS / "entropy_raw_n10.csv"
ENTROPY_CORRECTED_N10 = RESULTS / "entropy_corrected_n10.csv"

GENERAL_SCAN_FULL = RESULTS / "general_scan_n10_full.csv"
GENERAL_SCAN_FLAGGED = RESULTS / "general_scan_n10_flagged.csv"
THREE_PATTERNS = RESULTS / "three_confusion_patterns.csv"
LITERATURE_VALIDATED = RESULTS / "literature_validated_candidates.csv"
FULL_CROSS_DATASET = RESULTS / "full_cross_dataset_scan.csv"
CONFIRMED_HISTAMINERGIC = RESULTS / "confirmed_histaminergic_in_fafb.csv"
CONFIRMED_HISTAMINERGIC_SUMMARY = RESULTS / "confirmed_histaminergic_summary.csv"
CONNECTIVITY_SUMMARY = RESULTS / "connectivity_test_summary.csv"
VALIDATION_REPORT = RESULTS / "validation_report.txt"


def entropy_paths(min_members: int = 20) -> tuple[Path, Path]:
    """Return (raw_entropy_csv, corrected_entropy_csv) for a given min-members threshold."""
    if min_members <= 10:
        return ENTROPY_RAW_N10, ENTROPY_CORRECTED_N10
    return ENTROPY_RAW, ENTROPY_CORRECTED


def ensure_output_dirs() -> None:
    """Create output directories if they do not exist."""
    for directory in (RESULTS, CORRECTIONS, FIGURES):
        directory.mkdir(exist_ok=True)
