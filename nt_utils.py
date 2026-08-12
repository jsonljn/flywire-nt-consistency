"""Neurotransmitter verification helpers shared across correction scripts."""
from __future__ import annotations


def parse_verified_nts(value: str | None) -> set[str]:
    """Parse comma-separated verified NT codes from literature validation output."""
    if value is None or (isinstance(value, float) and str(value) == "nan"):
        return set()
    text = str(value).strip()
    if not text:
        return set()
    return {part.strip() for part in text.split(",") if part.strip()}


def prediction_needs_correction(current: str, verified_nts: set[str]) -> bool:
    """
    Return True if the classifier prediction should be flagged for review.

    Co-transmitter types (e.g. R8: ACH + HIST) only flag predictions outside the
    verified set — ACH-predicted R8 neurons are left alone per README logic.
    """
    if not verified_nts:
        return False
    return current not in verified_nts
