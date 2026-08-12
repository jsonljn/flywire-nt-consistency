"""Unit tests for core analysis and matching logic."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis import shannon_entropy, stratified_permutation_null, compute_entropy_per_type
from mcns_matching import aggregate_mcns_stats, build_mcns_nt_lookup, resolve_fafb_to_mcns
from name_matching import build_match_index, find_match
from nt_utils import parse_verified_nts, prediction_needs_correction


class TestShannonEntropy:
    def test_pure_distribution_has_zero_entropy(self):
        assert shannon_entropy([100]) == pytest.approx(0.0)

    def test_uniform_binary_is_one_bit(self):
        assert shannon_entropy([50, 50]) == pytest.approx(1.0)

    def test_ignores_zero_counts(self):
        assert shannon_entropy([0, 10, 10]) == pytest.approx(1.0)


class TestNameMatching:
    @pytest.mark.parametrize(
        "query,targets,expected",
        [
            ("R1-6", ["R1-R6", "R7"], "R1-R6"),
            ("R1-R6", ["R1-6", "R7"], "R1-6"),
            ("r7", ["R7", "R8"], "R7"),
            ("XYZ", ["R7"], None),
        ],
    )
    def test_find_match(self, query, targets, expected):
        canon_idx, range_idx = build_match_index(targets)
        matched, _ = find_match(query, canon_idx, range_idx)
        assert matched == expected


class TestMcnsMatching:
    def test_range_notation(self):
        names = ["R1-R6", "R7y", "R7p", "Dm1", "Dm12", "Dm1a"]
        assert resolve_fafb_to_mcns("R1-6", names) == (["R1-R6"], "range_notation")

    def test_dm1_does_not_match_dm12(self):
        names = ["Dm1", "Dm12", "Dm1a"]
        assert resolve_fafb_to_mcns("Dm1", names) == (["Dm1"], "exact")
        assert resolve_fafb_to_mcns("Dm1", names)[0] != ["Dm12"]

    def test_r7_subtypes(self):
        names = ["R7y", "R7p", "R7d", "R8"]
        matched, method = resolve_fafb_to_mcns("R7", names)
        assert method == "explicit_subtype_group"
        assert set(matched) == {"R7y", "R7p", "R7d"}

    def test_r7_plain_fallback(self):
        names = ["R7", "R8"]
        assert resolve_fafb_to_mcns("R7", names) == (["R7"], "exact")

    def test_aggregate_subtypes(self):
        lookup = {
            "R7y": {"n": 100, "majority_nt": "HIST", "majority_frac": 1.0},
            "R7p": {"n": 80, "majority_nt": "HIST", "majority_frac": 0.95},
        }
        stats = aggregate_mcns_stats(["R7y", "R7p"], lookup)
        assert stats["majority_nt"] == "HIST"
        assert stats["n"] == 180
        assert stats["all_subtypes_consistent"] is True


class TestNtUtils:
    def test_parse_verified_nts(self):
        assert parse_verified_nts("ACH,HIST") == {"ACH", "HIST"}

    def test_r8_co_transmitter(self):
        verified = {"ACH", "HIST"}
        assert prediction_needs_correction("ACH", verified) is False
        assert prediction_needs_correction("GLUT", verified) is True


class TestStratifiedNull:
    def test_perfectly_consistent_type_has_low_z(self):
        df = pd.DataFrame({
            "primary_type": ["A"] * 50 + ["B"] * 50,
            "nt_type": ["ACH"] * 50 + ["GABA"] * 50,
        })
        result = stratified_permutation_null(
            df, "primary_type", "nt_type", min_members=20, n_permutations=200, seed=0
        )
        for cell_type in ("A", "B"):
            row = result[result["cell_type"] == cell_type].iloc[0]
            assert row["entropy"] == pytest.approx(0.0)
            assert row["z_score"] < 1.0

    def test_mixed_type_has_positive_z(self):
        df = pd.DataFrame({
            "primary_type": ["Mixed"] * 60,
            "nt_type": (["ACH"] * 20 + ["GABA"] * 20 + ["GLUT"] * 20),
        })
        result = stratified_permutation_null(
            df, "primary_type", "nt_type", min_members=20, n_permutations=300, seed=1
        )
        row = result.iloc[0]
        assert row["entropy"] > 1.0
        assert row["z_score"] > 0
