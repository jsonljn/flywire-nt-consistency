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
from nt_simplex import entropy_from_count_matrix, js_divergence, parse_nt_distribution
from nt_utils import parse_verified_nts, prediction_needs_correction
from signature_scan import recovery_report, score_types


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


class TestNtSimplex:
    def test_parse_legacy_numpy_repr(self):
        raw = "{'GLUT': np.int64(209), 'GABA': np.int64(184), 'ACH': np.int64(69)}"
        parsed = parse_nt_distribution(raw)
        assert parsed == {"GLUT": 209, "GABA": 184, "ACH": 69}

    def test_parse_json(self):
        assert parse_nt_distribution('{"ACH": 10, "GABA": 2}') == {"ACH": 10, "GABA": 2}

    def test_identical_distributions_have_zero_js(self):
        p = np.array([0.5, 0.5, 0, 0, 0, 0], dtype=float)
        assert js_divergence(p, p) == pytest.approx(0.0, abs=1e-12)

    def test_entropy_from_count_matrix_matches_scalar(self):
        counts = np.array([[10, 0, 0], [5, 5, 0], [1, 1, 1]], dtype=float)
        got = entropy_from_count_matrix(counts)
        expected = [shannon_entropy(row) for row in counts]
        np.testing.assert_allclose(got, expected, atol=1e-12)


class TestSignatureScan:
    def _tiny_entropy_frame(self) -> pd.DataFrame:
        """Minimal table with HIST/ORN/Dm seeds plus a nearby unlabeled type."""
        rows = [
            ("R7", 100, "{'GLUT': 44, 'GABA': 39, 'ACH': 15, 'DA': 1, 'SER': 1, 'OCT': 0}"),
            ("R8", 100, "{'ACH': 58, 'GLUT': 32, 'GABA': 9, 'SER': 1, 'OCT': 0, 'DA': 0}"),
            ("R1-6", 200, "{'ACH': 164, 'GLUT': 30, 'GABA': 4, 'SER': 2, 'OCT': 0, 'DA': 0}"),
            ("ORN_V", 40, "{'SER': 20, 'ACH': 18, 'GABA': 2, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DL3", 40, "{'SER': 40, 'ACH': 0, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DL4", 40, "{'SER': 40, 'ACH': 0, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DM3", 40, "{'SER': 39, 'ACH': 1, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DM2", 40, "{'SER': 38, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DA4l", 40, "{'SER': 37, 'ACH': 3, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DA4m", 40, "{'SER': 36, 'ACH': 4, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_DA3", 40, "{'SER': 34, 'ACH': 6, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_VA2", 40, "{'SER': 33, 'ACH': 7, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("ORN_VM3", 40, "{'SER': 22, 'ACH': 18, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("Dm12", 80, "{'GABA': 48, 'GLUT': 32, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Dm19", 20, "{'GABA': 14, 'GLUT': 6, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Dm1", 40, "{'GABA': 35, 'GLUT': 5, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("CleanACH", 50, "{'ACH': 50, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("NearORN", 30, "{'SER': 28, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
        ]
        records = []
        for name, n, dist in rows:
            parsed = parse_nt_distribution(dist)
            from nt_simplex import counts_to_vector, normalize, shannon_entropy_from_probs
            p = normalize(counts_to_vector(parsed))
            records.append({
                "cell_type": name,
                "n_neurons": n,
                "entropy": shannon_entropy_from_probs(p),
                "dominant_nt": max(parsed, key=parsed.get),
                "dominant_frac": max(parsed.values()) / n,
                "counts": parsed,
                "p_vec": p,
            })
        return pd.DataFrame(records)

    def test_recovers_histamine_seeds(self):
        scored = score_types(self._tiny_entropy_frame())
        report = recovery_report(scored)
        assert report["histamine_blindspot"]["n_recovered"] == report["histamine_blindspot"]["n_seeds"]

    def test_clean_type_is_not_a_novel_candidate(self):
        scored = score_types(self._tiny_entropy_frame())
        row = scored[scored["cell_type"] == "CleanACH"].iloc[0]
        assert row["is_novel_candidate"] in (False, 0)

    def test_nearby_orn_is_flagged(self):
        scored = score_types(self._tiny_entropy_frame())
        row = scored[scored["cell_type"] == "NearORN"].iloc[0]
        assert row["best_pattern"] == "ORN_SER_confusion"
        assert row["in_neighborhood"] in (True, 1)
