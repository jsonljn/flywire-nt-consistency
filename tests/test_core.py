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
from nt_simplex import (
    batch_js_divergence,
    counts_to_vector,
    entropy_from_count_matrix,
    js_divergence,
    normalize,
    parse_nt_distribution,
    shannon_entropy_from_probs,
)
from nt_utils import parse_verified_nts, prediction_needs_correction
from signature_calibration import exact_p_value
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

    def test_r1_6_matches_combined_range_name(self):
        combined_name = ["R1-R6", "R7y", "R7p"]
        assert resolve_fafb_to_mcns("R1-6", combined_name) == (["R1-R6"], "range_notation")

    def test_exr7_exr8_excluded_from_r7_r8_subtype_groups(self):
        from mcns_matching import EXPLICIT_SUBTYPE_GROUPS
        assert "ExR7" not in EXPLICIT_SUBTYPE_GROUPS["R7"]
        assert "ExR8" not in EXPLICIT_SUBTYPE_GROUPS["R8"]

        names = ["R7y", "R7p", "R7d", "R7_unclear", "ExR7"]
        mcns_df = pd.DataFrame({
            "primary_type": ["R7y"] * 66 + ["R7p"] * 78 + ["R7d"] * 68
            + ["R7_unclear"] * 75 + ["ExR7"] * 4,
            "nt_type": ["HIST"] * (66 + 78 + 68 + 75) + ["ACH"] * 4,
        })
        lookup = build_mcns_nt_lookup(mcns_df)
        matched, method = resolve_fafb_to_mcns("R7", names)
        assert matched == ["R7y", "R7p", "R7d", "R7_unclear"]
        stats = aggregate_mcns_stats(matched, lookup)
        assert stats.get("all_subtypes_consistent", False) is True
        assert stats["majority_nt"] == "HIST"

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

    def test_bincount_matches_add_at(self):
        rng = np.random.default_rng(3)
        n_types, n_nt, n = 50, 6, 5000
        type_codes = rng.integers(0, n_types, size=n)
        nt_codes = rng.integers(0, n_nt, size=n)

        via_add_at = np.zeros((n_types, n_nt), dtype=np.int64)
        np.add.at(via_add_at, (type_codes, nt_codes), 1)

        flat_idx = type_codes.astype(np.int64) * n_nt + nt_codes
        via_bincount = np.bincount(flat_idx, minlength=n_types * n_nt).reshape(n_types, n_nt)

        np.testing.assert_array_equal(via_add_at, via_bincount)


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

    def test_batch_js_divergence_matches_scalar(self):
        rng = np.random.default_rng(0)
        seed = np.array([0.5, 0.3, 0.1, 0.05, 0.03, 0.02])
        P = rng.dirichlet(np.ones(6), size=25)
        got = batch_js_divergence(P, seed)
        expected = [js_divergence(row, seed) for row in P]
        np.testing.assert_allclose(got, expected, atol=1e-10)


class TestExactPValue:
    """Test reference-pool probabilities and continuity correction."""

    def test_nothing_close_gives_high_p(self):
        pool = np.full(400, 0.9)
        p, i_obs, M = exact_p_value(observed=0.01, pool_distances=pool, k=3)
        assert i_obs == 0
        assert p < 0.02

    def test_everything_close_gives_p_near_one(self):
        pool = np.full(400, 0.001)
        p, i_obs, M = exact_p_value(observed=0.5, pool_distances=pool, k=3)
        assert i_obs == M
        assert p > 0.99

    def test_more_active_seeds_makes_a_match_less_surprising(self):
        pool = np.concatenate([np.full(5, 0.02), np.full(395, 0.9)])
        p_k1, _, _ = exact_p_value(observed=0.05, pool_distances=pool, k=1)
        p_k9, _, _ = exact_p_value(observed=0.05, pool_distances=pool, k=9)
        assert p_k9 > p_k1

    def test_continuity_correction_avoids_hard_zero_at_small_m(self):
        pool = np.array([0.9, 0.8, 0.7])
        p, i_obs, M = exact_p_value(observed=0.01, pool_distances=pool, k=3)
        assert i_obs == 0
        assert p > 0.0

    def test_correction_is_negligible_at_real_dataset_scale(self):
        pool = np.full(388, 0.9)
        p, _, _ = exact_p_value(observed=0.01, pool_distances=pool, k=3)
        uncorrected = 1.0 - (1.0 - 0 / 388) ** 3
        assert abs(p - uncorrected) < 0.01


class TestSignatureScan:
    def _large_background_rows(self, n: int = 100, seed: int = 0) -> list[tuple]:
        """Generate reference profiles with varied dominant categories and purity.

        With k=10, at least 97 reference types are required for the minimum
        continuity-corrected p-value to fall below 0.05.
        """
        rng = np.random.default_rng(seed)
        order = ["ACH", "GABA", "GLUT", "DA", "SER", "OCT"]
        rows = []
        for i in range(n):
            dominant = rng.integers(0, 6)
            purity = rng.uniform(0.7, 1.0)
            n_neurons = int(rng.integers(20, 200))
            remainder = 1.0 - purity
            fracs = rng.dirichlet(np.ones(6)) * remainder
            fracs[dominant] += purity
            counts = {order[j]: int(round(fracs[j] * n_neurons)) for j in range(6)}
            counts = {k: v for k, v in counts.items() if v > 0} or {order[dominant]: n_neurons}
            rows.append((f"Rand_{i:03d}", sum(counts.values()), str(counts)))
        return rows

    def _background_rows(self) -> list[tuple]:
        """Provide 24 reference profiles spanning the transmitter categories."""
        return [
            ("Bg_ACH1", 45, "{'ACH': 45, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_ACH2", 90, "{'ACH': 80, 'GABA': 6, 'GLUT': 4, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_ACH3", 120, "{'ACH': 100, 'GLUT': 15, 'GABA': 5, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_ACH4", 33, "{'ACH': 30, 'GABA': 2, 'GLUT': 1, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GABA1", 60, "{'GABA': 55, 'GLUT': 5, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GABA2", 70, "{'GABA': 60, 'ACH': 10, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GABA3", 25, "{'GABA': 25, 'ACH': 0, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GLUT1", 55, "{'GLUT': 50, 'GABA': 5, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GLUT2", 80, "{'GLUT': 70, 'ACH': 10, 'GABA': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GLUT3", 40, "{'GLUT': 40, 'GABA': 0, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_DA1", 28, "{'DA': 26, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_DA2", 50, "{'DA': 45, 'GLUT': 5, 'ACH': 0, 'GABA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_OCT1", 30, "{'OCT': 28, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'SER': 0}"),
            ("Bg_SER1", 26, "{'SER': 24, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
            ("Bg_Mix1", 65, "{'ACH': 25, 'GABA': 25, 'GLUT': 15, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_Mix2", 44, "{'GABA': 20, 'GLUT': 20, 'ACH': 4, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_Mix3", 38, "{'ACH': 18, 'GLUT': 18, 'GABA': 2, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_Mix4", 52, "{'ACH': 30, 'DA': 15, 'GABA': 7, 'GLUT': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_ACH5", 200, "{'ACH': 180, 'GLUT': 15, 'GABA': 5, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GABA4", 90, "{'GABA': 80, 'GLUT': 10, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_ACH6", 21, "{'ACH': 21, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GLUT4", 22, "{'GLUT': 22, 'ACH': 0, 'GABA': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_GABA5", 35, "{'GABA': 33, 'ACH': 2, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
            ("Bg_ACH7", 48, "{'ACH': 44, 'GABA': 4, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}"),
        ]

    def _frame_from_rows(self, rows: list[tuple]) -> pd.DataFrame:
        """rows: (name, n, dist_str_or_dict[, z_score]) tuples -> scoreable frame."""
        records = []
        for row in rows:
            name, n, dist = row[0], row[1], row[2]
            z = row[3] if len(row) > 3 else 0.0
            parsed = dist if isinstance(dist, dict) else parse_nt_distribution(dist)
            p = normalize(counts_to_vector(parsed))
            records.append({
                "cell_type": name,
                "n_neurons": n,
                "entropy": shannon_entropy_from_probs(p),
                "dominant_nt": max(parsed, key=parsed.get),
                "dominant_frac": max(parsed.values()) / n,
                "counts": parsed,
                "p_vec": p,
                "z_score": z,
            })
        return pd.DataFrame(records)

    def _tiny_entropy_frame(self) -> pd.DataFrame:
        """Build a fixture for patterns with two or three active seeds.

        Include positive R7/R8 and negative R1-6 entropy z-scores. ORN tests require
        the larger reference pool supplied by _large_background_rows.
        """
        rows = [
            ("R7", 100, "{'GLUT': 44, 'GABA': 39, 'ACH': 15, 'DA': 1, 'SER': 1, 'OCT': 0}", 7.8),
            ("R8", 100, "{'ACH': 58, 'GLUT': 32, 'GABA': 9, 'SER': 1, 'OCT': 0, 'DA': 0}", 3.5),
            ("R1-6", 200, "{'ACH': 164, 'GLUT': 30, 'GABA': 4, 'SER': 2, 'OCT': 0, 'DA': 0}", -16.8),
            ("ORN_V", 40, "{'SER': 20, 'ACH': 18, 'GABA': 2, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DL3", 40, "{'SER': 40, 'ACH': 0, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DL4", 40, "{'SER': 40, 'ACH': 0, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DM3", 40, "{'SER': 39, 'ACH': 1, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DM2", 40, "{'SER': 38, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DA4l", 40, "{'SER': 37, 'ACH': 3, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DA4m", 40, "{'SER': 36, 'ACH': 4, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_DA3", 40, "{'SER': 34, 'ACH': 6, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_VA2", 40, "{'SER': 33, 'ACH': 7, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("ORN_VM3", 40, "{'SER': 22, 'ACH': 18, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}", 0.0),
            ("Dm12", 80, "{'GABA': 48, 'GLUT': 32, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}", -1.8),
            ("Dm19", 20, "{'GABA': 14, 'GLUT': 6, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}", -1.5),
            ("Dm1", 40, "{'GABA': 35, 'GLUT': 5, 'ACH': 0, 'DA': 0, 'SER': 0, 'OCT': 0}", -2.9),
            ("CleanACH", 50, "{'ACH': 50, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}", 0.0),
            ("CleanGABA", 50, "{'GABA': 50, 'ACH': 0, 'GLUT': 0, 'DA': 0, 'SER': 0, 'OCT': 0}", 0.0),
        ]
        rows += [(name, n, dist, 0.0) for name, n, dist in self._background_rows()]
        return self._frame_from_rows(rows)

    def test_recovers_seeds_via_at_least_one_channel(self):
        """Check Dm recovery with a small reference pool.

        Histamine recovery depends on reference composition; ORN recovery requires
        a larger pool. Separate tests evaluate those cases.
        """
        scored = score_types(self._tiny_entropy_frame())
        report = recovery_report(scored)
        assert report["Dm_GLUT_confusion"]["n_recovered"] >= 2

    def test_orn_cluster_recovers_via_simplex(self):
        """Use 120 reference types to support significance testing with ten ORN seeds."""
        rows = [
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
        ]
        rows += self._large_background_rows(n=120, seed=0)
        scored = score_types(self._frame_from_rows(rows))
        report = recovery_report(scored)
        assert report["ORN_SER_confusion"]["n_recovered"] >= 7

    def test_dual_channel_combines_both_sources_correctly(self):
        scored = pd.DataFrame([
            {"cell_type": "R7", "matched_patterns": "", "entropy_channel_significant": True},
            {"cell_type": "R8", "matched_patterns": "histamine_blindspot", "entropy_channel_significant": True},
            {"cell_type": "R1-6", "matched_patterns": "histamine_blindspot", "entropy_channel_significant": False},
        ])
        report = recovery_report(scored)
        info = report["histamine_blindspot"]
        assert info["recovered_via"]["R7"] == "entropy"
        assert info["recovered_via"]["R8"] == "simplex+entropy"
        assert info["recovered_via"]["R1-6"] == "simplex"
        assert info["n_recovered"] == 3

    def test_dual_channel_excludes_seeds_caught_by_neither(self):
        scored = pd.DataFrame([
            {"cell_type": "R7", "matched_patterns": "", "entropy_channel_significant": False},
        ])
        report = recovery_report(scored)
        info = report["histamine_blindspot"]
        assert "R7" not in info["recovered_via"]
        assert info["missed"] == ["R7"]

    def test_entropy_reconstruction_matches_real_zscores(self):
        from entropy_channel import validate_reconstruction_against_real_zscores
        from paths import ENTROPY_CORRECTED

        if not ENTROPY_CORRECTED.exists():
            pytest.skip("results/entropy_corrected.csv not present in this environment")

        real = pd.read_csv(ENTROPY_CORRECTED)
        real["counts"] = real["nt_distribution"].map(parse_nt_distribution)
        check = validate_reconstruction_against_real_zscores(real, n_permutations=2000)

        correlation = check["z_score"].corr(check["entropy_z_reconstructed"])
        assert correlation > 0.99
        assert check["abs_diff"].median() < 0.5

    def test_signature_scan_matches_validated_real_data_numbers(self):
        from signature_scan import load_entropy_table
        from paths import ENTROPY_RAW

        if not ENTROPY_RAW.exists():
            pytest.skip("results/entropy_raw.csv not present in this environment")

        df = load_entropy_table()
        scored = score_types(df)
        report = recovery_report(scored)

        novel_fraction = scored["is_novel_candidate"].mean()
        assert 0.01 < novel_fraction < 0.15, (
            f"got {novel_fraction:.1%} -- should be a small, reviewable minority, "
            "between 1% and 15% of scored types"
        )
        assert report["ORN_SER_confusion"]["n_recovered"] >= 7
        assert report["Dm_GLUT_confusion"]["n_recovered"] >= 1
        assert "R1-6" in set(scored["cell_type"])


    def test_r1_6_not_formally_significant_but_notably_close(self):
        """Check relative histamine proximity without requiring formal significance."""
        scored = score_types(self._tiny_entropy_frame())
        r16 = scored[scored["cell_type"] == "R1-6"].iloc[0]
        clean = scored[scored["cell_type"] == "CleanGABA"].iloc[0]
        assert r16["p_histamine_blindspot"] < clean["p_histamine_blindspot"]

    def test_clean_ach_type_is_not_a_novel_candidate(self):
        scored = score_types(self._tiny_entropy_frame())
        row = scored[scored["cell_type"] == "CleanACH"].iloc[0]
        assert row["is_novel_candidate"] in (False, 0)

    def test_clean_gaba_type_is_not_a_novel_candidate(self):
        scored = score_types(self._tiny_entropy_frame())
        row = scored[scored["cell_type"] == "CleanGABA"].iloc[0]
        assert row["is_novel_candidate"] in (False, 0)

    def test_nearby_orn_is_flagged(self):
        """Use the larger reference pool required by the ORN seed count."""
        rows = [
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
            ("NearORN", 30, "{'SER': 28, 'ACH': 2, 'GABA': 0, 'GLUT': 0, 'DA': 0, 'OCT': 0}"),
        ]
        rows += self._large_background_rows(n=120, seed=0)
        scored = score_types(self._frame_from_rows(rows))
        row = scored[scored["cell_type"] == "NearORN"].iloc[0]
        assert row["p_ORN_SER_confusion"] < 0.05

    def test_flagged_fraction_is_bounded(self):
        scored = score_types(self._tiny_entropy_frame())
        assert scored["is_novel_candidate"].mean() < 0.5

    def test_background_types_are_rarely_flagged(self):
        scored = score_types(self._tiny_entropy_frame())
        bg_names = [r[0] for r in self._background_rows()]
        bg = scored[scored["cell_type"].isin(bg_names)]
        assert bg["is_novel_candidate"].mean() < 0.20
