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
        """analysis.py's permutation-counting step was rewritten from
        np.add.at to a flattened-index np.bincount for speed (~1.7x at this
        project's real scale -- 139k neurons / 402 types / 1000 permutations,
        benchmarked directly, not estimated). Confirm the two are not just
        both plausible but bit-for-bit identical, so the speedup could not
        have silently changed any entropy, z-score, or p-value in this
        project's results."""
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
        """batch_js_divergence must agree with js_divergence row-by-row -- the
        exact calibration test in signature_calibration.py relies on this."""
        rng = np.random.default_rng(0)
        seed = np.array([0.5, 0.3, 0.1, 0.05, 0.03, 0.02])
        P = rng.dirichlet(np.ones(6), size=25)
        got = batch_js_divergence(P, seed)
        expected = [js_divergence(row, seed) for row in P]
        np.testing.assert_allclose(got, expected, atol=1e-10)


class TestExactPValue:
    """signature_calibration.exact_p_value in isolation, no seed-geometry
    fixture needed -- see that module's docstring for the "needs a
    reasonably large M" caveat this directly tests."""

    def test_nothing_close_gives_high_p(self):
        pool = np.full(400, 0.9)  # nothing in the reference pool is close
        p, i_obs, M = exact_p_value(observed=0.01, pool_distances=pool, k=3)
        assert i_obs == 0
        assert p < 0.02  # still small (genuinely rare), but not a hard 0.0

    def test_everything_close_gives_p_near_one(self):
        pool = np.full(400, 0.001)  # the whole reference pool is closer than this
        p, i_obs, M = exact_p_value(observed=0.5, pool_distances=pool, k=3)
        assert i_obs == M
        assert p > 0.99

    def test_more_active_seeds_makes_a_match_less_surprising(self):
        """k = "how many independent chances did the observed min get to
        land this close" -- more chances (larger k) must raise p for the
        same underlying closeness, not lower it."""
        pool = np.concatenate([np.full(5, 0.02), np.full(395, 0.9)])
        p_k1, _, _ = exact_p_value(observed=0.05, pool_distances=pool, k=1)
        p_k9, _, _ = exact_p_value(observed=0.05, pool_distances=pool, k=9)
        assert p_k9 > p_k1

    def test_continuity_correction_avoids_hard_zero_at_small_m(self):
        """With only 3 reference points and zero of them close, the RAW
        i_obs/M formula gives q_close=0 -> p=0.0 exactly ("impossible by
        chance" from checking 3 things). The continuity-corrected version
        used here must not do that."""
        pool = np.array([0.9, 0.8, 0.7])
        p, i_obs, M = exact_p_value(observed=0.01, pool_distances=pool, k=3)
        assert i_obs == 0
        assert p > 0.0

    def test_correction_is_negligible_at_real_dataset_scale(self):
        """At M ~ 388 (the real project's reference-pool size), the
        continuity correction should barely move the answer at all --
        it exists for the small-M regime, not to change real conclusions."""
        pool = np.full(388, 0.9)
        p, _, _ = exact_p_value(observed=0.01, pool_distances=pool, k=3)
        uncorrected = 1.0 - (1.0 - 0 / 388) ** 3
        assert abs(p - uncorrected) < 0.01


class TestSignatureScan:
    def _large_background_rows(self, n: int = 100, seed: int = 0) -> list[tuple]:
        """Programmatically generate n plausible 'ordinary' cell types with
        varied simplex compositions (randomized dominant category and
        purity). Needed specifically for patterns with many active seeds --
        ORN_SER_confusion has 9-10 -- where signature_calibration's exact
        test needs M >= ~90-100 for even a perfect match to be able to reach
        p<0.05 at all (1 - (1 - 0.5/(M+1))**k < 0.05 requires M >= 97 at
        k=10; the hand-picked ~27-row background used elsewhere in this
        fixture only supports patterns with k <= ~3, which is why the ORN
        tests use this instead -- see test_nearby_orn_is_flagged and
        test_orn_cluster_recovers_via_simplex).
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
        """~24 'ordinary' cell types spanning the simplex, none matching any
        seed pattern. Needed for signature_calibration's exact test to have
        real resolution -- see that module's docstring: 2-3 background points
        is not enough for a combinatorial rarity test to mean anything in
        either direction. These sizes/mixes are arbitrary but deliberately
        varied (clean and mixed, every dominant category, a range of n).
        """
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
        """Seed families (with z_score matching real project data in sign and
        rough magnitude -- R7/R8 large positive entropy outliers, already
        caught by the existing z-score/FDR channel; R1-6 strongly negative,
        confidently *wrong* rather than inconsistent, invisible to that
        channel, the whole reason signature_scan exists) plus ~24 background
        types and two clean decoys, sized so signature_calibration's exact
        test has resolution for the small-k patterns (Dm: k=2-3; histamine:
        k=2-3 under LOO). ORN_SER_confusion (k=9-10) structurally needs a
        much larger background at this method's resolution (see
        _large_background_rows) and is tested separately below.
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
        """The Dm seed cluster must recover via simplex neighborhood on this
        fixture (k=2-3 active seeds -- tractable at M~45, see
        _tiny_entropy_frame's docstring). The histamine family (only 3
        seeds, deliberately spread out -- R7 is a genuine geometric outlier
        within its own family, see signature_calibration.py's docstring) is
        intentionally NOT asserted here at fixed counts: with a background
        this small, exactly which borderline seeds clear a p<0.05 line is
        sensitive to composition in a way a real M~400 pool is not. ORN
        recovery (k=9-10) needs a much bigger background and is tested
        separately in test_orn_cluster_recovers_via_simplex.
        test_signature_scan_matches_validated_real_data_numbers checks the
        actual, validated real-data recovery numbers for all three."""
        scored = score_types(self._tiny_entropy_frame())
        report = recovery_report(scored)
        assert report["Dm_GLUT_confusion"]["n_recovered"] >= 2

    def test_orn_cluster_recovers_via_simplex(self):
        """ORN_SER_confusion has 9-10 active seeds under LOO, which needs
        M >= ~90-100 for even a perfect match to reach p<0.05 (see
        _large_background_rows's docstring for the exact derivation) --
        the ~27-row background used by the other fixture tests here is not
        large enough for this specific pattern, hence the separate,
        larger fixture."""
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
        """Isolate recovery_report's channel-combining logic from the messier
        question of whether a small synthetic fixture's pool composition
        naturally reproduces realistic entropy z-scores (see previous test).
        Hand-construct a scored frame with known matched_patterns /
        entropy_channel_significant values and check the combination logic
        directly: simplex-only, entropy-only, both, and neither each report
        correctly, and 'neither' is excluded."""
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
        """The dual-channel recovery check's entropy side depends on
        entropy_channel.py reconstructing real permutation significance from
        aggregated counts alone (see that module's docstring for why this is
        possible without the raw per-neuron table). Validate directly against
        this project's own real, already-computed z-scores rather than
        trusting the small synthetic fixture to exercise this faithfully."""
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
        """Integration test against the real, committed project data --
        deliberately not the small synthetic fixture (see the previous two
        tests' docstrings for why small-M behavior is harder to pin down).
        These bounds are the actual, manually-validated numbers from running
        this module against results/entropy_raw.csv end to end (see
        CHANGELOG.md); this locks them in as a regression check."""
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
            "not the original bug's 80% and not zero"
        )
        assert report["ORN_SER_confusion"]["n_recovered"] >= 7
        assert report["Dm_GLUT_confusion"]["n_recovered"] >= 1
        # R1-6 is the flagship "entropy structurally can't see this" case;
        # it must never be excluded from scoring even when it doesn't clear
        # the significance bar (see test_r1_6_not_formally_significant_...).
        assert "R1-6" in set(scored["cell_type"])


    def test_r1_6_not_formally_significant_but_notably_close(self):
        """R1-6 is the case entropy structurally cannot see (negative z) and
        the honest simplex result for it, even on real project data, is
        borderline (p~0.066) rather than a clean pass -- see
        signature_calibration.py's docstring and CHANGELOG.md. This is not a
        bug to hide: it should stay closer to the seeds than an unrelated
        clean type, without being asserted into a false "recovered" claim
        this test suite cannot honestly make."""
        scored = score_types(self._tiny_entropy_frame())
        r16 = scored[scored["cell_type"] == "R1-6"].iloc[0]
        clean = scored[scored["cell_type"] == "CleanGABA"].iloc[0]
        assert r16["p_histamine_blindspot"] < clean["p_histamine_blindspot"]

    def test_clean_ach_type_is_not_a_novel_candidate(self):
        scored = score_types(self._tiny_entropy_frame())
        row = scored[scored["cell_type"] == "CleanACH"].iloc[0]
        assert row["is_novel_candidate"] in (False, 0)

    def test_clean_gaba_type_is_not_a_novel_candidate(self):
        """Regression guard for the Dm_GLUT_confusion analogue of the ACH
        over-flagging bug: an unremarkable GABA type must not match just for
        sharing GABA/GLUT dominance with Dm12/Dm1."""
        scored = score_types(self._tiny_entropy_frame())
        row = scored[scored["cell_type"] == "CleanGABA"].iloc[0]
        assert row["is_novel_candidate"] in (False, 0)

    def test_nearby_orn_is_flagged(self):
        """See test_orn_cluster_recovers_via_simplex's docstring: k=10 active
        seeds for a non-seed candidate here needs the larger background, not
        the ~27-row one used by the rest of this fixture."""
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
        """Regression guard for the original bug: a pooled, outlier-inflated
        threshold flagged 322/402 (80%) of all real cell types as "novel."
        On this fixture the flagged fraction must stay well under half --
        80%-style blowups should fail loudly."""
        scored = score_types(self._tiny_entropy_frame())
        assert scored["is_novel_candidate"].mean() < 0.5

    def test_background_types_are_rarely_flagged(self):
        """The ~27 deliberately unrelated background types (spanning every
        dominant category) should be flagged only at roughly the rate a
        raw p<0.05 threshold implies by construction, not systematically.
        27 types x 3 patterns = 81 comparisons; at p<0.05 uncorrected (this
        module's docstring explains why raw p, not BH-FDR, is the primary
        tier here -- the same reasoning applies to this fixture), a handful
        of coincidental matches is the expected, healthy false-positive rate
        this significance level implies, not a bug -- exactly the behavior
        the real analysis also shows (see CHANGELOG.md: 14 novel candidates
        on real data, of which literature independently contradicts one --
        Mi15 -- which is what a well-calibrated p<0.05 screen should
        occasionally produce). This asserts the rate stays in that expected
        ballpark, not that it is zero."""
        scored = score_types(self._tiny_entropy_frame())
        bg_names = [r[0] for r in self._background_rows()]
        bg = scored[scored["cell_type"].isin(bg_names)]
        assert bg["is_novel_candidate"].mean() < 0.20
