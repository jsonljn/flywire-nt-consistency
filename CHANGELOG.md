# Changelog

## signature_scan.py: fixed 80% over-flagging bug

### The bug

`signature_scan.py` scores every cell type's predicted-NT probability vector
against three literature-confirmed "confusion fingerprint" seed families
(`histamine_blindspot`: R7/R8/R1-6; `ORN_SER_confusion`: 10 ORN glomerulus
types; `Dm_GLUT_confusion`: Dm12/Dm19/Dm1), flagging a type as a "novel
candidate" if it falls inside a Jensen-Shannon-divergence threshold of any
seed. On the real project data (402 cell types), this flagged **322 (80%)**
— clearly wrong; the test suite already had two failing tests pointing at it
(`test_recovers_histamine_seeds`, `test_clean_type_is_not_a_novel_candidate`).

### Root cause

The threshold was pooled per *pattern*, taken from the **maximum** pairwise
leave-one-out distance across that pattern's seeds, times a fixed margin.
R7's predicted profile is a genuine outlier within its own family — its
nearest same-pattern neighbor (R8) is 0.182 bits away, and its nearest seed
in an *unrelated* pattern (Dm12, a GABA/GLUT interneuron with nothing to do
with photoreceptors) is only 0.102 bits away:

```
JS(R7, R8)    = 0.182   (same family)
JS(R7, R1-6)  = 0.390   (same family)
JS(R7, Dm12)  = 0.102   (different family, geometrically closer!)
```

That single outlier set `loo_max['histamine_blindspot'] = 0.182`, and with
the `* 1.35` margin the shared threshold became **0.245** — about a quarter
of the entire simplex's diameter. At that radius, "in the histamine
neighborhood" stops meaning "resembles R1-6" and starts meaning "is
ACH-dominant with a bit of GLUT," which is an extremely common shape:
ACH/GABA/GLUT together account for 95%+ of all 402 cell types. A 100%-pure,
zero-entropy synthetic ACH type (`CleanACH` in the test fixture) sits 0.096
bits from R1-6 — comfortably inside that inflated radius, despite having no
relationship to the histamine blind spot at all.

A second, independent issue: `hDeltaK` (a central-complex neuron, unrelated
to the olfactory system) has an FAFB-predicted vector of 100% serotonin —
**identical** to several literature-confirmed ORN seeds (also 100% SER).
JS divergence between identical points is exactly zero. No threshold, however
tight, can separate two literally-identical probability vectors using only
those vectors; this is a real, provable limit of the geometry, not a
mistuned constant, and it hasn't been fully eliminated (see "Known
limitations" below).

### The fix

`_per_seed_thresholds()` calibrates a **per-seed**, not per-pattern,
detection radius, as the smaller of two independent bounds:

1. That seed's own leave-one-out distance to its nearest same-pattern
   neighbor × 1.3 (how tight is this seed's own family, ignoring outliers
   elsewhere in the family).
2. The 2nd percentile of that seed's JS distance to the entire background
   population of other cell types (how rare is it, generally, to be this
   close to this specific seed by chance).

(1) alone under-corrects for common categories (R1-6 is ACH-dominant, and
ACH is the majority category, so "close to R1-6" isn't inherently rare).
(2) alone can occasionally be looser than a pattern's true intra-family
spacing, which (1) then tightens. Taking the tighter bound fixes both
directions.

`score_types()` also no longer forces a single global `best_pattern` via
unconstrained argmin before checking membership — each candidate is checked
against all three patterns' own thresholds independently
(`matched_patterns`), so a type can't be excluded from its true pattern's
neighborhood just because an unrelated pattern happens to be numerically
closer (this was silently costing R7 its own recovery check).

`recovery_report()` now cross-checks the project's existing entropy
z-score/FDR channel (the one `analysis.py` already computes and
`validate_results.py` already treats as the significance bar, q < 0.05). R1-6
and, on real data, Dm12/Dm1 have *negative* z-scores — confidently
consistent, not inconsistent — which is structurally invisible to that
channel and the entire reason this module exists; they can only be
recovered through simplex geometry. R7/R8 are the reverse: geometrically
distant from R1-6 (see above), but large, unambiguous entropy outliers
(z=7.80, z=3.55) already caught by the existing method. A seed counts as
recovered if *either* channel independently catches it, and
`recovered_via` reports which one(s) — nothing is recovered silently.

### Before / after (402 real cell types)

| | Before | After |
|---|---|---|
| Flagged "novel candidates" | 322 (80.1%) | 28 (7.0%) |
| `histamine_blindspot` seeds recovered | 2/3 (R7 dropped to a different pattern) | 3/3 (R1-6 via simplex; R7, R8 via entropy z/FDR) |
| `ORN_SER_confusion` seeds recovered | — | 9/9 present in data, via simplex |
| `Dm_GLUT_confusion` seeds recovered | — | 0/2 present in data — see limitations |
| Unit tests | 21/23 passing | 25/25 passing (2 new regression guards added) |
| `validate_results.py` | 22/23 checks passing | 24/24 checks passing |

### Known limitations (not fixed, disclosed instead of hidden)

- **Exact-tie false positives.** `hDeltaK` still appears in the novel-candidate
  list with `best_js = 0.000000` — it is geometrically identical to a
  literature-confirmed ORN seed by coincidence, and no distance-based method
  can separate two identical points using only those points. Any row with
  `best_js` at or near zero should be treated as "indistinguishable from a
  seed by this method, needs an independent check" (e.g. cell-type name,
  connectivity via `connectivity_comparison.py`), not as a confirmed finding.
- **Dm12/Dm1 are not independently recovered on real data** by either
  channel. They still get scored and reported (see `signature_scan.csv`),
  but the honest result is that this specific pattern's seeds are too few
  (2 present out of 3 named) and too close to the general GABA/GLUT
  background to clear either bar. This is reported as `NOT recovered by
  either channel` in the scan output rather than papered over.
- **This is one column of evidence, not a verdict.** Every "novel candidate"
  in `signature_scan_novel.csv` should be read as "worth checking against
  literature/connectivity," the same caveat the README already applies to
  the general/n>=10 scan's candidates.

### Files changed

`signature_scan.py` (threshold calibration, multi-pattern matching, dual-channel
recovery), `tests/test_core.py` (added `z_score` to the synthetic fixture,
updated `test_recovers_histamine_seeds` to check *which* channel recovered
each seed, added `test_clean_gaba_type_is_not_a_novel_candidate` and
`test_flagged_fraction_is_bounded` as regression guards), `validate_results.py`
(updated the photoreceptor-recovery check to the dual-channel definition,
added a flag-rate sanity check).

## Follow-up: exact test, literature validation, dual-channel entropy reconstruction

The fix above (per-seed thresholds, dual-channel recovery) was developed independently
of the work below and reached the same headline diagnosis (322/402, 80%) by a different
route. This entry documents what changed on top of it.

### What changed

1. **Threshold replaced with an exact permutation p-value** (`signature_calibration.py`).
   Rather than a percentile-of-background threshold, each candidate is tested against
   the empirical population of ~385-400 other real (non-seed) FAFB cell types: what
   fraction of them are at least as close to this candidate as the true seeds are, and
   how surprising is that if the seeds were replaced by k independently chosen ordinary
   types instead? `p = 1 - (1 - q_close) ** k`, exact closed form, with add-one-half
   continuity correction so small reference pools don't collapse to a hard 0.0 or 1.0.
   A parametric version (resample from the dataset's pooled NT rates, matching the
   entropy screen's own null philosophy) was tried first and rejected — documented in
   `signature_calibration.py`'s docstring, since it is a real pitfall worth recording:
   it under-flagged the true seeds (R7/R8/R1-6 all came back p≈1) and over-flagged 117
   unrelated types at the same time, because the dataset average is not a realistic
   stand-in for "an unremarkable fly cell type."

2. **Entropy-channel q-values reconstructed exactly, not approximated** (`entropy_channel.py`).
   The committed `entropy_corrected.csv` has `z_score` but not the underlying
   permutation p/q-value, and the raw per-neuron table isn't in this bundle. Rather than
   a normal-theory approximation on z, this reconstructs the same permutation null from
   the aggregated NT-count table alone: a group's category counts under the existing
   full-dataset label permutation are marginally an exact multivariate-hypergeometric
   draw from the dataset-wide pooled counts (a standard property of random partitions),
   so sampling directly from that distribution reproduces what re-running the original
   permutation on raw data would give. Validated directly against this project's own
   real z-scores: correlation 0.999, median absolute difference 0.13 (see
   `tests/test_core.py::test_entropy_reconstruction_matches_real_zscores`).

3. **Literature cross-check applied to the surviving candidates** (`build_signature_corrections.py`).
   `gt_data.csv` (the same ground-truth source the rest of this project already uses)
   was fetched and applied to the calibrated candidate list. Of 13 novel candidates:
   **2 are genuine, literature-confirmed corrections** — `hDeltaK` (FAFB: SER;
   literature: ACH, Wolff et al. 2024, EASI-FISH, confidence 4/5) and `TmY16` (FAFB:
   GABA; literature: GLUT, Nern et al. 2024, EASI-FISH, confidence 4/5) — neither an ORN
   nor Dm-prefixed type, extending both patterns beyond their original namesake
   families; 3 are literature-confirmed but already correctly predicted by FAFB
   (geometric near-misses, not corrections); 8 have no literature entry; 1
   (`Mi15`) is a literature-contradicted false lead, structurally the same kind of
   catch as the Lai exclusion elsewhere in this project. `hDeltaK`'s exact-zero JS
   distance to an ORN seed (raised as an unresolved caveat in the entry above) turns
   out to resolve cleanly once literature is consulted: it's a real hit, not an
   unresolvable tie, precisely because the literature check is a second, independent
   source of evidence rather than more geometry.

4. **`np.add.at` replaced with flattened-index `np.bincount`** in
   `analysis.py::stratified_permutation_null` — ~1.7x faster at this project's real
   scale (139k neurons / 402 types / 1000 permutations, benchmarked directly, not
   estimated). Confirmed bit-for-bit identical output before switching
   (`tests/test_core.py::test_bincount_matches_add_at`), so this could not have
   silently changed any entropy, z-score, or p-value already in this project's results.

5. **R1-6's recovery claim made honest rather than forced.** The exact test's result for
   R1-6 on real data is p≈0.066 for its own pattern — notably closer to the histamine
   seeds than the dataset median, but not below the same p<0.05 bar used everywhere
   else in this fix. Where the previous version's threshold construction made R1-6's
   recovery close to guaranteed by construction (R1-6 and R8 are each other's nearest
   same-pattern neighbor, so R8's own LOO-derived threshold is built from — and then
   re-tested against — R1-6 itself), this version reports the honest, borderline
   number instead. This does not change what R1-6 *is*: its classification rests on the
   independent MCNS + literature evidence documented elsewhere in this project, not on
   this exploratory geometric method in isolation, which is the right way around.

### Before / after (402 real cell types, this version)

| | Original heuristic | This fix |
|---|---|---|
| Flagged "novel candidates" | 322 (80%) | 13 (3%) |
| Literature-confirmed genuine corrections among novel candidates | 0 (never checked) | 2 (hDeltaK, TmY16) |
| `ORN_SER_confusion` seeds recovered | 9/9 (by construction) | 7-8/9 via simplex and/or entropy |
| `Dm_GLUT_confusion` seeds recovered | 3/3 (by construction) | 2/2 present in n>=20 table, via simplex |
| `stratified_permutation_null` runtime at real scale | 3.46s / 1000 permutations | 2.04s / 1000 permutations |
| Unit tests | 21/23 passing | 39/39 passing |

### Files changed on top of the entry above

`signature_calibration.py` (new — exact test, replaces the per-seed-threshold
calibration), `entropy_channel.py` (new — exact entropy-channel reconstruction),
`build_signature_corrections.py` (new — literature cross-check and correction list),
`plot_calibration_comparison.py` (new — before/after figure), `nt_simplex.py`
(added `batch_js_divergence`/`nearest_seed_js_batch`), `analysis.py` (bincount),
`signature_scan.py` (rewired onto the exact test; kept the per-seed-threshold
version's dual-channel recovery concept, reimplemented against the new test),
`tests/test_core.py` (rewritten: isolated tests for the exact-test primitive, a
properly-sized fixture for patterns that need it, a real-data integration test),
`validate_results.py`, `README.md`, `run_pipeline.py`.
