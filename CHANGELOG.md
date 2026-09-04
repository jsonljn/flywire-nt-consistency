# Changelog

## Verified fix not propagated to committed results; two more blockers found

Follow-up to the entry below. The ExR7/ExR8 fix is real and verified live
against a real Codex MCNS v1.0 download (`data_mcns/neurons.csv`,
`normalize_mcns.py` → `general_scan_n10.py`): no drop warning, R7, R8,
R1-6, and Lai all resolve correctly, matching the project's headline
numbers exactly (R1-6: n=4090, 81.9% ACH, MCNS 750/750 HIST).

That result was **not** committed into `results/general_scan_n10_full.csv`
/ `general_scan_n10_flagged.csv`. Reason: this repo is also missing
`results/entropy_raw_n10.csv` (the FAFB n>=10 sensitivity screen), which
needs real `data/merged_annotations.csv` to rebuild, also unavailable here.
Without it, `general_scan_n10.py` silently falls back to the coarser
n>=20 `entropy_raw.csv`. Checked directly: the resulting match list drops
from 381 types (original, committed) to 245 (this session's regeneration).
It gains R1-6, and loses 136 legitimate n=10-19 matches, including
`ORN_DA3` and `Dm19`, both already present in the committed correction
lists. That's a large net loss for one gain. Reverted
(`git checkout db1dd76 -- results/general_scan_n10_full.csv
results/general_scan_n10_flagged.csv`) rather than keep the smaller file.
`results/three_confusion_patterns.csv`, `results/full_cross_dataset_scan.csv`,
`results/signature_scan.csv`, and `results/signature_scan_novel.csv` were
also regenerated this session and also reverted for the same reason (the
signature_scan files additionally lose their `gt_verified_nt` /
`gt_agrees_with_pattern` columns without `data/gt_data.csv`, a separate,
also-unavailable-here file).

### Full blocker list, for whoever runs this next

1. `data_mcns/neurons.csv` -- now resolved (Codex, MCNS v1.0, "Neuron
   Attributes"). Note the version: this repo's own instructions and
   citations say MCNS v0.9 throughout; v1.0 is what was actually used to
   verify the fix above. Not checked: whether v0.9 had the same ExR7/ExR8
   contamination, or whether R1-6's real MCNS name changed between
   versions. The `EXPLICIT_SUBTYPE_GROUPS` fix itself doesn't depend on
   which version turns out to be true.
2. `data/merged_annotations.csv` (FAFB, from `neurons.csv` +
   `cell_types.csv` + `classification.csv` via `merge_data.py`) -- still
   needed, to rebuild `results/entropy_raw_n10.csv` via
   `analysis.py data/merged_annotations.csv --min-members 10`.
3. `data/gt_data.csv` -- still needed, for `validate_against_literature.py`
   and `build_corrections.py` to regenerate `corrections_fafb.csv` and the
   1,390/18 totals against the fix above. Public source:
   github.com/flyconnectome/drosophila_neurotransmitters (also archived on
   Zenodo, DOI 10.5281/zenodo.20818142). Not fetched here: no network
   access in this environment to pull external files into the repo.

Once both `data/merged_annotations.csv` and `data/gt_data.csv` are in
place: `python run_pipeline.py` should regenerate everything end to end,
correctly, in one pass.

## R1-6/R7/R8 resolved against real MCNS data; two earlier hypotheses were wrong

Real `data_mcns/neurons.csv` (Codex, MCNS v1.0, "Neuron Attributes") became
available and was run through the actual pipeline. Result: R1-6 was never
broken. R7 and R8 were, for an unrelated reason. Both prior guesses in this
file (below) were wrong; recorded here rather than deleted, since they were
reasonable given the evidence available at the time and the way they were
wrong is itself useful.

### What was actually wrong

`EXPLICIT_SUBTYPE_GROUPS["R7"]` and `["R8"]` (`mcns_matching.py`) included
`"ExR7"`/`"ExR8"`. These are Extrinsic Ring neurons of the ellipsoid body
(Hanesch et al. 1989), an unrelated cell class that happens to share a name
substring with the photoreceptors. Verified: ExR7 (n=4, real MCNS data) is
100% ACH, while true R7 photoreceptor subtypes (R7y/R7p/R7d/R7_unclear) are
100% HIST where predicted. Mixing them into one group meant
`all_subtypes_consistent` failed for every R7 and R8 lookup, so both were
silently dropped, every run, regardless of threshold values.

R1-6 needed no fix. It resolves via plain range notation (`"R1-6" ->
"R1-R6"`) exactly as originally written, and reproduces the project's
headline numbers exactly against real data: n=4,090, 81.9% predicted ACH
(matches the "82%" cited throughout), MCNS confirms 750/750 HIST.

### Why the two earlier hypotheses in this file were wrong

1. The original claim (below, "Known issue" entry) guessed the drop happened
   *after* a successful match, at the `mcns_frac < 0.9` / `mcns_n < 10`
   threshold check. Wrong: R1-6 never had a threshold problem, its real
   `mcns_frac` is 1.0 and `mcns_n` is 750, both well past threshold. The
   drop (of R7/R8, not R1-6) happened one step earlier, at the
   `all_subtypes_consistent` check, which that hypothesis didn't consider.

2. The revised claim (below, "Known issue: R1-6 missing" entry) guessed R1-6
   itself needed a subtype-group fallback because it has none, unlike
   R7/R8. Wrong in the opposite direction: R7/R8's *having* an (incorrect)
   subtype group is what broke them. R1-6 having no subtype group was never
   the problem; it doesn't need one. Both guesses were made without the raw
   MCNS file and turned out to be reasonable-sounding but incorrect; this is
   why the "not fixed here, do not guess" note in that entry existed, and
   why it was right to leave the actual `EXPLICIT_SUBTYPE_GROUPS` code
   untouched until real data was available to check against.

### What was fixed

Removed `"ExR7"`/`"ExR8"` from `EXPLICIT_SUBTYPE_GROUPS`. Re-ran
`general_scan_n10.py` against real data: no `WARNING` for R7, R8, R1-6, or
Lai; all four now appear in `results/general_scan_n10_full.csv` with
`is_categorical_blindspot=True`, matching the headline claims. Not yet done:
re-running the rest of the pipeline (`signature_scan.py`,
`validate_against_literature.py`, `build_corrections.py`) to regenerate
`corrections_fafb.csv` and the committed 1,390/18 totals against this fix,
and normalizing FAFB's own raw download the same way MCNS's was (FAFB's
`neurons.csv`/`cell_types.csv`/`classification.csv` use Title-Case column
names from Codex; `merge_data.py` still expects snake_case, the same class
of bug this entry just fixed on the MCNS side).

## Known issue: R1-6 missing from corrections_fafb.csv

R1-6 (4,090 neurons, the project's flagship finding) is present in
`results/signature_scan.csv`, correctly flagged as `histamine_blindspot`
with strong evidence (z=-16.8, JS-distance to the histamine seed family well
below the calibrated threshold), but is **absent** from
`results/general_scan_n10_full.csv`, `three_confusion_patterns.csv`,
`literature_validated_candidates.csv`, and ultimately `corrections_fafb.csv`.
Verified directly: zero matches for "R1-6" in all four files. It is silently
dropped somewhere inside `general_scan_n10.py`'s per-type loop.

### What was checked

`resolve_fafb_to_mcns("R1-6", ...)` only reaches a match via range notation
(`"R1-6" -> "R1-R6"`), which requires MCNS's real `primary_type` column to
contain that exact combined string. This was previously "confirmed correct"
by `tests/test_core.py::test_range_notation` — but that test uses a synthetic
name list that already assumes `"R1-R6"` is the real spelling; it does not,
and cannot, prove that against real data, since the raw MCNS file
(`data/data_mcns/merged_annotations.csv`) is gitignored and not in this
bundle.

R7 and R8 do not depend on this fragile a path: both have an explicit entry
in `EXPLICIT_SUBTYPE_GROUPS` (`mcns_matching.py`) because MCNS is known to
split them into named subtypes (`R7y`, `R7p`, `R7d`, ...). R1-6 has no such
entry. If MCNS also splits R1-6's subtypes into individual names rather than
one combined `"R1-R6"` string, the same way it does for R7/R8, range
notation cannot find it, and there is no subtype-group fallback to catch
that case for R1-6 the way there is for R7/R8. Verified directly: a
synthetic MCNS name list using split subtypes (`R1`, `R2`, ..., `R6`)
produces the same silent `(None, None)` this project's real data does;
only the exact single combined string `"R1-R6"` succeeds
(`tests/test_core.py::test_r1_6_has_no_subtype_fallback_unlike_r7_r8`).

A second, independent lead was ruled out: `results/full_cross_dataset_scan.csv`
records R1-6's match method as `"variant (R1-6 -> R1-R6)"` — a label the
current `resolve_fafb_to_mcns` can never produce (it only ever returns
`explicit_subtype_group`, `exact`, `exact_case_insensitive`,
`range_notation`, or `subtype_group`). That CSV must have been generated by
an earlier, more permissive matcher that no longer exists in this codebase.
It is stale and should not be read as evidence that the current matcher
succeeds against real data.

### What was fixed, and what was not

Fixed: `general_scan_n10.py` and `full_histamine_scan.py` now track a
`WATCH_TYPES` set (`R7`, `R8`, `R1-6`, `Lai`) and print a loud `WARNING`
with the specific drop reason (no match / empty lookup / threshold miss) if
any of them fails to reach `results`. This class of failure can no longer
happen silently. Verified against a synthetic MCNS dataframe reproducing
the suspected failure mode (no-match case): the warning fires with the
correct reason string.

Not fixed: the actual root cause requires the real
`data/data_mcns/merged_annotations.csv` to confirm the exact MCNS spelling
for R1-6's subtypes, which is not available in this environment. Do not
guess at the fix (e.g. hardcoding a subtype list into
`EXPLICIT_SUBTYPE_GROUPS`) without first checking the real
`primary_type` values — a wrong guess would silently produce a wrong
correction, the exact failure mode `name_matching.py` is designed to avoid.
Run `general_scan_n10.py` against the real data and read its `WARNING`
output first; it will now state the exact reason R1-6 (or any watched
type) is dropped.

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
