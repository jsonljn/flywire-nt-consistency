# Neurotransmitter Consistency Within Cell Types Across FlyWire Connectomes

FlyWire summer internship project. Tests whether neurons sharing a cell type label agree on their predicted neurotransmitter, as Dale's law would suggest — and traces the strongest outliers to their actual cause.

## Headline result

An unbiased screen of 402 FAFB cell types found only 2 statistically significant outliers by neurotransmitter inconsistency: **R7 and R8**, the fly's photoreceptors. Cross-dataset comparison against MCNS then revealed why: **FAFB's neurotransmitter classifier (Eckstein et al., Cell 2024) only predicts six neurotransmitters, and histamine is not one of them** — despite photoreceptors being canonically histaminergic. As a result, every histaminergic neuron in FAFB gets force-classified into one of the six other categories, and which category varies neuron to neuron, producing exactly the high-entropy signal this project was built to detect.

This was cross-checked against MCNS (whose classifier does include histamine): R7, R8, Lai, and R1-6 all sit in the extreme right tail of FAFB's entropy distribution (93rd-99th percentile), and FAFB predicts zero histamine anywhere in its 139,255-neuron dataset. A subsequent check against literature ground truth confirmed the histamine explanation for R7 and R8 (and R1-6 at the n>=20 threshold), but **overturned it for Lai**: MCNS's classifier says Lai is histaminergic, but literature (Davis et al. 2020) verifies it as glutamatergic. Lai is excluded from the final correction list -- see "Deliverable" section below.

![Histamine blindspot](figures/histamine_blindspot.png)


## Reproducing

**FAFB:**
1. Download from Codex (codex.flywire.ai, FAFB v783, Download Data page):
   - Neurotransmitter Type Predictions -> `data/neurons.csv`
   - Cell Types -> `data/cell_types.csv`
   - Classification / Hierarchical Annotations -> `data/classification.csv`
2. Run:
   ```
   python merge_data.py
   python analysis.py data/merged_annotations.csv
   ```

**MCNS:**
1. Download from Codex (switch dataset to MCNS v0.9, Download Data page):
   - Neuron Attributes -> `data_mcns/neurons.csv`
2. Run:
   ```
   python normalize_mcns.py
   python analysis.py data_mcns/merged_annotations.csv
   ```

**Cross-dataset check:**
```
python histamine_pattern_check.py
```
Requires both `data/merged_annotations.csv` and `data_mcns/merged_annotations.csv` to already exist, plus `results/entropy_raw.csv` from the FAFB run.

Requires: pandas, numpy, scipy, matplotlib.

## Method

For each cell type with at least 20 labeled members:
1. Compute Shannon entropy of the neurotransmitter prediction distribution
2. Build a stratified permutation null (1000 shuffles) that preserves both the overall neurotransmitter counts and each cell type's group size, to avoid flagging small groups as inconsistent just by chance
3. Compute a size-corrected z-score per cell type from this null

For the cross-dataset check:
1. From MCNS, find cell types the classifier labels >=90% consistently as histaminergic
2. Match those cell types to FAFB by name (exact match, or stripping MCNS subtype suffixes like `R7y` -> `R7`)
3. Compare FAFB's raw entropy for those matched types against the full FAFB entropy distribution

## Results in detail

### FAFB, single-dataset screen (402 cell types, >=20 members)

Only 2 significant outliers at z > 2:

| Cell type | n | entropy | z-score | dominant NT (wrong) |
|---|---|---|---|---|
| R7 | 474 | 1.62 | 7.80 | GLUT (44%) |
| R8 | 475 | 1.38 | 3.55 | ACH (58%) |

### Cross-dataset validation

MCNS subtypes R7y, R7p, R7d, R8y, R8p, R8d, and their "unclear" bins are all **100% histamine**, zero entropy -- matching well-established fly photoreceptor biology. FAFB predicts **zero histamine anywhere** in its 139,255-neuron dataset, consistent with the classifier's fixed six-category output (ACh, GABA, glutamate, dopamine, octopamine, serotonin; histamine explicitly listed as unsolved in Eckstein et al., Cell 2024).

Four FAFB cell types independently confirmed histaminergic via MCNS:

| Cell type | FAFB entropy | Percentile vs. all 402 FAFB types | FAFB dominant (wrong) NT |
|---|---|---|---|
| R7 | 1.625 | ~99th | GLUT (44%) |
| R8 | 1.376 | ~98th | ACH (58%) |
| Lai | 1.239 | ~97th | GABA (50%) |
| R1-6 | 0.853 | 93rd | ACH (82%) |

Mean entropy of the 4 confirmed-histaminergic types: 1.27. Mean entropy across all 402 FAFB cell types: 0.16.

## Interpretation

The original R7/R8 finding is not genuine biological heterogeneity (e.g. distinct functional subtypes disagreeing on transmitter) and not random annotation noise. It is a systematic, predictable consequence of a documented classifier limitation: neurons that are truly histaminergic cannot ever be correctly classified by FAFB's 6-category model, and instead get distributed across whichever of the 6 available categories best matches their synapse ultrastructure -- inconsistently across neurons of the same true type, which is exactly what produces high measured entropy.

This reframes the "correction" question. Rather than flagging individual neurons as mislabeled, the useful output is: **any FAFB cell type with unexpectedly high NT entropy is a candidate for being a real histaminergic (or other missing-category) population**, not necessarily a genuinely inconsistent one. This is a testable, generalizable flag the FlyWire annotation team could use to find other histamine-blind-spot cell types beyond the 4 confirmed here.

## Connectivity comparison: does the wrong label carry real structure?

Knowing FAFB's classifier can never predict histamine for R7/R8 raises a follow-up question: when it's forced to guess among its 6 wrong categories, is that guess random noise, or does it correlate with something real -- such as the neuron's true anatomical subtype? Built per-neuron output connectivity profiles (normalized synapse weight to the top 30 downstream partner cell types) and tested whether neurons sharing the same (wrong) NT label are more similar to each other than to neurons with a different label.

![Connectivity PCA](figures/connectivity_pca.png)

| Cell type | n | within-label sim | between-label sim | permutation p-value | classifier CV accuracy | majority baseline |
|---|---|---|---|---|---|---|
| R7 | 473 | 0.557 | 0.494 | < 0.0005 | 54.1% | 44.7% |
| R8 | 474 | 0.726 | 0.716 | 0.095 | 58.0% | 58.4% |

**R7**: the wrong label carries real signal. Same-label neurons are significantly more similar in connectivity than different-label neurons (2000-permutation test, p < 0.0005), and a cross-validated classifier predicts the label from connectivity alone better than the majority-class baseline. This suggests the classifier's confusion for R7 is not arbitrary -- it may be tracking real anatomical/wiring variation (consistent with documented R7 subtypes such as pale/yellow ommatidia, though this project does not have direct subtype labels in FAFB to confirm that specific correspondence).

**R8**: no such structure. Within- and between-label similarity are statistically indistinguishable (p = 0.095), and the classifier cannot beat baseline. R8's misprediction looks like noise, not signal.

Method: `connectivity_comparison.py` (profile construction, permutation test, classifier test), `connectivity_plots.py` (PCA visualization). Uses FAFB's Connections (Filtered) edge list, synapses filtered to >=5 per connection.

## Sensitivity check: lowering the cell-type size threshold to n>=10

The primary result above uses a minimum cell type size of 20 members, matching the original project spec. Re-running the full pipeline at a more permissive minimum of 10 members (702 cell types, vs. 402) is a robustness check that also nearly doubles coverage. Combined with the same MCNS-based validation strategy, generalized to catch any single-transmitter-confirmed FAFB type with elevated entropy (not just histamine-specific), this surfaced two additional systematic confusion patterns beyond the histamine blind spot -- both involving categories the classifier *can* predict, so these are genuine confusion errors, not categorical gaps.

**Pattern 1 -- Categorical blind spot (histamine).** R7, R8, Lai (and R1-6 at n>=20 already reported above). Confirmed structural: HIST is not a predictable output category.

**Pattern 2 -- ORN serotonin confusion.** Of 53 ORN (olfactory receptor neuron) glomerulus types matched to MCNS, all 53 are confirmed cholinergic by MCNS (consistent with the well-established fact that all Drosophila ORNs are cholinergic), and 43 are correctly and cleanly predicted ACH by FAFB with near-zero entropy. But 10 specific glomerulus types (ORN_V, ORN_VM3, ORN_VA2, ORN_DA3, ORN_DA4m, ORN_DA4l, ORN_DM2, ORN_DM3, ORN_DL4, ORN_DL3) are instead predicted predominantly serotonergic. Serotonin is a valid output category for the classifier, so this is not a structural gap -- it looks like genuine, glomerulus-specific classifier confusion.

**Pattern 3 -- Dm glutamate confusion.** Of 13 Dm-numbered (distal medulla) cell types matched to MCNS and confirmed glutamatergic, only 6 are correctly predicted GLUT by FAFB. The other 7 (Dm12, Dm16, Dm20, Dm19, Dm1, Dm6, Dm9) are predicted predominantly GABA or ACH instead -- roughly half of this cell-type family is systematically misclassified.

| Pattern | Flagged types | Denominator | True/confirmed NT | FAFB's wrong guess |
|---|---|---|---|---|
| Categorical blind spot | R7, R8, Lai | -- | HIST (unpredictable) | GLUT / ACH / GABA |
| ORN confusion | 10 | 53 ORN types checked | ACH | SER |
| Dm confusion | 7 | 13 GLUT-confirmed Dm types | GLUT | GABA / ACH |

Method: `general_scan_n10.py` (broad MCNS-based consistency check across all n>=10 FAFB types, catching any confirmed-consistent type with elevated FAFB entropy, not just histamine), `three_patterns_summary.py` (pattern extraction and categorization). Full results in `results/general_scan_n10_full.csv` and `results/three_confusion_patterns.csv`.

Caveat: MCNS's own predictions are also a classifier output, not ground truth from wet-lab literature -- treated here as a strong, independent cross-check rather than absolute truth, consistent with how R7/R8/Lai/R1-6 were validated above. Some flagged cell types have small FAFB sample sizes (n as low as 10-17), so these should be treated as candidates for further investigation, not settled conclusions.

## Deliverable: correction lists validated against literature ground truth

Following team guidance, every candidate above was cross-checked against the literature-curated ground truth database ([flyconnectome/drosophila_neurotransmitters](https://github.com/flyconnectome/drosophila_neurotransmitters), `gt_data.csv`), which links cell types to verified neurotransmitters from published sources with confidence scores. This is authoritative in a way MCNS's predictions are not, since MCNS's classifier can itself be wrong.

**Cell type name matching.** Different datasets name the same cell type differently (e.g. FAFB's `R1-6` vs. MCNS's `R1-R6`). `name_matching.py` handles this safely: exact case/whitespace-insensitive matches, plus an explicit range-notation normalizer for the `X1-6` / `X1-X6` pattern. Anything that isn't a safe match is left unmatched rather than guessed at, since a wrong automated match could produce a wrong correction. R7/R8's split into MCNS subtypes (`R7y`, `R7p`, `R7d`, etc.) is handled as an explicit, documented aggregation rather than a generic prefix rule, to avoid incorrectly conflating e.g. `Dm1` with `Dm12`/`Dm19`.

**Result of literature cross-check** (`validate_against_literature.py`, `results/literature_validated_candidates.csv`):

| Outcome | Count | Cell types |
|---|---|---|
| Confirmed by literature | 15 | R7, R8, 10 ORN types, Dm12, Dm19, Dm1 |
| Contradicted by literature -- excluded | 1 | Lai (MCNS said histamine; literature/Davis et al. 2020 verifies glutamate, confidence 4/5) |
| No literature match -- unconfirmed, not corrected | 4 | Dm16, Dm20, Dm6, Dm9 |

The Lai exclusion is a real finding in its own right: MCNS's own classifier appears to be wrong for this specific type, which is exactly the kind of error this literature cross-check is designed to catch before it turns into a bad correction.

R8 is confirmed as a genuine co-transmitter (ACH **and** histamine both verified present, matching Xiao et al. 2023's finding independently) -- so R8 neurons already predicted ACH are not flagged as wrong, only the ones predicted something else.

**Final correction lists** (`corrections/corrections_fafb.csv`, `corrections/corrections_mcns.csv`): one row per neuron whose current prediction doesn't match the literature-verified transmitter(s), with the source citation, confidence score, and proposed action attached.

| Dataset | Neurons flagged | Cell types covered |
|---|---|---|
| FAFB | 1,118 | 15 |
| MCNS | 4 | 1 (ExR7, a related subtype not part of the original 15) |

`corrections/excluded_unconfirmed_candidates.csv` lists the 5 excluded/unconfirmed types for transparency, so nothing is silently dropped.

## Status / next steps

- [x] FAFB entropy analysis with size-corrected null
- [x] MCNS cross-dataset comparison
- [x] Traced R7/R8 outlier signal to a documented FAFB classifier limitation (no histamine category)
- [x] Validated the explanation against 4 independently confirmed histaminergic cell types
- [x] Scanned all 402 FAFB cell types against MCNS for additional histamine-blind-spot candidates (`full_histamine_scan.py`, `results/full_cross_dataset_scan.csv`) -- confirms R7, R8, Lai, R1-6 are the complete set found by this method; 246/402 FAFB types had any MCNS name match; one additional candidate (T1) was excluded for insufficient MCNS sample size (n=1)
- [x] Connectivity comparison for R7 and R8: tested whether the classifier's wrong-category guess correlates with real connectivity structure. R7 shows significant structure (p < 0.0005); R8 does not (p = 0.095)
- [x] Sensitivity check at n>=10 (702 cell types): surfaced two additional systematic confusion patterns (ORN serotonin confusion, Dm glutamate confusion) beyond the histamine blind spot
- [x] Cross-referenced all flagged candidates against the literature-curated ground-truth NT database (`flyconnectome/drosophila_neurotransmitters`) -- confirmed 15, excluded 1 (Lai, contradicted), left 4 unconfirmed
- [x] Built robust cell-type name matching handling FAFB/MCNS naming differences (e.g. R1-6 vs R1-R6), per team guidance
- [x] Produced final deliverable: per-dataset, per-neuron correction lists with literature citations and confidence scores (`corrections/corrections_fafb.csv`, `corrections/corrections_mcns.csv`)
- [ ] Add MANC as a third cross-dataset check where applicable (note: MANC is nerve-cord only, so it won't contain R7/R8/Lai/R1-6, but may still be useful for other cell types)
- [ ] Investigate the 4 unconfirmed Dm types (Dm16, Dm20, Dm6, Dm9) further -- no literature match found yet, may need a different ground-truth source or direct follow-up
- [ ] Write final report and presentation

## Data source

FlyWire Codex (codex.flywire.ai), FAFB v783 and MCNS v0.9. See [FlyWire citation guidelines](https://codex.flywire.ai) for attribution requirements.

## Key reference

Eckstein, N. et al. Neurotransmitter classification from electron microscopy images at synaptic sites in Drosophila melanogaster. *Cell* (2024). doi:10.1016/j.cell.2024.03.016 -- explicitly lists histamine as an unsolved prediction target, consistent with this project's finding.
