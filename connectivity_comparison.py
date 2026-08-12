"""
Connectivity comparison: does the WRONG neurotransmitter prediction that
FAFB assigns to histaminergic neurons (R7, R8) encode real structure --
e.g. correlating with true anatomical subtype (R7y/R7p/R7d) via connectivity
patterns -- or is it essentially noise with no connectivity signal?

Approach:
1. For each R7 (or R8) neuron with an NT prediction, build a normalized
   output connectivity profile: fraction of its output synapses going to
   each of the top-K most common downstream partner cell types.
2. Test whether neurons sharing the same (wrong) NT label have more similar
   connectivity profiles than neurons with different labels, using:
   - A permutation test on mean within-label vs between-label cosine similarity
   - A cross-validated classifier: can connectivity profile predict NT label
     better than chance?
3. Repeat for R8.

If there's real structure, that's evidence the FAFB classifier's mispredictions
aren't just noise -- they may be picking up on the same real variation that
distinguishes anatomical subtypes. If there's no structure, the mispredictions
carry no extra information and should just be treated as a blanket "this cell
type is subject to the histamine blind spot" flag.
"""
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

from paths import (
    CONNECTIVITY_SUMMARY,
    FAFB_CONNECTIONS,
    FAFB_MERGED,
    RESULTS,
    ensure_output_dirs,
)

ensure_output_dirs()

MIN_SYNAPSES = 5   # ignore very weak connections
TOP_K_PARTNERS = 30  # dimensionality of the connectivity profile
N_PERMUTATIONS = 2000
RANDOM_SEED = 42

print("Loading FAFB annotations and connections...")
annotations = pd.read_csv(FAFB_MERGED)
connections = pd.read_csv(FAFB_CONNECTIONS)
print(f"  Annotations: {len(annotations)} neurons")
print(f"  Connections: {len(connections)} edges")

# Filter weak connections
connections = connections[connections['syn_count'] >= MIN_SYNAPSES]
print(f"  Connections after >= {MIN_SYNAPSES} synapse filter: {len(connections)}")

# Map root_id -> cell type, for labeling downstream partners
id_to_type = dict(zip(annotations['root_id'], annotations['primary_type']))


def build_connectivity_profiles(cell_type_name, top_k=TOP_K_PARTNERS):
    """
    For each neuron of the given cell type with an NT prediction, build a
    normalized output connectivity profile vector over the top_k most common
    downstream partner cell types (computed within this population).
    """
    neurons = annotations[
        (annotations['primary_type'] == cell_type_name) & annotations['nt_type'].notna()
    ][['root_id', 'nt_type']].copy()

    neuron_ids = set(neurons['root_id'])
    print(f"\n{cell_type_name}: {len(neurons)} neurons with NT prediction")

    # Get all outgoing edges from these neurons
    out_edges = connections[connections['pre_root_id'].isin(neuron_ids)].copy()
    out_edges['post_type'] = out_edges['post_root_id'].map(id_to_type)
    out_edges = out_edges.dropna(subset=['post_type'])
    print(f"  Outgoing edges to typed partners: {len(out_edges)}")

    # Determine top-K partner types by total synapse weight
    partner_weight = out_edges.groupby('post_type')['syn_count'].sum().sort_values(ascending=False)
    top_partners = partner_weight.head(top_k).index.tolist()
    print(f"  Top {len(top_partners)} partner types cover "
          f"{partner_weight.head(top_k).sum() / partner_weight.sum():.1%} of output weight")

    # Build profile matrix: one row per neuron, one column per top partner type
    profiles = defaultdict(lambda: np.zeros(len(top_partners)))
    partner_idx = {p: i for i, p in enumerate(top_partners)}

    for _, row in out_edges.iterrows():
        if row['post_type'] in partner_idx:
            profiles[row['pre_root_id']][partner_idx[row['post_type']]] += row['syn_count']

    # Only keep neurons that actually have output edges to typed partners
    valid_ids = [nid for nid in neurons['root_id'] if nid in profiles and profiles[nid].sum() > 0]
    print(f"  Neurons with usable connectivity profile: {len(valid_ids)}")

    X = np.array([profiles[nid] / profiles[nid].sum() for nid in valid_ids])  # normalize per neuron
    labels = neurons.set_index('root_id').loc[valid_ids, 'nt_type'].values

    return X, labels, valid_ids, top_partners


def cosine_sim_matrix(X):
    norm = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-12)
    return norm @ norm.T


def within_vs_between_test(X, labels, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED):
    """
    Permutation test: is mean within-label cosine similarity higher than
    mean between-label cosine similarity, more than expected by chance?
    """
    sim = cosine_sim_matrix(X)
    n = len(labels)
    labels = np.array(labels)

    iu = np.triu_indices(n, k=1)
    same_label = labels[iu[0]] == labels[iu[1]]

    observed_within = sim[iu][same_label].mean()
    observed_between = sim[iu][~same_label].mean()
    observed_diff = observed_within - observed_between

    rng = np.random.default_rng(seed)
    null_diffs = []
    for _ in range(n_permutations):
        shuffled = rng.permutation(labels)
        same_label_perm = shuffled[iu[0]] == shuffled[iu[1]]
        null_within = sim[iu][same_label_perm].mean()
        null_between = sim[iu][~same_label_perm].mean()
        null_diffs.append(null_within - null_between)

    null_diffs = np.array(null_diffs)
    p_value = (null_diffs >= observed_diff).mean()

    return {
        'observed_within_sim': observed_within,
        'observed_between_sim': observed_between,
        'observed_diff': observed_diff,
        'null_mean_diff': null_diffs.mean(),
        'null_std_diff': null_diffs.std(),
        'p_value': p_value,
    }


def classifier_test(X, labels, seed=RANDOM_SEED):
    """
    Cross-validated random forest: can connectivity profile predict the
    (wrong) NT label better than chance?
    """
    le = LabelEncoder()
    y = le.fit_transform(labels)

    # Only keep classes with enough members for stratified CV
    vc = pd.Series(y).value_counts()
    valid_classes = vc[vc >= 5].index
    mask = np.isin(y, valid_classes)
    X_f, y_f = X[mask], y[mask]

    if len(np.unique(y_f)) < 2:
        return None

    clf = RandomForestClassifier(n_estimators=200, random_state=seed, class_weight='balanced')
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    scores = cross_val_score(clf, X_f, y_f, cv=cv, scoring='accuracy')

    baseline = vc[valid_classes].max() / vc[valid_classes].sum()  # majority-class baseline

    return {
        'n_neurons_used': len(y_f),
        'n_classes': len(np.unique(y_f)),
        'cv_accuracy_mean': scores.mean(),
        'cv_accuracy_std': scores.std(),
        'majority_class_baseline': baseline,
    }


# ─────────────────────────────────────────────
# Run for R7 and R8
# ─────────────────────────────────────────────

results_summary = {}

for cell_type in ['R7', 'R8']:
    print("\n" + "=" * 70)
    print(f"CELL TYPE: {cell_type}")
    print("=" * 70)

    X, labels, ids, partners = build_connectivity_profiles(cell_type)

    if len(X) < 20:
        print(f"  Too few neurons with usable profiles ({len(X)}), skipping.")
        continue

    print(f"\n  Label distribution among profiled neurons:")
    print(pd.Series(labels).value_counts())

    print("\n  Running permutation test (within- vs between-label similarity)...")
    perm_result = within_vs_between_test(X, labels)
    print(f"    Within-label mean cosine similarity:  {perm_result['observed_within_sim']:.4f}")
    print(f"    Between-label mean cosine similarity: {perm_result['observed_between_sim']:.4f}")
    print(f"    Observed difference: {perm_result['observed_diff']:.4f}")
    print(f"    Null distribution: mean={perm_result['null_mean_diff']:.4f}, std={perm_result['null_std_diff']:.4f}")
    print(f"    P-value: {perm_result['p_value']:.4f}")

    print("\n  Running classifier test (can connectivity predict NT label?)...")
    clf_result = classifier_test(X, labels)
    if clf_result:
        print(f"    Neurons used: {clf_result['n_neurons_used']}, classes: {clf_result['n_classes']}")
        print(f"    CV accuracy: {clf_result['cv_accuracy_mean']:.3f} +/- {clf_result['cv_accuracy_std']:.3f}")
        print(f"    Majority-class baseline: {clf_result['majority_class_baseline']:.3f}")
    else:
        print("    Not enough classes with sufficient samples for classifier test.")

    results_summary[cell_type] = {
        'n_profiled': len(X),
        'permutation_test': perm_result,
        'classifier_test': clf_result,
    }

    # Save profile data for this cell type
    profile_df = pd.DataFrame(X, columns=[f'partner_{p}' for p in partners])
    profile_df.insert(0, 'root_id', ids)
    profile_df.insert(1, 'nt_label', labels)
    profile_df.to_csv(RESULTS / f"connectivity_profiles_{cell_type}.csv", index=False)

# ─────────────────────────────────────────────
# Summary table + plot
# ─────────────────────────────────────────────

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

summary_rows = []
for ct, res in results_summary.items():
    row = {
        'cell_type': ct,
        'n_profiled': res['n_profiled'],
        'within_sim': res['permutation_test']['observed_within_sim'],
        'between_sim': res['permutation_test']['observed_between_sim'],
        'diff': res['permutation_test']['observed_diff'],
        'p_value': res['permutation_test']['p_value'],
    }
    if res['classifier_test']:
        row['cv_accuracy'] = res['classifier_test']['cv_accuracy_mean']
        row['baseline'] = res['classifier_test']['majority_class_baseline']
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))
summary_df.to_csv(CONNECTIVITY_SUMMARY, index=False)
print(f"\nSaved {CONNECTIVITY_SUMMARY}")
