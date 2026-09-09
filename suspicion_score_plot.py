"""Visualize the suspicion score distribution and ranking.

Usage:
    python suspicion_score_plot.py
Requires corrections/corrections_fafb_scored.csv (run suspicion_score.py first).
Writes figures/suspicion_score.png.
"""
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('corrections/corrections_fafb_scored.csv')

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax = axes[0]
colors = {
    'ORN_SER_confusion': '#1f77b4', 'Dm_GLUT_confusion': '#ff7f0e',
    'categorical_blindspot_HIST': '#d62728',
    'Dm9_GLUT_confirmed_via_direct_literature_search': '#2ca02c',
    'hDeltaK_signature_scan_novel_candidate': '#9467bd',
    'TmY16_signature_scan_novel_candidate': '#8c564b',
}
for p in df['pattern'].unique():
    subset = df[df['pattern'] == p]
    ax.hist(subset['suspicion_score'], bins=25, alpha=0.5, label=p, color=colors.get(p, 'gray'))
ax.set_xlabel('Suspicion score')
ax.set_ylabel('Number of neurons')
ax.set_title('Suspicion score by pattern\n(categorical blind-spot uses E1 only -- not directly\ncomparable in scale to the full E1xE2 cases)')
ax.legend(fontsize=7)

ax2 = axes[1]
full_cases = df[df['score_type'] == 'full (E1 x E2)']
sc = ax2.scatter(full_cases['E1_confident_wrongness'], full_cases['E2_evidence_for_correct'],
                 c=full_cases['suspicion_score'], cmap='Reds', alpha=0.6, s=20)
plt.colorbar(sc, ax=ax2, label='Suspicion score')
ax2.set_xlabel('E1: classifier confidence in wrong answer')
ax2.set_ylabel('E2: probability mass on correct answer')
ax2.set_title('Non-categorical confusion cases only\n(top-right = most suspicious: confident + wrong)')

plt.tight_layout()
plt.savefig('figures/suspicion_score.png', dpi=150)
plt.close()
print("Saved figures/suspicion_score.png")

print("\nMean/median by score_type:")
print(df.groupby('score_type')['suspicion_score'].agg(['mean', 'median', 'count']).round(3))
