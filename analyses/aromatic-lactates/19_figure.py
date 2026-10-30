"""Forest plot of per-cohort case-control effects (SD units) for the main replicated items. Output fig_forest_key.png"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mwlib import *
A = pd.read_csv(OUT + 'key_items_per_study.csv'); A = A[A.model == 'base']
K = pd.read_csv(OUT + 'key_items_sensitivity.csv', index_col=0)
items = ['MODULE aromatic lactates', 'phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'imidazole lactate', 'indolelactate', 'indolepropionate (IPA)',
         'MODULE collagen turnover', 'trans-4-hydroxyproline', 'methionine sulfoxide / methionine', 'PC(14:0/16:0)', 'plasmalogen PC(P-16:0/18:2)',
         'glucuronate', '5-oxoproline / glutamate', 'carnitine', 'kynurenine', 'hippurate', 'caffeine (lifestyle reference)']
fig, ax = plt.subplots(figsize=(8, 9), dpi=150); fig.patch.set_facecolor('#fcfcfb'); ax.set_facecolor('#fcfcfb')
for i, it in enumerate(items):
    y = len(items) - i
    g = A[A.item == it]
    ax.scatter(g.est, y + np.linspace(-0.15, 0.15, len(g)), s=14, color='#9a9994', zorder=2, label='single cohort' if i == 0 else None)
    e, p, k = K.loc[it, 'base_est'], K.loc[it, 'base_p'], int(K.loc[it, 'base_k'])
    se = abs(e) / max(1e-9, __import__('scipy').stats.norm.isf(p / 2))
    ax.plot([e - 1.96 * se, e + 1.96 * se], [y, y], color='#2a78d6', lw=2, zorder=3)
    ax.scatter([e], [y], marker='D', s=40, color='#2a78d6', zorder=4, label='pooled (random effects, 95% CI)' if i == 0 else None)
    ax.text(1.02, y, f'k={k}  {e:+.2f}  p={p:.1g}', transform=ax.get_yaxis_transform(), va='center', fontsize=7.5, color='#52514e')
ax.axvline(0, color='#52514e', lw=0.8)
ax.set_yticks(range(len(items), 0, -1)); ax.set_yticklabels(items, fontsize=8, color='#0b0b0b')
ax.set_xlabel('ME/CFS minus controls (SD units, adjusted where covariates exist)', fontsize=8, color='#52514e')
for s in ['top', 'right']: ax.spines[s].set_visible(False)
ax.grid(axis='x', color='#e6e5e0', lw=0.6); ax.set_axisbelow(True)
ax.legend(loc='lower left', fontsize=7, frameon=False)
ax.set_title('Replicated blood metabolite differences across ME/CFS cohorts', fontsize=10, loc='left', color='#0b0b0b')
plt.tight_layout(); plt.savefig(OUT + 'fig_forest_key.png', facecolor=fig.get_facecolor(), bbox_inches='tight')
