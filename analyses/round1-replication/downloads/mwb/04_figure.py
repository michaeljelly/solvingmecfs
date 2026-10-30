"""Forest plot: per-study patient-vs-control effects (SD units) and random-effects pooled estimates for the main leads (blood)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = '/mnt/project-files/mecfs-research/downloads/mwb/'
E = pd.read_csv(OUT + 'per_study_effects.csv'); M = pd.read_csv(OUT + 'meta_blood.csv').set_index('feature')
E = E[E.tissue.isin(['plasma', 'serum'])]
F = ['IPA (indolepropionate)', 'indolelactate', 'kynurenine', 'kynurenine/tryptophan ratio', 'carnitine', 'COMP long-chain acylcarnitines',
     'COMP long-chain free fatty acids', '3-hydroxybutyrate', 'COMP primary bile acids', 'COMP secondary bile acids', 'VMA (vanillylmandelate)', 'caffeine']
short = {s: s.split('_')[0] for s in E.study.unique()}
fig, axes = plt.subplots(3, 4, figsize=(15, 11), sharex=True)
INK, MUTED, DISC, POOL = '#1f2328', '#8c959f', '#b35900', '#0b5cad'
for ax, f in zip(axes.ravel(), F):
    b = E[E.feature == f].reset_index(drop=True)
    y = np.arange(len(b))[::-1] + 2
    for yi, r in zip(y, b.itertuples()):
        c = DISC if r.study.startswith('BioMapAI') else MUTED
        ax.plot([r.est - 1.96 * r.se, r.est + 1.96 * r.se], [yi, yi], color=c, lw=2, solid_capstyle='round')
        ax.plot(r.est, yi, 'o', color=c, ms=8, mec='white', mew=1.5)
    m = M.loc[f]
    ax.plot([m.ext_lo, m.ext_hi], [0.6, 0.6], color=POOL, lw=2.5)
    ax.plot(m.ext_est, 0.6, 'D', color=POOL, ms=9, mec='white', mew=1.5)
    ax.set_yticks(list(y) + [0.6]); ax.set_yticklabels([short[s] for s in b.study] + ['Pooled external'], fontsize=8, color=INK)
    ax.axvline(0, color='#d0d7de', lw=1, zorder=0)
    ax.set_title(f"{f.replace('COMP ', '')}\nexternal {m.ext_est:+.2f} SD, p={m.ext_p:.2g}, q={m.ext_q:.2g}", fontsize=9.5, color=INK, loc='left')
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)
for ax in axes[-1]: ax.set_xlabel('patients minus controls (SD units, 95% CI)', fontsize=9, color=INK)
fig.suptitle('ME/CFS vs controls in blood: BioMapAI (orange, discovery) and independent cohorts (grey); blue = random-effects pool of independent cohorts only', fontsize=11, color=INK, x=0.01, ha='left')
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(OUT + 'fig_forest_leads.png', dpi=130)
print('saved')
