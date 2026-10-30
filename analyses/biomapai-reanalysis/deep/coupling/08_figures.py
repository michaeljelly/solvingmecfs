"""Figures: (1) gut-species predictability of each metabolite, controls vs patients (matched n);
(2) case-control effect sizes for microbe-derived metabolites in this cohort and two external cohorts."""
from cc import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e', 'ytick.color': '#52514e'})
BLUE, ORANGE, GRAY = '#2a78d6', '#eb6834', '#b5b4ad'

d = pd.read_csv(OUT + 'pred_metab_species_all.csv', index_col=0)
fig, ax = plt.subplots(figsize=(5.2, 5))
pr = (d.ctrl_q < .05) | (d.pat_q < .05)
ax.scatter(d.ctrl_R2[~pr], d.pat_R2[~pr], s=9, color=GRAY, alpha=.6, lw=0, label='not gut-predictable (821 metabolites total)')
ax.scatter(d.ctrl_R2[pr], d.pat_R2[pr], s=36, color=BLUE, edgecolor='#fcfcfb', lw=1.5, label='gut-predictable in either group (q<0.05)')
for n in d[pr].index:
    if d.loc[n, ['ctrl_R2', 'pat_R2']].max() > .09:
        ax.annotate(n.replace(' (hydrocinnamate)', ''), (d.ctrl_R2[n], d.pat_R2[n]), fontsize=7.5, color='#0b0b0b',
                    xytext=(4, -3), textcoords='offset points')
lim = (-0.12, 0.32); ax.plot(lim, lim, color='#52514e', lw=.8, ls='--'); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('Controls: out-of-sample R² from gut species'); ax.set_ylabel('Patients (subsampled to same n): R²')
ax.set_title('Gut species predict few plasma metabolites,\nand about equally in patients and controls', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=7.5, loc='upper left')
fig.tight_layout(); fig.savefig(OUT + 'fig1_predictability_species_metabolome.png', dpi=160)

L = pd.read_csv(OUT + 'levels_mixed.csv').set_index('metabolite')
E = pd.read_csv(OUT + 'external_replication.csv')
mets = ['indolepropionate', 'indolelactate', 'hippurate', 'catechol sulfate', 'cinnamoylglycine', 'p-cresol glucuronide*',
        'glycodeoxycholate 3-sulfate', 'taurodeoxycholate', 'deoxycholate', 'kynurenate', '3-indoxyl sulfate']
fig, ax = plt.subplots(figsize=(6.2, 5.4))
y = 0; ticks = []; labels = []
for m in mets[::-1]:
    ests = [('This cohort (BioMapAI)', L.loc[m, 'beta_SD'], L.loc[m, 'lo'], L.loc[m, 'hi'], BLUE)]
    g = E[(E.cohort == 'Germain2022_D1PRE') & (E.metabolite == m)]
    if len(g): ests.append(('Germain 2022', g.beta_SD.iloc[0], g.lo.iloc[0], g.hi.iloc[0], ORANGE))
    if m == 'indolepropionate':
        c = E[E.metabolite == 'indolepropionate [Indole-3-propionic acid]']
        ests.append(('Che 2022', c.beta_SD.iloc[0], c.lo.iloc[0], c.hi.iloc[0], '#1baf7a'))
    for j, (nm, b, lo, hi, col) in enumerate(ests):
        yy = y + j * 0.28
        ax.plot([lo, hi], [yy, yy], color=col, lw=2, solid_capstyle='round')
        ax.scatter([b], [yy], color=col, s=30, zorder=3, edgecolor='#fcfcfb', lw=1.2, label=nm if m == 'indolepropionate' else None)
    ticks.append(y + 0.14 * (len(ests) - 1)); labels.append(m); y += 1.2
ax.axvline(0, color='#52514e', lw=.8)
ax.set_yticks(ticks); ax.set_yticklabels(labels)
ax.set_xlabel('Patients minus controls (SD units, adjusted; 95% CI)')
ax.set_title('Indolepropionate is lower in ME/CFS in all 3 cohorts;\nmost other shifts do not replicate', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=7.5, loc='upper center', bbox_to_anchor=(0.45, -0.1), ncol=3)
fig.tight_layout(); fig.savefig(OUT + 'fig2_microbial_metabolites_levels_3cohorts.png', dpi=160)
print('ok')
