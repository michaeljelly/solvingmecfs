"""Forest plot of replication (standardised effects, 95% CI) and symptom between/within plot."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from data import OUT

COL = {'BioMapAI (all visits)': '#2a78d6', 'Raijmakers2020_stool': '#eb6834', 'Guo2023_Columbia_stool': '#1baf7a',
       'Che2022_Columbia_plasma': '#eda100', 'Germain2022_plasma': '#e87ba4', 'External meta (FE)': '#333333'}
MK = {'BioMapAI (all visits)': 'o', 'Raijmakers2020_stool': 's', 'Guo2023_Columbia_stool': 'D', 'Che2022_Columbia_plasma': '^',
      'Germain2022_plasma': 'v', 'External meta (FE)': 'P'}
E = pd.read_csv(OUT + 'replication_per_cohort.csv')
R = pd.read_csv(OUT + 'casecontrol_all_models.csv.gz'); R = R[(R.model == 'LMM_all') & (R.contrast == 'patient_vs_control')].set_index('feature')
T = pd.read_csv(OUT + 'replication_table.csv').set_index('feature')
feats = ['Faecalibacterium_prausnitzii', 'Roseburia_faecis', 'Roseburia_inulinivorans', 'Coprococcus_eutactus', 'Coprococcus_comes', 'Anaerostipes_hadrus',
         'COMPOSITE butyrate producers (log10 rel. abund.)', 'Ruminococcus_gnavus', 'Flavonifractor_plautii', 'Dysosmobacter_welbionis', 'Enterocloster_bolteae',
         'Clostridium_symbiosum', 'kynurenine', 'vanillylmandelate (VMA)', 'choline', 'paraxanthine', 'quinate', 'COMPOSITE primary bile acids',
         'COMPOSITE secondary bile acids', 'COMPOSITE ceramides', 'COMPOSITE sphingomyelins', 'COMPOSITE long-chain free fatty acids', 'COMPOSITE long-chain acylcarnitines']
# BioMapAI composites are not in R -> take from replication table estimate with se unknown; recompute se from p
fig, ax = plt.subplots(figsize=(8, 11))
y = 0; yt = []; yl = []
for f in feats[::-1]:
    rows = []
    if f in R.index: rows.append(('BioMapAI (all visits)', R.loc[f, 'est'], R.loc[f, 'se']))
    elif f in T.index and not np.isnan(T.loc[f, 'BioMapAI_est']):
        from scipy import stats
        z = abs(stats.norm.isf(T.loc[f, 'BioMapAI_p'] / 2)); rows.append(('BioMapAI (all visits)', T.loc[f, 'BioMapAI_est'], abs(T.loc[f, 'BioMapAI_est']) / max(z, 1e-6)))
    for _, x in E[E.feature == f].iterrows(): rows.append((x.cohort, x.est, x.se))
    if f in T.index and 'ext_meta_FE' in T and not np.isnan(T.loc[f, 'ext_meta_FE']):
        rows.append(('External meta (FE)', T.loc[f, 'ext_meta_FE'], (T.loc[f, 'ext_meta_FE_hi'] - T.loc[f, 'ext_meta_FE']) / 1.96))
    n = len(rows)
    for i, (c, e, s) in enumerate(rows):
        yy = y + (i - (n - 1) / 2) * 0.13
        ax.errorbar(e, yy, xerr=1.96 * s, fmt=MK[c], color=COL[c], ms=5, lw=1.2, capsize=0)
    yt.append(y); yl.append(f.replace('COMPOSITE ', '').replace('_', ' ').replace(' (log10 rel. abund.)', ''))
    y += 1
ax.axvline(0, color='#999999', lw=1)
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=8)
ax.set_xlabel('Patient vs control, standardised difference (SD units, 95% CI)')
ax.set_title('Replication of BioMapAI case-control signals in independent cohorts', fontsize=10)
for k in COL: ax.plot([], [], MK[k], color=COL[k], label=k)
ax.legend(fontsize=7, loc='lower right', frameon=False)
ax.grid(axis='x', color='#eeeeee'); ax.set_axisbelow(True)
for s in ['top', 'right']: ax.spines[s].set_visible(False)
plt.tight_layout(); plt.savefig(OUT + 'fig_replication_forest.png', dpi=130)

# symptom composites: between vs within, standardised by patient SD of the symptom
C = pd.read_csv(OUT + 'symptom_targeted_composites.csv')
sd = {'orthostatic_NLT': 0.117, 'sleep_PSQI': 0.189, 'fatigue_energy': 0.119, 'cognition_DANA': 0.156, 'overall_severity': 0.107}
sel = [('FFA_long_chain', 'orthostatic_NLT'), ('acylcarnitines_long_chain', 'orthostatic_NLT'), ('hydroxy_FA', 'orthostatic_NLT'), ('BHB', 'orthostatic_NLT'),
       ('kynurenine', 'overall_severity'), ('kyn_trp_ratio_log2', 'overall_severity'), ('secondary_bile_acids', 'cognition_DANA'), ('butyrate_producers_log10', 'overall_severity')]
fig, ax = plt.subplots(figsize=(7, 4.5))
for i, (f, s) in enumerate(sel[::-1]):
    r = C[(C.feature == f) & (C.symptom == s)].iloc[0]
    for j, (comp, col, mk) in enumerate([('between', '#2a78d6', 'o'), ('within', '#eb6834', 's')]):
        e = r[comp + '_est'] / sd[s]; se = r[comp + '_se'] / sd[s]
        ax.errorbar(e, i + (0.12 if j == 0 else -0.12), xerr=1.96 * se, fmt=mk, color=col, ms=5, lw=1.2)
ax.set_yticks(range(len(sel))); ax.set_yticklabels([f'{f} -> {s}' for f, s in sel[::-1]], fontsize=8)
ax.axvline(0, color='#999999', lw=1)
ax.plot([], [], 'o', color='#2a78d6', label='between-person (person mean)'); ax.plot([], [], 's', color='#eb6834', label='within-person (visit deviation)')
ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.set_xlabel('Symptom SDs per 1 SD of metabolite (+ = worse), 95% CI')
ax.set_title('Patients, all visits: associations are between-person only', fontsize=10)
for s_ in ['top', 'right']: ax.spines[s_].set_visible(False)
plt.tight_layout(); plt.savefig(OUT + 'fig_symptom_between_within.png', dpi=130)
print('ok')
