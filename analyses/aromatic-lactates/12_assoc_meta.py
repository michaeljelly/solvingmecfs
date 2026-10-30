"""Metabolome-wide case-control effects per cohort, random-effects meta-analysis, and discovery/replication designs.
Outputs: per_study_all.csv, meta_all.csv, replication_designs.csv
"""
import numpy as np, pandas as pd
from mwlib import *

S, M = load()
rows = []
for st, d in S.items():
    L, C = d['L'], d['C']
    r = ols_effect(L, C, 'case', covs(st, C))
    r['study'] = st; r['feature'] = r.index; r['n_case'] = int(C.case.sum()); r['n_ctrl'] = int((1 - C.case).sum())
    rows.append(r)
# Hoel summary
h = pd.read_csv(W + 'hoel_summary.csv')
dd, se = d_from_fc_p(h.FC.astype(float).values, h.p.astype(float).values)
rows.append(pd.DataFrame({'est': dd, 'se': se, 'p': h.p.values, 'n': sum(HOEL_N), 'study': 'Hoel2021', 'feature': h.Met.values, 'n_case': HOEL_N[0], 'n_ctrl': HOEL_N[1]}))
E = pd.concat(rows, ignore_index=True)
E = E.merge(M[['hkey', 'label', 'study', 'feature', 'super', 'sub']], on=['study', 'feature'], how='inner')
E['family'] = E.study.map(FAMILY)
E['q_within'] = E.groupby('study').p.transform(bh)
E.to_csv(OUT + 'per_study_all.csv', index=False)
print(E.groupby('study').agg(n_feat=('hkey', 'size'), n_q05=('q_within', lambda q: (q < 0.05).sum())))

lab = M.drop_duplicates('hkey').set_index('hkey')[['label', 'super', 'sub']]


def run_meta(E, studies, tag):
    e = E[E.study.isin(studies)]
    out = []
    for hk, g in e.groupby('hkey'):
        r = meta(g.est, g.se); r['hkey'] = hk; r['studies'] = ','.join(g.study)
        r['n_pos'] = int((g.est > 0).sum()); r['n_neg'] = int((g.est < 0).sum()); r['N'] = int((g.n_case + g.n_ctrl).sum())
        out.append(r)
    R = pd.DataFrame(out).set_index('hkey'); R.columns = [f'{tag}_{c}' for c in R.columns]
    return R


ALL = list(S) + ['Hoel2021']
MET = [s for s in ALL if FAMILY[s] == 'Metabolon']
OTH = [s for s in ALL if FAMILY[s] == 'Other']
NOHANSON = [s for s in ALL if s not in ('Germain2018', 'Germain2020')]
R = lab.join(run_meta(E, ALL, 'all'), how='inner').join(run_meta(E, MET, 'met')).join(run_meta(E, OTH, 'oth')).join(run_meta(E, NOHANSON, 'noH')) \
       .join(run_meta(E, [s for s in ALL if s != 'BioMapAI'], 'exBM'))
for t in ['all', 'met', 'oth', 'noH', 'exBM']:
    m = R[f'{t}_k'] >= 2
    R[f'{t}_q'] = np.nan; R.loc[m, f'{t}_q'] = bh(R.loc[m, f'{t}_p'])
# BioMapAI alone
bmE = E[E.study == 'BioMapAI'].set_index('hkey')[['est', 'p']].rename(columns={'est': 'BM_est', 'p': 'BM_p'})
R = R.join(bmE); R['BM_k'] = R.BM_p.notna().astype(int)
R['BM_q'] = np.nan; m = R.BM_p.notna(); R.loc[m, 'BM_q'] = bh(R.loc[m, 'BM_p'])


def replicate(R, disc, rep, tag):
    """disc: prefix of discovery meta (need k>=2 unless BM); rep: prefix of replication meta (k>=1)."""
    kd = R[f'{disc}_k'] >= 2 if disc != 'BM' else R.BM_p.notna()
    D = R[kd & (R[f'{disc}_q'] < 0.05) & (R[f'{rep}_k'] >= 1)].copy()
    same = np.sign(D[f'{disc}_est']) == np.sign(D[f'{rep}_est'])
    D['rep_p1'] = np.where(same, D[f'{rep}_p'] / 2, 1 - D[f'{rep}_p'] / 2)
    D['rep_q'] = bh(D.rep_p1.values)
    D['replicated'] = (D.rep_q < 0.05)
    D['design'] = tag
    return D


cols_keep = ['label', 'super', 'sub']
des = []
for disc, rep, tag in [('met', 'oth', 'A: discover Metabolon -> replicate non-Metabolon'), ('oth', 'met', 'B: discover non-Metabolon -> replicate Metabolon'),
                       ('BM', 'exBM', 'C: discover BioMapAI -> replicate all other cohorts')]:
    D = replicate(R, disc, rep, tag)
    D = D[cols_keep + ['design', f'{disc}_est', f'{disc}_q', f'{disc}_k', f'{rep}_est', f'{rep}_se', f'{rep}_k', f'{rep}_N', 'rep_p1', 'rep_q', 'replicated']] \
        .rename(columns={f'{disc}_est': 'disc_est', f'{disc}_q': 'disc_q', f'{disc}_k': 'disc_k', f'{rep}_est': 'rep_est', f'{rep}_se': 'rep_se', f'{rep}_k': 'rep_k', f'{rep}_N': 'rep_N'})
    print(tag, 'discovered', len(D), 'replicated', int(D.replicated.sum()))
    des.append(D)
DES = pd.concat(des)
DES.to_csv(OUT + 'replication_designs.csv')
R['rep_A'] = R.index.isin(DES[(DES.design.str.startswith('A')) & DES.replicated].index)
R['rep_B'] = R.index.isin(DES[(DES.design.str.startswith('B')) & DES.replicated].index)
R['rep_C'] = R.index.isin(DES[(DES.design.str.startswith('C')) & DES.replicated].index)
R.to_csv(OUT + 'meta_all.csv')
print('all-cohort meta: k>=3 tested', (R.all_k >= 3).sum(), 'q<0.05', ((R.all_q < 0.05) & (R.all_k >= 3)).sum())
show = R[(R.all_k >= 3)].sort_values('all_p').head(60)
pd.set_option('display.width', 250)
print(show[['label', 'sub', 'all_k', 'all_N', 'all_est', 'all_p', 'all_q', 'all_I2', 'met_est', 'oth_est', 'oth_k', 'BM_est', 'rep_A', 'rep_B', 'rep_C']].round(3).to_string())
print(DES[DES.replicated][['label', 'design', 'disc_est', 'disc_q', 'rep_est', 'rep_k', 'rep_q']].round(3).to_string())
