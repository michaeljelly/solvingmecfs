"""Sensitivity for the pre-specified key items: filtration-marker adjustment, full confounder model (M3),
sex-stratified effects and sex x case interaction, leave-one-cohort-out, per-cohort estimates.
Output: key_items_per_study.csv, key_items_sensitivity.csv"""
import pickle
import numpy as np, pandas as pd
from mwlib import *

S, M = load()
CADJ = pickle.load(open(W + 'C_adj.pkl', 'rb'))
rows = []
for st in S:
    Y = build_items(st, S, M)
    if Y.empty: continue
    C, MOD = CADJ[st]; C = C.copy()
    f = filtration(st, S, M)
    base = covs(st, C)
    specs = [('base', base, None), ('M3_caf_drug_diet_ket', [c for c in MOD['M3'][0] if C[c].std() > 0], None)]
    if f is not None:
        C['filt'] = f; specs.append(('+filtration', base + ['filt'], None))
    if 'male' in C:
        specs += [('men', [c for c in base if c != 'male'], C.male == 1), ('women', [c for c in base if c != 'male'], C.male == 0)]
        C['caseXmale'] = C.case * C.male
    for tag, cv, keep in specs:
        idx = C.index if keep is None else C.index[keep.values]
        Cc = C.loc[idx]
        if Cc.case.sum() < 4 or (1 - Cc.case).sum() < 4: continue
        cv = [c for c in cv if Cc[c].std() > 0]
        r = ols_effect(Y.loc[idx], Cc, 'case', cv)
        r['study'] = st; r['model'] = tag; r['item'] = r.index; r['n_case'] = int(Cc.case.sum()); r['n_ctrl'] = int((1 - Cc.case).sum()); rows.append(r)
    if 'male' in C and C.male.std() > 0:
        r = ols_effect(Y, C, 'caseXmale', base + ['case']); r['study'] = st; r['model'] = 'interaction caseXmale'; r['item'] = r.index
        r['n_case'] = int(C.case.sum()); r['n_ctrl'] = int((1 - C.case).sum()); rows.append(r)
# Hoel summary for single labels
h = pd.read_csv(W + 'hoel_summary.csv'); hm = dict(zip(M[M.study == 'Hoel2021'].label, M[M.study == 'Hoel2021'].feature))
for name, (kind, spec) in KEY.items():
    if kind == 'lab' and spec in hm:
        rr = h[h.Met == hm[spec]]
        if len(rr):
            d, se = d_from_fc_p(float(rr.FC.iloc[0]), float(rr.p.iloc[0]))
            rows.append(pd.DataFrame([dict(est=d, se=se, p=float(rr.p.iloc[0]), n=118, study='Hoel2021', model='base', item=name, n_case=83, n_ctrl=35)]))
A = pd.concat(rows, ignore_index=True); A.to_csv(OUT + 'key_items_per_study.csv', index=False)
out = []
for (it, mdl), g in A.groupby(['item', 'model']):
    r = meta(g.est, g.se); r.update(item=it, model=mdl, N=int((g.n_case + g.n_ctrl).sum()), n_neg=int((g.est < 0).sum()))
    if mdl == 'base' and len(g) >= 3:
        loo = [meta(g.drop(i).est, g.drop(i).se) for i in g.index]
        r['loo_max_p'] = max(x['p'] for x in loo); r['loo_est_range'] = f"{min(x['est'] for x in loo):+.2f} to {max(x['est'] for x in loo):+.2f}"
        g2 = g[~g.study.isin(['Germain2018', 'Germain2020'])]; r2 = meta(g2.est, g2.se); r['noHanson_est'] = r2['est']; r['noHanson_p'] = r2['p']
        g3 = g[g.study.map(FAMILY) == 'Other']; r3 = meta(g3.est, g3.se); r['nonMetabolon_est'] = r3['est']; r['nonMetabolon_p'] = r3['p']; r['nonMetabolon_k'] = r3['k']
    out.append(r)
R = pd.DataFrame(out)
W_ = R.pivot(index='item', columns='model', values='est').add_suffix('_est').join(R.pivot(index='item', columns='model', values='p').add_suffix('_p')).join(R.pivot(index='item', columns='model', values='k').add_suffix('_k'))
extra = R[R.model == 'base'].set_index('item')[['N', 'I2', 'n_neg', 'loo_max_p', 'loo_est_range', 'noHanson_est', 'noHanson_p', 'nonMetabolon_est', 'nonMetabolon_p', 'nonMetabolon_k']]
W_ = W_.join(extra)
W_['base_q'] = bh(W_['base_p'].values)
W_ = W_.sort_values('base_p')
W_.to_csv(OUT + 'key_items_sensitivity.csv')
pd.set_option('display.width', 300); pd.set_option('display.max_colwidth', 36)
print(W_[['base_k', 'N', 'base_est', 'base_p', 'base_q', 'I2', 'loo_max_p', 'noHanson_est', 'nonMetabolon_est', 'nonMetabolon_k', 'nonMetabolon_p', '+filtration_est', '+filtration_p', 'M3_caf_drug_diet_ket_est', 'M3_caf_drug_diet_ket_p']].round(3).to_string())
print(W_[['men_est', 'men_p', 'men_k', 'women_est', 'women_p', 'interaction caseXmale_est', 'interaction caseXmale_p']].round(3).to_string())
