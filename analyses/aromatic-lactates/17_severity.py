"""Within-patient severity: do the replicated metabolites track how ill people are?
BioMapAI (person means over visits): composite severity (mean of 12 scaled scores, higher = worse), orthostatic intolerance (NLT),
physical functioning (RAND-36, sign flipped so higher = worse), energy/fatigue (flipped), cognitive efficiency (flipped).
Germain 2022: Bell disability score (flipped: higher = worse). Partial Spearman adjusted for age, sex, BMI.
Metabolome-wide: BioMapAI severity + Germain Bell, Fisher-z meta, BH. Also illness duration (short <4y vs long) in BioMapAI and men-only.
Outputs: severity_key_items.csv, severity_metabolome_wide.csv"""
import numpy as np, pandas as pd
from scipy import stats
from mwlib import *

S, M = load()


def pspear(Y, x, Z):
    """Y: DataFrame, x: Series, Z: covariate DataFrame -> r, p, n (rank-residualised)."""
    ok = x.notna() & Z.notna().all(axis=1)
    Yr = Y[ok].rank(); xr = x[ok].rank(); Zc = np.column_stack([np.ones(ok.sum()), Z[ok].values.astype(float)])
    ry = Yr.values - Zc @ np.linalg.lstsq(Zc, Yr.values, rcond=None)[0]; rx = xr.values - Zc @ np.linalg.lstsq(Zc, xr.values, rcond=None)[0]
    r = (ry * rx[:, None]).sum(0) / np.sqrt((ry ** 2).sum(0) * (rx ** 2).sum())
    n = ok.sum(); k = Zc.shape[1] - 1
    t = r * np.sqrt((n - 2 - k) / (1 - r ** 2)); p = 2 * stats.t.sf(np.abs(t), n - 2 - k)
    return pd.DataFrame({'r': r, 'p': p, 'n': n}, index=Y.columns)


C = S['BioMapAI']['C']; pt = C.case == 1
SC = {'severity (composite, higher=worse)': C['severity'], 'orthostatic intolerance NLT': C['NLT.Orthostatic Intolerance'],
      'physical function (flipped)': -C['RAND36.Physical Functioning'], 'energy/fatigue (flipped)': -C['RAND36.Energy/Fatigue'],
      'cognitive efficiency (flipped)': -C['DANA.Cognitive Efficiency']}
rows = []
Yb = build_items('BioMapAI', S, M)
for nm, x in SC.items():
    for sub, m in [('patients', pt), ('male patients', pt & (C.male == 1))]:
        cov = ['age', 'male', 'bmi'] if sub == 'patients' else ['age', 'bmi']
        r = pspear(Yb[m], x[m], C.loc[m, cov]); r['cohort'] = 'BioMapAI'; r['score'] = nm; r['subset'] = sub; r['item'] = r.index; rows.append(r)
# duration
d = ols_effect(Yb[pt], C[pt].assign(long=1 - C[pt].short), 'long', ['age', 'male', 'bmi'])
d['cohort'] = 'BioMapAI'; d['score'] = 'long (>4y) vs short illness, SD'; d['subset'] = 'patients'; d['item'] = d.index; d = d.rename(columns={'est': 'r'}); rows.append(d)
Cg = S['Germain2022']['C']; pg = Cg.case == 1
Yg = build_items('Germain2022', S, M)
for sub, m in [('patients', pg), ('male patients', pg & (Cg.male == 1))]:
    cov = ['age', 'male', 'bmi'] if sub == 'patients' else ['age', 'bmi']
    r = pspear(Yg[m], -Cg.bell[m], Cg.loc[m, cov]); r['cohort'] = 'Germain2022'; r['score'] = 'Bell (flipped, higher=worse)'; r['subset'] = sub; r['item'] = r.index; rows.append(r)
K = pd.concat(rows, ignore_index=True)
K['q_within_score'] = K.groupby(['cohort', 'score', 'subset']).p.transform(bh)
K.to_csv(OUT + 'severity_key_items.csv', index=False)
pd.set_option('display.width', 250)
P = K[K.subset == 'patients'].pivot_table(index='item', columns='score', values='r').round(2)
print(P.to_string())
print('nominal p<0.05 (patients):'); print(K[(K.p < 0.05) & (K.subset == 'patients')][['cohort', 'score', 'item', 'r', 'p', 'q_within_score', 'n']].round(3).to_string())
print('men:'); print(K[(K.subset == 'male patients')].pivot_table(index='item', columns='score', values='r').round(2).to_string())
print('n male patients BioMapAI', int((pt & (C.male == 1)).sum()), 'Germain', int((pg & (Cg.male == 1)).sum()))

# metabolome-wide severity (BioMapAI composite) + Germain Bell, Fisher-z meta over harmonised features
rb = pspear(S['BioMapAI']['L'][pt], C.severity[pt], C.loc[pt, ['age', 'male', 'bmi']])
rg = pspear(S['Germain2022']['L'][pg], -Cg.bell[pg], Cg.loc[pg, ['age', 'male', 'bmi']])
mb = M[M.study == 'BioMapAI'].set_index('feature').hkey; mg = M[M.study == 'Germain2022'].set_index('feature').hkey
rb['hkey'] = mb.reindex(rb.index).values; rg['hkey'] = mg.reindex(rg.index).values
J = rb.dropna(subset=['hkey']).set_index('hkey').join(rg.dropna(subset=['hkey']).set_index('hkey'), lsuffix='_BM', rsuffix='_G22', how='left')
zb = np.arctanh(J.r_BM); sb = 1 / np.sqrt(J.n_BM - 6); zg = np.arctanh(J.r_G22); sg = 1 / np.sqrt(J.n_G22 - 6)
both = J.r_G22.notna()
mz = np.where(both, (zb / sb ** 2 + zg / sg ** 2) / (1 / sb ** 2 + 1 / sg ** 2), zb); ms = np.where(both, np.sqrt(1 / (1 / sb ** 2 + 1 / sg ** 2)), sb)
J['r_meta'] = np.tanh(mz); J['p_meta'] = 2 * stats.norm.sf(np.abs(mz / ms)); J['q_meta'] = bh(J.p_meta.values)
J['q_BM'] = bh(J.p_BM.values)
lab = M.drop_duplicates('hkey').set_index('hkey').label
J['label'] = lab.reindex(J.index).values
J['same_dir'] = np.sign(J.r_BM) == np.sign(J.r_G22)
J = J.sort_values('p_meta'); J.to_csv(OUT + 'severity_metabolome_wide.csv')
print('metabolome-wide severity: features', len(J), 'q_meta<0.05', int((J.q_meta < 0.05).sum()), 'q_BM<0.05', int((J.q_BM < 0.05).sum()),
      'replicate (q_BM<0.05 & G22 same dir p<0.05):', int(((J.q_BM < 0.05) & J.same_dir & (J.p_G22 < 0.05)).sum()))
print(J[['label', 'r_BM', 'p_BM', 'q_BM', 'r_G22', 'p_G22', 'r_meta', 'q_meta']].head(25).round(3).to_string())
