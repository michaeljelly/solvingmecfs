"""Sensitivity analyses for the leads that showed any external signal.
a) adjust for coffee markers (caffeine/paraxanthine) where measured per sample
b) leave-one-study-out and drop-summary-level (Hoel) meta
c) men only
d) Baraniuk 2025 CSF: disease x prior-exercise interaction (between-person)
Output: sensitivity.csv, exertion_csf.csv
"""
import sys, pickle, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
warnings.filterwarnings('ignore')
OUT = '/mnt/project-files/mecfs-research/downloads/mwb/'
sys.path.insert(0, OUT)
exec(open(OUT + '02_effects.py').read().split('# ---------------------------------------------------------------- Hoel')[0].split('rows = []')[0])  # imports + bh
from importlib import util
S = pickle.load(open(OUT + 'work_studies.pkl', 'rb'))
E = pd.read_csv(OUT + 'per_study_effects.csv')


def meta(e, s):
    e, s = np.asarray(e, float), np.asarray(s, float)
    w = 1 / s ** 2; fe = np.sum(w * e) / w.sum(); k = len(e); Q = np.sum(w * (e - fe) ** 2)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (s ** 2 + tau2); m = np.sum(wr * e) / wr.sum(); se = np.sqrt(1 / wr.sum())
    return m, se, 2 * stats.norm.sf(abs(m / se))


def fit(st, f, extra=(), sub=None):
    d = S[st]; X = d['X']; C = d['C'].copy()
    if f not in X: return None
    Z = (X - X.mean()) / X.std()
    covs = [c for c in C.columns if c not in ('case', 'person')]
    for c in covs:
        if C[c].nunique() > 2: C[c] = (C[c] - C[c].mean()) / C[c].std()
    D = C[['case'] + covs].astype(float)
    for e in extra:
        if e not in Z: return None
        D[e] = Z[e]
    if sub is not None:
        D = D[sub(C)]; D = D.loc[:, D.nunique() > 1] if len(D) else D
        if 'case' not in D: return None
    y = Z[f].loc[D.index]; ok = y.notna()
    if st == 'BioMapAI_Xiong2025':
        md = sm.MixedLM(y[ok], sm.add_constant(D[ok]), groups=C.loc[D.index[ok], 'person']).fit()
    else:
        md = sm.OLS(y[ok], sm.add_constant(D[ok])).fit()
    return md.params['case'], md.bse['case'], md.pvalues['case'], int(D.case[ok].sum()), int((1 - D.case[ok]).sum())


LEADS = ['IPA (indolepropionate)', 'indolelactate', 'kynurenine', 'kynurenine/tryptophan ratio', 'carnitine', 'COMP long-chain acylcarnitines', 'palmitate (16:0)', 'COMP primary bile acids', 'VMA (vanillylmandelate)']
BLOOD = [s for s, d in S.items() if d['tissue'] in ('plasma', 'serum')]
rows = []
for f in LEADS:
    ext = E[(E.feature == f) & E.tissue.isin(['plasma', 'serum']) & (E.study != 'BioMapAI_Xiong2025')]
    m, se, p = meta(ext.est, ext.se); rows.append(dict(feature=f, analysis='main external', k=len(ext), est=m, lo=m - 1.96 * se, hi=m + 1.96 * se, p=p))
    nh = ext[ext.study != 'Hoel2021_JCI_summary']
    m, se, p = meta(nh.est, nh.se); rows.append(dict(feature=f, analysis='external, per-sample data only (no Hoel)', k=len(nh), est=m, lo=m - 1.96 * se, hi=m + 1.96 * se, p=p))
    for st in ext.study:
        o = ext[ext.study != st]
        if len(o) >= 2:
            m, se, p = meta(o.est, o.se); rows.append(dict(feature=f, analysis=f'leave out {st}', k=len(o), est=m, lo=m - 1.96 * se, hi=m + 1.96 * se, p=p))
    # coffee adjustment (per-sample studies with caffeine or paraxanthine)
    rr = []
    for st in BLOOD:
        X = S[st]['X']
        adj = [c for c in ['paraxanthine', 'caffeine'] if c in X][:1]
        r = fit(st, f, adj)
        if r: rr.append((st, adj, *r))
    if rr:
        e = pd.DataFrame(rr, columns=['study', 'adj', 'est', 'se', 'p', 'nc', 'nh'])
        eo = e[e.study != 'BioMapAI_Xiong2025']
        m, se, p = meta(eo.est, eo.se)
        rows.append(dict(feature=f, analysis='external, adjusted for coffee marker (' + '; '.join(f"{s.split('_')[0]}:{'/'.join(a) or 'none'}" for s, a in zip(eo.study, eo.adj)) + ')', k=len(eo), est=m, lo=m - 1.96 * se, hi=m + 1.96 * se, p=p))
        b = e[e.study == 'BioMapAI_Xiong2025']
        if len(b): rows.append(dict(feature=f, analysis='BioMapAI adjusted for ' + '/'.join(b.adj.iloc[0]), k=1, est=b.est.iloc[0], lo=np.nan, hi=np.nan, p=b.p.iloc[0]))
    # men only
    rr = []
    for st in BLOOD:
        if 'male' not in S[st]['C']: continue
        r = fit(st, f, sub=lambda C: C.male == 1)
        if r and r[3] >= 5 and r[4] >= 5: rr.append((st, *r))
    if rr:
        e = pd.DataFrame(rr, columns=['study', 'est', 'se', 'p', 'nc', 'nh'])
        for lab, sub_ in [('men only, all cohorts incl. BioMapAI', e), ('men only, external', e[e.study != 'BioMapAI_Xiong2025'])]:
            if len(sub_) >= 2:
                m, se, p = meta(sub_.est, sub_.se)
                rows.append(dict(feature=f, analysis=f'{lab} (' + ', '.join(f"{s.split('_')[0]} {a}/{b}" for s, a, b in zip(sub_.study, sub_.nc, sub_.nh)) + ')', k=len(sub_), est=m, lo=m - 1.96 * se, hi=m + 1.96 * se, p=p))
R = pd.DataFrame(rows)
R.to_csv(OUT + 'sensitivity.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 120); pd.set_option('display.max_rows', 300)
print(R.round(3).to_string())

# ---------------------------------------------------------------- Baraniuk 2025 CSF: does prior exercise change patient-control difference?
d = S['Baraniuk2025_CSF']; X, C = d['X'], d['C']
rows = []
for f in X.columns:
    D = C[['case', 'postexercise', 'age', 'male']].astype(float).copy(); D['case_x_post'] = D.case * D.postexercise
    y = X[f]; ok = y.notna()
    r = sm.OLS(y[ok], sm.add_constant(D[ok])).fit()
    rows.append(dict(feature=f, case_nonexercise=r.params['case'], p_case=r.pvalues['case'], post_in_controls=r.params['postexercise'], p_post=r.pvalues['postexercise'],
                     interaction=r.params['case_x_post'], p_interaction=r.pvalues['case_x_post'], n=int(ok.sum())))
P = pd.DataFrame(rows); P.to_csv(OUT + 'exertion_csf.csv', index=False); print(P.round(3).to_string())
