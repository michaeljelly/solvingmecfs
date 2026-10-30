"""Shared helpers: vectorised OLS case effect, DerSimonian-Laird meta-analysis, BH-FDR, study covariates."""
import pickle
import numpy as np, pandas as pd
from scipy import stats

W = '/tmp/work_mw/'
OUT = '/mnt/project-files/mecfs-research/deep2/metabolome_wide/'
COV = {'BioMapAI': ['age', 'male', 'bmi', 'bmi_missing'], 'Germain2022': ['age', 'male', 'bmi'],
       'Che2022': ['age', 'male', 'bmi', 'site_*'], 'NagySzakal2018': ['age_over50', 'male', 'bmi'],
       'Naviaux2016': ['male'], 'Germain2018': [], 'Germain2020': [], 'NCNP2026': [], 'Armstrong2015': []}
FAMILY = {'BioMapAI': 'Metabolon', 'Germain2022': 'Metabolon', 'Germain2018': 'Metabolon', 'Germain2020': 'Metabolon', 'Hoel2021': 'Metabolon',
          'Che2022': 'Other', 'NagySzakal2018': 'Other', 'Naviaux2016': 'Other', 'NCNP2026': 'Other', 'Armstrong2015': 'Other'}
HOEL_N = (83, 35)


def load():
    S = pickle.load(open(W + 'studies_full.pkl', 'rb'))
    M = pickle.load(open(W + 'hmap.pkl', 'rb'))
    return S, M


def covs(study, C, extra=()):
    cols = []
    for c in COV[study]:
        cols += [x for x in C.columns if x.startswith('site_')] if c == 'site_*' else [c]
    return cols + [c for c in extra if c not in cols]


def ols_effect(Y, C, test='case', cov=()):
    """Y: n x p DataFrame (already any scale; z-scored here). Returns DataFrame est, se, p for `test` column, per feature."""
    Y = Y.loc[:, Y.std() > 1e-12]
    Z = (Y - Y.mean()) / Y.std()
    X = pd.concat([pd.Series(1.0, index=C.index, name='const'), C[[test] + list(cov)].astype(float)], axis=1)
    ok = X.notna().all(axis=1) & Z.notna().all(axis=1)
    X, Z = X[ok].values, Z[ok].values
    n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Z
    R = Z - X @ B
    s2 = (R ** 2).sum(0) / (n - k)
    se = np.sqrt(s2 * XtXi[1, 1])
    t = B[1] / se
    p = 2 * stats.t.sf(np.abs(t), n - k)
    return pd.DataFrame({'est': B[1], 'se': se, 'p': p, 'n': n}, index=Y.columns)


def ols_masked(Y, C, test='case', cov=()):
    """Per-feature OLS allowing missing values in Y (slower)."""
    rows = {}
    for c in Y.columns:
        y = Y[c]; m = y.notna() & C[[test] + list(cov)].notna().all(axis=1)
        if m.sum() < 10 or y[m].std() == 0: continue
        r = ols_effect(y[m].to_frame(), C[m], test, cov)
        rows[c] = r.iloc[0]
    return pd.DataFrame(rows).T


def d_from_fc_p(fc, p, n1=HOEL_N[0], n0=HOEL_N[1]):
    vd = 1 / n1 + 1 / n0
    z = stats.norm.isf(np.clip(p, 1e-300, 1) / 2) * np.sign(np.log(fc))
    d = z * np.sqrt(vd)
    return d, np.sqrt(vd + d ** 2 / (2 * (n1 + n0)))


def meta(e, s):
    e, s = np.asarray(e, float), np.asarray(s, float)
    k = len(e)
    if k == 0: return dict(est=np.nan, se=np.nan, p=np.nan, I2=np.nan, k=0)
    w = 1 / s ** 2; fe = np.sum(w * e) / w.sum(); Q = np.sum(w * (e - fe) ** 2)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (s ** 2 + tau2); re_ = np.sum(wr * e) / wr.sum(); rse = np.sqrt(1 / wr.sum())
    I2 = max(0, (Q - (k - 1)) / Q) if (k > 1 and Q > 0) else 0
    return dict(est=re_, se=rse, p=2 * stats.norm.sf(abs(re_ / rse)), I2=I2, k=k)


def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); m = ~np.isnan(p)
    if m.sum() == 0: return q
    pv = p[m]; o = np.argsort(pv); r = pv[o] * m.sum() / (np.arange(m.sum()) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]; qq = np.empty_like(r); qq[o] = np.minimum(r, 1); q[m] = qq
    return q


# ---------------------------------------------------------------- pre-specified key items (chosen from 12/14 discovery output)
AROM = ['phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'indolelactate', 'imidazole lactate']
COLL = ['trans-4-hydroxyproline', '5-(galactosylhydroxy)-lysine', '4-hydroxyglutamate', '4-hydroxy-2-oxoglutaric acid', '5-hydroxylysine']
LCAC = ['palmitoylcarnitine (C16)', 'stearoylcarnitine (C18)', 'oleoylcarnitine (C18:1)', 'linoleoylcarnitine (C18:2)*', 'myristoylcarnitine (C14)', 'arachidonoylcarnitine (C20:4)']
FILT = ['pseudouridine', 'C-glycosyltryptophan', 'N-formylmethionine', 'creatinine']
KEY = {
    'MODULE aromatic lactates': ('mod', AROM),
    'phenyllactate (PLA)': ('lab', 'phenyllactate (PLA)'), '3-(4-hydroxyphenyl)lactate': ('lab', '3-(4-hydroxyphenyl)lactate'),
    'indolelactate': ('lab', 'indolelactate'), 'imidazole lactate': ('lab', 'imidazole lactate'), 'vanillactate': ('lab', 'vanillactate'),
    'alpha-hydroxyisovalerate': ('lab', 'alpha-hydroxyisovalerate'), 'alpha-hydroxyisocaproate': ('lab', 'alpha-hydroxyisocaproate'),
    '2-hydroxy-4-(methylthio)butanoate': ('lab', '2-hydroxy-4-(methylthio)butanoic acid'),
    'indolepropionate (IPA)': ('lab', 'indolepropionate'), 'hippurate': ('lab', 'hippurate'),
    'MODULE collagen turnover': ('mod', COLL), 'trans-4-hydroxyproline': ('lab', 'trans-4-hydroxyproline'),
    '5-(galactosylhydroxy)-lysine': ('lab', '5-(galactosylhydroxy)-lysine'), '4-hydroxyglutamate': ('lab', '4-hydroxyglutamate'),
    'methionine sulfoxide / methionine': ('ratio', (['methionine sulfoxide'], ['methionine'])),
    'PC(14:0/16:0)': ('lab', '1-myristoyl-2-palmitoyl-GPC (14:0/16:0)'),
    'plasmalogen PC(P-16:0/18:2)': ('lab', '1-(1-enyl-palmitoyl)-2-linoleoyl-GPC (P-16:0/18:2)*'),
    'glucuronate': ('lab', 'glucuronate'), '5-oxoproline / glutamate': ('ratio', (['5-oxoproline'], ['glutamate'])),
    'kynurenine': ('lab', 'kynurenine'), 'carnitine': ('lab', 'carnitine'), 'tiglylcarnitine (C5:1-DC)': ('lab', 'tiglylcarnitine (C5:1-DC)'),
    'MODULE long-chain acylcarnitines': ('mod', LCAC), '1,5-anhydroglucitol': ('lab', '1,5-anhydroglucitol (1,5-AG)'),
    'MTA': ('lab', '5-methylthioadenosine (MTA)'), 'phenylacetylglutamine': ('lab', 'phenylacetylglutamine'),
    'TMAO / carnitine': ('ratio', (['trimethylamine N-oxide'], ['carnitine'])),
    'dihydroceramide (d18:0/16:0)': ('lab', 'N-palmitoyl-sphinganine (d18:0/16:0)'),
    'caffeine (lifestyle reference)': ('lab', 'caffeine'),
}


def build_items(st, S, M, items=KEY, min_mod=2):
    L = S[st]['L']; fm = dict(zip(M[M.study == st].label, M[M.study == st].feature))
    g = lambda l: fm[l] if (l in fm and fm[l] in L.columns) else None
    Y = {}
    for name, (kind, spec) in items.items():
        if kind == 'lab':
            c = g(spec)
            if c: Y[name] = L[c]
        elif kind == 'ratio':
            nc = [g(x) for x in spec[0]]; dc = [g(x) for x in spec[1]]
            if all(nc) and all(dc): Y[name] = np.log2(np.exp2(L[nc]).sum(axis=1)) - np.log2(np.exp2(L[dc]).sum(axis=1))
        elif kind == 'mod':
            cs = [c for c in (g(x) for x in spec) if c]
            if len(cs) >= min_mod:
                Z = (L[cs] - L[cs].mean()) / L[cs].std(); Y[name] = Z.mean(axis=1)
    return pd.DataFrame(Y, index=L.index)


def filtration(st, S, M):
    L = S[st]['L']; fm = dict(zip(M[M.study == st].label, M[M.study == st].feature))
    cs = [fm[l] for l in FILT if l in fm and fm[l] in L.columns]
    if len(cs) < 2: return None
    Z = (L[cs] - L[cs].mean()) / L[cs].std(); return Z.mean(axis=1)
