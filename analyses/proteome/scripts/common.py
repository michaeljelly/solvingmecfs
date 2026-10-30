import numpy as np, pandas as pd
from scipy import stats
W = '/tmp/work_proteome/'; D = W + 'data/'; R = W + 'res/'

def bh(p):
    p = np.asarray(p, float); out = np.full(p.shape, np.nan); ok = ~np.isnan(p)
    if ok.sum() == 0: return out
    q = p[ok]; n = len(q); o = np.argsort(q); r = q[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]; tmp = np.empty(n); tmp[o] = np.minimum(r, 1); out[ok] = tmp
    return out

def load(ds):
    X = pd.read_parquet(D + f'{ds}_X.parquet'); m = pd.read_csv(D + f'{ds}_meta.csv', index_col=0)
    f = pd.read_csv(D + f'{ds}_feat.csv', index_col=0); m.index = m.index.astype(str); X.index = X.index.astype(str)
    return X, m, f

def gene_collapse(X, feat, min_frac=0.5):
    """z-score each feature, average features mapping to the same gene symbol (first symbol if several)."""
    ok = X.notna().mean() >= min_frac; X = X.loc[:, ok]
    Z = (X - X.mean()) / X.std(ddof=1)
    g = feat.loc[X.columns, 'gene'].astype(str).str.split(r'[ |;,]').str[0]
    g = g.where(~g.isin(['nan', 'None', '', 'NA']))
    Z = Z.loc[:, g.notna().values]; g = g[g.notna()]
    G = Z.T.groupby(g.values).mean().T
    return G

def ols_many(Y, design, test_col):
    """Y: samples x features DataFrame (NaN allowed -> per-feature fit with complete rows when needed).
    design: DataFrame of covariates incl. test_col. Returns beta, se, t, p, n for test_col.
    Y is standardised per feature so beta is in SD units."""
    Xd = design.copy(); Xd.insert(0, 'const', 1.0); Xd = Xd.astype(float)
    j = list(Xd.columns).index(test_col)
    res = {}
    Yv = Y.values.astype(float)
    full = ~np.isnan(Yv).any(0)
    out = pd.DataFrame(index=Y.columns, columns=['beta', 'se', 't', 'p', 'n'], dtype=float)
    A = Xd.values
    def fit(A, Yb):
        n, k = A.shape
        XtX_inv = np.linalg.pinv(A.T @ A); B = XtX_inv @ A.T @ Yb
        resid = Yb - A @ B; df = n - k; s2 = (resid ** 2).sum(0) / df
        se = np.sqrt(s2 * XtX_inv[j, j]); t = B[j] / se; p = 2 * stats.t.sf(np.abs(t), df)
        return B[j], se, t, p, n
    if full.any():
        Yb = Yv[:, full]; Yb = (Yb - Yb.mean(0)) / Yb.std(0, ddof=1)
        b, se, t, p, n = fit(A, Yb)
        out.loc[Y.columns[full], ['beta', 'se', 't', 'p']] = np.c_[b, se, t, p]; out.loc[Y.columns[full], 'n'] = n
    for i in np.where(~full)[0]:
        y = Yv[:, i]; ok = ~np.isnan(y)
        if ok.sum() < A.shape[1] + 5 or np.nanstd(y) == 0: continue
        Ab = A[ok]
        if np.linalg.matrix_rank(Ab) < Ab.shape[1]: continue
        yb = (y[ok] - y[ok].mean()) / y[ok].std(ddof=1)
        b, se, t, p, n = fit(Ab, yb[:, None])
        out.iloc[i] = [b[0], se[0], t[0], p[0], n]
    out['q'] = bh(out.p.values)
    return out

def meta_re(d, se):
    """DerSimonian-Laird random effects + fixed effects. d, se arrays (NaN dropped)."""
    d = np.asarray(d, float); se = np.asarray(se, float); ok = ~(np.isnan(d) | np.isnan(se)); d = d[ok]; se = se[ok]
    k = len(d)
    if k == 0: return dict(k=0)
    w = 1 / se ** 2; fe = (w * d).sum() / w.sum(); fese = np.sqrt(1 / w.sum())
    Q = (w * (d - fe) ** 2).sum(); c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0, (Q - (k - 1)) / c) if k > 1 else 0
    wr = 1 / (se ** 2 + tau2); re = (wr * d).sum() / wr.sum(); rese = np.sqrt(1 / wr.sum())
    I2 = max(0, (Q - (k - 1)) / Q) if Q > 0 and k > 1 else 0
    return dict(k=k, fe=fe, fe_se=fese, fe_p=2 * stats.norm.sf(abs(fe / fese)), re=re, re_se=rese,
                re_p=2 * stats.norm.sf(abs(re / rese)), I2=I2, Qp=stats.chi2.sf(Q, k - 1) if k > 1 else np.nan)
