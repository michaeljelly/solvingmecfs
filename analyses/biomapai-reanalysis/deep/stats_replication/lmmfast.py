"""Vectorised random-intercept linear mixed model (REML), fitted for many features at once.

Model per feature k:  y_k = C b_k + F_k g_k + u_person + e,  u ~ N(0, s2*lam), e ~ N(0, s2)
The variance ratio lam is chosen per feature by maximising the REML likelihood over a fine grid
of intra-class correlations (ICC = lam/(1+lam)); given lam the model is exact GLS. The
transformation v* = v - c_g * mean_g(v), c_g = 1 - 1/sqrt(1+n_g*lam) whitens the covariance.
Wald t-tests use df = n_persons - n_fixed (conservative for within-person terms).
Validated against statsmodels MixedLM (see validate_lmm.py).
"""
import numpy as np
from scipy import stats

ICC_GRID = np.r_[0.0, np.linspace(0.01, 0.97, 97)]


def _group_ops(groups):
    g, inv = np.unique(np.asarray(groups), return_inverse=True)
    ng = np.bincount(inv)
    return inv, ng


def _transform(A, inv, ng, c):
    """A: (n, m) -> A - c_g * groupmean(A). Works on 2-D arrays."""
    sums = np.zeros((len(ng),) + A.shape[1:])
    np.add.at(sums, inv, A)
    means = sums / ng.reshape((-1,) + (1,) * (A.ndim - 1))
    return A - (c[inv].reshape((-1,) + (1,) * (A.ndim - 1))) * means[inv]


def fit_many(Y, C, groups, L=None, F=None, grid=ICC_GRID):
    """Y: (n,K) outcomes (one column per feature) OR (n,) shared outcome when F given.
    C: (n,p) shared fixed design (include intercept).
    L: (r,p) contrast rows on C coefficients to report (used when F is None).
    F: (K,n,q) feature-specific design columns (effects reported are those of F).
    Returns dict with est (K,r|q), se, t, p, icc, df.
    """
    inv, ng = _group_ops(groups)
    n = len(inv); G = len(ng)
    C = np.asarray(C, float)
    p = C.shape[1]
    shared_y = F is not None and np.ndim(Y) == 1
    if F is not None:
        F = np.asarray(F, float)
        K, _, q = F.shape
        Yf = np.asarray(Y, float)
        if shared_y:
            Yf = np.repeat(Yf[:, None], K, axis=1)
    else:
        Yf = np.asarray(Y, float)
        K = Yf.shape[1]; q = 0
    ptot = p + q
    df = G - ptot
    r = L.shape[0] if F is None else q
    best_ll = np.full(K, -np.inf); est = np.zeros((K, r)); se = np.zeros((K, r)); icc = np.zeros(K)
    for rho in grid:
        lam = rho / (1 - rho)
        c = 1 - 1 / np.sqrt(1 + ng * lam)
        logdetV = np.sum(np.log1p(ng * lam))
        Cs = _transform(C, inv, ng, c)
        Ys = _transform(Yf, inv, ng, c)
        CtC = Cs.T @ Cs
        CtC_inv = np.linalg.inv(CtC)
        sign, logdetC = np.linalg.slogdet(CtC)
        B = CtC_inv @ (Cs.T @ Ys)              # (p,K)
        R = Ys - Cs @ B                        # residual on C
        if F is None:
            rss = np.einsum('ij,ij->j', R, R)
            s2 = rss / (n - ptot)
            ll = -0.5 * ((n - ptot) * np.log(s2) + logdetV + logdetC)
            e = (L @ B).T                      # (K,r)
            v = np.einsum('ij,jk,ik->i', L, CtC_inv, L)  # (r,)
            s = np.sqrt(np.outer(s2, v))
        else:
            Fs = _transform(F.transpose(1, 0, 2).reshape(n, K * q), inv, ng, c).reshape(n, K, q)
            # residualise F on C
            BF = CtC_inv @ (Cs.T @ Fs.reshape(n, K * q))
            Fr = Fs - (Cs @ BF).reshape(n, K, q)
            FtF = np.einsum('nkq,nkr->kqr', Fr, Fr)
            FtY = np.einsum('nkq,nk->kq', Fr, R)
            FtF_inv = np.linalg.inv(FtF + 1e-12 * np.eye(q))
            g = np.einsum('kqr,kr->kq', FtF_inv, FtY)
            rss = np.einsum('nk,nk->k', R, R) - np.einsum('kq,kq->k', g, FtY)
            s2 = rss / (n - ptot)
            _, ldF = np.linalg.slogdet(FtF + 1e-12 * np.eye(q))
            ll = -0.5 * ((n - ptot) * np.log(s2) + logdetV + logdetC + ldF)
            e = g
            s = np.sqrt(s2[:, None] * np.einsum('kqq->kq', FtF_inv))
        upd = ll > best_ll
        best_ll[upd] = ll[upd]; est[upd] = e[upd]; se[upd] = s[upd]; icc[upd] = rho
    t = est / se
    pv = 2 * stats.t.sf(np.abs(t), df)
    return dict(est=est, se=se, t=t, p=pv, icc=icc, df=df, ll=best_ll)


def bh(p):
    p = np.asarray(p, float); out = np.full_like(p, np.nan)
    ok = ~np.isnan(p); pp = p[ok]; m = len(pp)
    if m == 0: return out
    o = np.argsort(pp); ranked = pp[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(m); res[o] = np.minimum(q, 1); out[ok] = res
    return out
