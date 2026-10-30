"""Task 1: how well does the gut (species CLR or KEGG log-CPM) predict each plasma metabolite / immune feature,
separately in controls and patients, with person-level sampling (one random visit per person per repeat),
patients subsampled to the control n, 5-fold CV, ridge (alpha tuned per target by GCV inside the training fold),
covariates (age, sex, BMI) regressed out of the target inside the fold.
Nulls: (a) Y rows shuffled vs X (no-signal null), (b) label swap: pool matched controls+patients, split at random
into two pseudo-groups of the same size (tests whether the control-patient difference exceeds chance).
usage: python3 01_predictability.py <target: metab|immune> <pred: species|kegg> <pgroup: all|Short|Long> [R] [K]
Nulls use a FIXED person-level permutation/partition per null draw (K draws x R0=3 visit resamples), so person-level
sampling noise is represented in the null.
"""
import sys, time
from cc import *
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

target, pred, pgroup = sys.argv[1], sys.argv[2], sys.argv[3]
R = int(sys.argv[4]) if len(sys.argv) > 4 else 30
rng = np.random.default_rng(1)

if target == 'metab':
    Y, ann = metab()
else:
    Y = immune()
if pred == 'species':
    X, _ = species_clr()
else:
    X, raw = kegg_logcpm()
    X = X.loc[:, (raw > 0).mean() > 0.2]

common = X.index.intersection(Y.index).intersection(meta.index)
M = meta.loc[common]
X, Y = X.loc[common], Y.loc[common]
Y = Y.loc[:, Y.std() > 0]
ctrl = M[M.y == 0]
pat = M[M.y == 1] if pgroup == 'all' else M[(M.y == 1) & (M.illness_duration == pgroup)]
cs, ps = ctrl.subject.unique(), pat.subject.unique()
n = min(len(cs), len(ps))
print(target, pred, pgroup, 'persons ctrl', len(cs), 'pat', len(ps), '-> n', n, 'X', X.shape, 'Y', Y.shape, flush=True)
alphas = np.logspace(-1, 6, 15)


def one_visit(Mg, subs):
    idx = []
    for s in subs:
        v = Mg.index[Mg.subject == s]
        idx.append(v[rng.integers(len(v))])
    return idx


def cv_r2(idx, shuffle=False, extra_cov=None):
    Xi = X.loc[idx].values
    Yi = Y.loc[idx].values.copy()
    C = M.loc[idx, COV].values.astype(float)
    if extra_cov is not None:
        C = np.column_stack([C, extra_cov])
    if shuffle:
        Yi = Yi[rng.permutation(len(Yi))]
        C = C  # covariates stay with X side; target is scrambled relative to both
    oof = np.zeros_like(Yi); yres = np.zeros_like(Yi)
    for tr, te in KFold(5, shuffle=True, random_state=int(rng.integers(1e9))).split(Xi):
        Ct = np.column_stack([np.ones(len(tr)), C[tr]]); Ce = np.column_stack([np.ones(len(te)), C[te]])
        b, *_ = np.linalg.lstsq(Ct, Yi[tr], rcond=None)
        ytr = Yi[tr] - Ct @ b; yte = Yi[te] - Ce @ b
        mu, sd = ytr.mean(0), ytr.std(0) + 1e-12
        ytr = np.clip((ytr - mu) / sd, -4, 4); yte = np.clip((yte - mu) / sd, -4, 4)
        xm, xs = Xi[tr].mean(0), Xi[tr].std(0) + 1e-9
        xtr = (Xi[tr] - xm) / xs; xte = (Xi[te] - xm) / xs
        m = RidgeCV(alphas=alphas, alpha_per_target=True).fit(xtr, ytr)
        oof[te] = m.predict(xte); yres[te] = yte
    sse = ((yres - oof) ** 2).sum(0); sst = ((yres - yres.mean(0)) ** 2).sum(0)
    r2 = 1 - sse / sst
    yc = yres - yres.mean(0); oc = oof - oof.mean(0)
    r = (yc * oc).sum(0) / np.sqrt((yc ** 2).sum(0) * (oc ** 2).sum(0) + 1e-12)
    return r2, r


K = int(sys.argv[5]) if len(sys.argv) > 5 else 30   # null partitions
R0 = 3                                              # visit-resampling repeats per null partition (fewer than R -> conservative)


def vis(Mg, subs):
    return one_visit(Mg, subs)


res = {k: [] for k in ['c_r2', 'p_r2', 'c_r', 'p_r']}
t0 = time.time()
for rep in range(R):
    cidx = vis(ctrl, rng.choice(cs, n, replace=False))
    pidx = vis(pat, rng.choice(ps, n, replace=False))
    a, b = cv_r2(cidx); res['c_r2'].append(a); res['c_r'].append(b)
    a, b = cv_r2(pidx); res['p_r2'].append(a); res['p_r'].append(b)
    if rep % 10 == 0:
        print(rep, round(time.time() - t0), 'median R2 ctrl %.3f pat %.3f' % (np.median(a), np.median(res['c_r2'][-1])), flush=True)


def cv_generic(Xi, Yi, C):
    oof = np.zeros_like(Yi); yres = np.zeros_like(Yi)
    for tr, te in KFold(5, shuffle=True, random_state=int(rng.integers(1e9))).split(Xi):
        Ct = np.column_stack([np.ones(len(tr)), C[tr]]); Ce = np.column_stack([np.ones(len(te)), C[te]])
        b, *_ = np.linalg.lstsq(Ct, Yi[tr], rcond=None)
        ytr = Yi[tr] - Ct @ b; yte = Yi[te] - Ce @ b
        mu, sd = ytr.mean(0), ytr.std(0) + 1e-12
        ytr = np.clip((ytr - mu) / sd, -4, 4); yte = np.clip((yte - mu) / sd, -4, 4)
        xm, xs = Xi[tr].mean(0), Xi[tr].std(0) + 1e-9
        m = RidgeCV(alphas=alphas, alpha_per_target=True).fit((Xi[tr] - xm) / xs, ytr)
        oof[te] = m.predict((Xi[te] - xm) / xs); yres[te] = yte
    return 1 - ((yres - oof) ** 2).sum(0) / ((yres - yres.mean(0)) ** 2).sum(0)


# (a) shuffle null: fixed person-level permutation (Y of person pi(i) paired with X of person i) across R0 visit draws
nullc, nullp = [], []
for k in range(K):
    for Mg, subs, store in [(ctrl, cs, nullc), (pat, ps, nullp)]:
        sel = rng.choice(subs, n, replace=False); perm = rng.permutation(sel)
        acc = []
        for r in range(R0):
            xi = vis(Mg, sel); yi = vis(Mg, perm)
            acc.append(cv_generic(X.loc[xi].values, Y.loc[yi].values, M.loc[yi, COV].values.astype(float)))
        store.append(np.mean(acc, 0))
    if k % 10 == 0:
        print('null', k, round(time.time() - t0), flush=True)
# (b) swap null: fixed random split of the pooled (ctrl + subsampled patients) persons into two pseudo-groups
swd = []
Mall = pd.concat([ctrl, pat])
for k in range(K):
    pool = np.r_[rng.choice(cs, n, replace=False), rng.choice(ps, n, replace=False)]
    grp = pd.Series(np.r_[np.zeros(n), np.ones(n)], index=pool)
    perm = rng.permutation(pool); g1, g2 = perm[:n], perm[n:]
    acc = []
    for r in range(R0):
        d = []
        for g in (g1, g2):
            vi = vis(Mall, g)
            C = np.column_stack([M.loc[vi, COV].values.astype(float), grp.loc[g].values])
            d.append(cv_generic(X.loc[vi].values, Y.loc[vi].values, C))
        acc.append(d[0] - d[1])
    swd.append(np.mean(acc, 0))
    if k % 10 == 0:
        print('swap', k, round(time.time() - t0), flush=True)
nullc, nullp, swd = np.array(nullc), np.array(nullp), np.array(swd)
np.savez_compressed(OUT + f'pred_nulls_{target}_{pred}_{pgroup}.npz', nullc=nullc, nullp=nullp, swd=swd,
                    c_r2=np.array(res['c_r2']), p_r2=np.array(res['p_r2']), cols=np.array(Y.columns))

A = {k: np.array(v) for k, v in res.items()}  # R x targets
out = pd.DataFrame(index=Y.columns)
out['ctrl_R2'] = A['c_r2'].mean(0); out['pat_R2'] = A['p_r2'].mean(0)
out['ctrl_r'] = A['c_r'].mean(0); out['pat_r'] = A['p_r'].mean(0)
out['diff_R2'] = out.ctrl_R2 - out.pat_R2
# null of per-target mean R2 (pooled over targets): each target's mean of shuffled repeats
pooled_null = np.r_[nullc.ravel(), nullp.ravel()]
out['ctrl_p'] = [(1 + (pooled_null >= v).sum()) / (1 + len(pooled_null)) for v in out.ctrl_R2]
out['pat_p'] = [(1 + (pooled_null >= v).sum()) / (1 + len(pooled_null)) for v in out.pat_R2]
from scipy.stats import norm
mu0, sd0 = pooled_null.mean(), pooled_null.std()
# empirical p is floored at 1/(1+2*targets); use a Gaussian tail approximation of the pooled null for q-values
from scipy.stats import t as tdist
# per-target null (K draws of that target's own null) -> t-approximation; robust to target-specific R2 spread
out['ctrl_pz'] = tdist.sf((out.ctrl_R2 - nullc.mean(0)) / (nullc.std(0, ddof=1) * np.sqrt(1 + 1 / K)), K - 1)
out['pat_pz'] = tdist.sf((out.pat_R2 - nullp.mean(0)) / (nullp.std(0, ddof=1) * np.sqrt(1 + 1 / K)), K - 1)
out['ctrl_q'] = bh(out.ctrl_pz); out['pat_q'] = bh(out.pat_pz)
swd_t = swd  # K x targets
swd = swd.ravel()
out['swap_sd'] = swd.std()
out['diff_p'] = [(1 + (np.abs(swd) >= abs(v)).sum()) / (1 + len(swd)) for v in out.diff_R2]
out['diff_pz'] = 2 * tdist.sf(np.abs(out.diff_R2 - swd_t.mean(0)) / (swd_t.std(0, ddof=1) * np.sqrt(1 + 1 / K)), K - 1)
out['diff_q'] = bh(out.diff_pz)
# per-repeat sd of the difference (Monte-Carlo stability)
out['diff_R2_sd_rep'] = (A['c_r2'] - A['p_r2']).std(0)
if target == 'metab':
    out['super'] = ann.SUPER_PATHWAY.reindex(out.index).values
    out['sub'] = ann.SUB_PATHWAY.reindex(out.index).values
tag = f'{target}_{pred}_{pgroup}'
out.to_csv(OUT + f'pred_{tag}.csv')
# global summary with swap null on the global statistic
glob = dict(tag=tag, n_persons_per_group=n, R=R, n_targets=len(out),
            ctrl_median_R2=np.median(out.ctrl_R2), pat_median_R2=np.median(out.pat_R2),
            ctrl_mean_R2=out.ctrl_R2.mean(), pat_mean_R2=out.pat_R2.mean(),
            null_q99=np.quantile(pooled_null, .99),
            ctrl_n_q05=(out.ctrl_q < .05).sum(), pat_n_q05=(out.pat_q < .05).sum(),
            both_q05=((out.ctrl_q < .05) & (out.pat_q < .05)).sum(),
            ctrl_only=((out.ctrl_q < .05) & (out.pat_q >= .05)).sum(),
            pat_only=((out.pat_q < .05) & (out.ctrl_q >= .05)).sum(),
            n_diff_q05=(out.diff_q < .05).sum(),
            corr_ctrl_pat_R2=np.corrcoef(out.ctrl_R2, out.pat_R2)[0, 1])
# global test: mean R2 diff per repeat vs swap mean diff per repeat
gd = out.ctrl_R2.mean() - out.pat_R2.mean(); gs = swd_t.mean(1)
glob['global_meanR2_diff'] = gd; glob['swap_null_sd'] = gs.std()
glob['global_p_swap'] = (1 + (np.abs(gs - gs.mean()) >= abs(gd - gs.mean())).sum()) / (1 + len(gs))
# number of targets predictable (q<.05) difference vs swap null is not computed (would need per-partition q); use mean R2 among targets predictable in either group
pr = (out.ctrl_q < .05) | (out.pat_q < .05)
glob['n_predictable_either'] = pr.sum()
glob['meanR2_ctrl_among_predictable'] = out.ctrl_R2[pr].mean(); glob['meanR2_pat_among_predictable'] = out.pat_R2[pr].mean()
pd.Series(glob).to_csv(OUT + f'pred_{tag}_global.csv')
print(pd.Series(glob).to_string())
top = out.sort_values('ctrl_R2', ascending=False).head(12)[['ctrl_R2', 'pat_R2', 'ctrl_q', 'pat_q', 'diff_q']]
print(top.round(3).to_string())
