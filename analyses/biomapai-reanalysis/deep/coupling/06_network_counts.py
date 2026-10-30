"""Paper-style 'network' test done honestly: number of significant gut-species x blood-feature Spearman correlations
(BH q<0.05 and |r|>0.35, the authors' thresholds) in controls vs patients, with patients subsampled to the control
person count, one value per person (mean over paired visits), covariates regressed out within group.
Null: random person-level split of pooled (controls + subsampled patients) into two pseudo-groups."""
import sys
from cc import *
from scipy import stats

rng = np.random.default_rng(11)
spc, _ = species_clr()
Xm, _ = metab()
IM = immune()


def corr_counts(A, B):
    ra = A.rank().values; rb = B.rank().values
    ra = (ra - ra.mean(0)) / ra.std(0); rb = (rb - rb.mean(0)) / rb.std(0)
    r = ra.T @ rb / len(ra)
    n = len(ra)
    t = r * np.sqrt((n - 2) / np.clip(1 - r ** 2, 1e-12, None))
    p = 2 * stats.t.sf(np.abs(t), n - 2)
    q = bh(p.ravel()).reshape(p.shape)
    return dict(n_q05=int((q < .05).sum()), n_q05_r35=int(((q < .05) & (np.abs(r) > .35)).sum()), mean_abs_r=float(np.abs(r).mean()))


rows = []
for blood_name, Bl in [('metabolome', Xm), ('immune', IM)]:
    for pg in ['all', 'Short', 'Long']:
        common = spc.index.intersection(Bl.index).intersection(meta.index)
        M = meta.loc[common]
        D = pd.concat([spc.loc[common].add_prefix('sp|'), Bl.loc[common].add_prefix('bl|')], axis=1)
        per = D.groupby(M.subject).mean()
        pm = M.groupby('subject')[COV + ['y']].first(); pm['dur'] = M.groupby('subject').illness_duration.first()
        cs = pm.index[pm.y == 0]; ps = pm.index[(pm.y == 1) & ((pm.dur == pg) | (pg == 'all'))]
        n = min(len(cs), len(ps))
        spcols = [c for c in per.columns if c.startswith('sp|')]; blcols = [c for c in per.columns if c.startswith('bl|')]

        def cnt(idx):
            R = resid(per.loc[idx], pm.loc[idx, COV]); return corr_counts(R[spcols], R[blcols])
        obs = []
        for rep in range(30):
            c = rng.choice(cs, n, replace=False); p = rng.choice(ps, n, replace=False)
            a, b = cnt(c), cnt(p)
            obs.append({**{'c_' + k: v for k, v in a.items()}, **{'p_' + k: v for k, v in b.items()}})
        null = []
        for rep in range(100):
            c = rng.choice(cs, n, replace=False); p = rng.choice(ps, n, replace=False)
            pool = rng.permutation(np.r_[c, p]); a, b = cnt(pool[:n]), cnt(pool[n:])
            null.append({k: a[k] - b[k] for k in a})
        O = pd.DataFrame(obs); N = pd.DataFrame(null)
        r = dict(blood=blood_name, patients=pg, n_per_group=n)
        for k in ['n_q05', 'n_q05_r35', 'mean_abs_r']:
            r['ctrl_' + k] = O['c_' + k].mean(); r['pat_' + k] = O['p_' + k].mean()
            d = (O['c_' + k] - O['p_' + k]).mean(); r['diff_' + k] = d
            r['p_' + k] = (1 + (np.abs(N[k]) >= abs(d)).sum()) / (1 + len(N))
        rows.append(r); print(pd.Series(r).to_string(), flush=True)
pd.DataFrame(rows).to_csv(OUT + 'network_counts_matched.csv', index=False)
