"""Step 4: symptom-score prediction within patients only (the paper's core 'symptom-specific biomarker' claim).
(a) Between-person: person-grouped 5-fold CV (2 repeats), elastic net (inner-tuned) and LightGBM (fixed, shallow),
    per layer and fused (mean of available layer predictions), plus a demographics+illness-duration baseline.
    Metrics: pooled out-of-fold R^2 and Spearman rho (95% CI by bootstrapping persons).
(b) Within-person: person-mean-centred features vs person-mean-centred scores (patients with >=2 visits),
    same CV. Tests whether visit-to-visit changes in a layer track visit-to-visit changes in symptoms."""
import sys; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from prep import *
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.metrics import r2_score
from scipy.stats import spearmanr

MODE = sys.argv[1] if len(sys.argv) > 1 else 'between'
REPS = 2
os.makedirs(OUT + 'results/s4', exist_ok=True)
pat = meta[(meta.y == 1) & meta.index.isin(score.index)].copy()
pat['short'] = (pat.illness_duration == 'Short').astype(int)


def boot_rho(y, p, g, n=500, seed=0):
    rng = np.random.RandomState(seed); ug = np.unique(g); idx = {k: np.where(g == k)[0] for k in ug}; bs = []
    for _ in range(n):
        ii = np.concatenate([idx[k] for k in rng.choice(ug, len(ug))]); bs.append(spearmanr(y[ii], p[ii])[0])
    return np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)


def centre(df, groups):
    return df - df.groupby(groups).transform('mean')


def layer_matrix(L):
    if L == 'demo':
        return pat[['age', 'male', 'bmi', 'short']].copy()
    X = load_layer(L); X = X.loc[X.index.intersection(pat.index)]
    if L in COMP:  # stateless CLR here so that centring (within-person mode) is on the transformed scale
        cf = CompFilter('clr', 0.0, COMP[L]).fit(X.values)
        X = pd.DataFrame(cf.transform(X.values), index=X.index)
        X = X.loc[:, (load_layer(L).loc[X.index] > 0).mean(0).values >= 0.1]
    return X


def pipes(L, M):
    st = [('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler())]
    if M == 'enet':
        if L == 'kegg': st.append(('sel', SelectKBest(f_regression, k=300)))
        st.append(('m', ElasticNet(l1_ratio=0.5, max_iter=5000)))
        return Pipeline(st), {'m__alpha': [0.03, 0.1, 0.3, 1.0]}
    st.append(('m', lgb.LGBMRegressor(n_estimators=200, learning_rate=0.03, num_leaves=4, min_child_samples=15,
                                      subsample=0.8, subsample_freq=1, colsample_bytree=0.3, reg_lambda=1.0,
                                      n_jobs=1, verbose=-1, random_state=0)))
    return Pipeline(st), None


if __name__ == '__main__':
    rows = []; allpred = []
    layers = ['demo', 'immune', 'specie', 'kegg', 'metabolome', 'quest'] if MODE == 'between' else ['immune', 'specie', 'kegg', 'metabolome', 'quest']
    SUF = ''
    if len(sys.argv) > 2: layers = sys.argv[2].split(','); SUF = '_' + '_'.join(layers)
    for L in layers:
        X = layer_matrix(L); S = score.loc[X.index]; g = pat.loc[X.index, 'person'].values
        if MODE == 'within':
            keep = pd.Series(g).map(pd.Series(g).value_counts()).values >= 2
            X = X[keep]; S = S[keep]; g = g[keep]
            X = centre(X.fillna(X.median()), g); S = centre(S, g)
        for M in ['enet', 'lgbm']:
            if L == 'demo' and M == 'lgbm': continue
            for sc in SCORES:
                y = S[sc].values; P = np.zeros((REPS, len(y)))
                for rep in range(REPS):
                    pf = person_folds(g, np.zeros(len(g)), 5, seed=rep); fold = np.array([pf[p] for p in g])
                    for k in range(5):
                        tr = fold != k; te = fold == k
                        pipe, gr = pipes(L, M)
                        if gr:
                            mdl = GridSearchCV(pipe, gr, cv=GroupKFold(3), scoring='r2', n_jobs=1).fit(X.values[tr], y[tr], groups=g[tr])
                        else:
                            mdl = pipe.fit(X.values[tr], y[tr])
                        P[rep, te] = mdl.predict(X.values[te])
                p = P.mean(0)
                r2s = [r2_score(y, P[r]) for r in range(REPS)]
                rho = spearmanr(y, p)[0]; lo, hi = boot_rho(y, p, g)
                rows.append(dict(mode=MODE, layer=L, model=M, score=sc, n=len(y), n_persons=len(np.unique(g)),
                                 r2=np.mean(r2s), rho=rho, rho_lo=lo, rho_hi=hi))
                allpred.append(pd.DataFrame({'sample': X.index, 'person': g, 'layer': L, 'model': M, 'score': sc, 'y': y, 'pred': p}))
            sub = pd.DataFrame(rows); sub = sub[(sub.layer == L) & (sub.model == M)]
            print(MODE, L, M, 'median R2 %.3f  max R2 %.3f (%s)  median rho %.2f' % (sub.r2.median(), sub.r2.max(),
                  sub.loc[sub.r2.idxmax(), 'score'], sub.rho.median()), flush=True)
            pd.DataFrame(rows).to_csv(OUT + f'results/s4_{MODE}{SUF}.csv', index=False)
    ap = pd.concat(allpred); ap.to_csv(OUT + f'results/s4/pred_{MODE}{SUF}.csv', index=False)
    # fusion: mean of available omics-layer predictions (demo excluded)
    for M in ['enet', 'lgbm']:
        for sc in SCORES:
            a = ap[(ap.model == M) & (ap.score == sc) & (ap.layer != 'demo')]
            f = a.groupby('sample').agg(person=('person', 'first'), y=('y', 'first'), pred=('pred', 'mean'))
            lo, hi = boot_rho(f.y.values, f.pred.values, f.person.values)
            rows.append(dict(mode=MODE, layer='fused', model=M, score=sc, n=len(f), n_persons=f.person.nunique(),
                             r2=r2_score(f.y, f.pred), rho=spearmanr(f.y, f.pred)[0], rho_lo=lo, rho_hi=hi))
    pd.DataFrame(rows).to_csv(OUT + f'results/s4_{MODE}{SUF}.csv', index=False)
    print('done')
