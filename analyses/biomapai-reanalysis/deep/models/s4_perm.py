"""Permutation null for the strongest within-person symptom hit (visits shuffled within each person).
usage: python s4_perm.py LAYER MODEL SCORE NPERM"""
import sys; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from s4_symptoms import *

L, M, SC, NP = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
X = layer_matrix(L); S = score.loc[X.index]; g = pat.loc[X.index, 'person'].values
keep = pd.Series(g).map(pd.Series(g).value_counts()).values >= 2
X = X[keep]; S = S[keep]; g = g[keep]
X = centre(X.fillna(X.median()), g); y0 = centre(S, g)[SC].values


def cv_rho(y):
    P = np.zeros((REPS, len(y)))
    for rep in range(REPS):
        pf = person_folds(g, np.zeros(len(g)), 5, seed=rep); fold = np.array([pf[p] for p in g])
        for k in range(5):
            tr = fold != k; te = fold == k; pipe, gr = pipes(L, M)
            mdl = GridSearchCV(pipe, gr, cv=GroupKFold(3), scoring='r2').fit(X.values[tr], y[tr], groups=g[tr]) if gr else pipe.fit(X.values[tr], y[tr])
            P[rep, te] = mdl.predict(X.values[te])
    return spearmanr(y, P.mean(0))[0]


obs = cv_rho(y0); rng = np.random.RandomState(0); null = []
for i in range(NP):
    yp = y0.copy()
    for p in np.unique(g):
        ii = np.where(g == p)[0]; yp[ii] = yp[rng.permutation(ii)]
    null.append(cv_rho(yp))
null = np.array(null)
print(f'{L} {M} {SC}: observed rho {obs:.3f}; null mean {null.mean():.3f} sd {null.std():.3f}; perm p = {(1 + (null >= obs).sum()) / (NP + 1):.3f}')
pd.DataFrame({'null_rho': null}).assign(observed=obs).to_csv(OUT + f'results/s4_perm_{L}_{M}.csv', index=False)
