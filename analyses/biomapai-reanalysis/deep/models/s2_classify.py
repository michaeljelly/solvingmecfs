"""Step 2: per-layer classifiers under person-grouped 5-fold CV (3 repeats) with nested tuning (inner 3-fold
StratifiedGroupKFold by person, AUC). Saves held-out predictions plus inner out-of-fold predictions for the
outer-training set (used for stacking in s2_fusion.py).
usage: python s2_classify.py TASK LAYER[,LAYER] MODEL[,MODEL]   TASK in all|short|long"""
import sys, time, json; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from prep import *
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, cross_val_predict
from sklearn.base import clone

REPS = 3
os.makedirs(OUT + 'results/s2', exist_ok=True)


def task_samples(task):
    m = meta
    if task == 'short': m = meta[meta.illness_duration.isin(['Short', 'Control'])]
    if task == 'long': m = meta[meta.illness_duration.isin(['Long', 'Control'])]
    return m


def task_folds(task, rep):
    m = task_samples(task)
    return person_folds(m.person.values, m.y.values, 5, seed=rep)


if __name__ == '__main__':
    TASK = sys.argv[1]; LAYS = sys.argv[2].split(','); MODELS = sys.argv[3].split(',')
    m = task_samples(TASK)
    for L in LAYS:
        Xall = load_layer(L)
        idx = Xall.index.intersection(m.index)
        X = Xall.loc[idx].values; y = m.loc[idx, 'y'].values; g = m.loc[idx, 'person'].values
        for M in MODELS:
            t0 = time.time(); rows = []
            for rep in range(REPS):
                pf = task_folds(TASK, rep); fold = np.array([pf[p] for p in g])
                for k in range(5):
                    tr = np.where(fold != k)[0]; te = np.where(fold == k)[0]
                    inner = StratifiedGroupKFold(3, shuffle=True, random_state=rep * 10 + k)
                    gs = GridSearchCV(make_pipe(L, M), grid(L, M), cv=inner, scoring='roc_auc', n_jobs=1, refit=True)
                    gs.fit(X[tr], y[tr], groups=g[tr])
                    pte = gs.predict_proba(X[te])[:, 1]
                    ptr = cross_val_predict(clone(gs.best_estimator_), X[tr], y[tr], groups=g[tr], cv=inner,
                                            method='predict_proba')[:, 1]
                    bp = json.dumps(gs.best_params_)
                    for i, p in zip(te, pte): rows.append([idx[i], g[i], y[i], rep, k, 'test', p, bp])
                    for i, p in zip(tr, ptr): rows.append([idx[i], g[i], y[i], rep, k, 'trainoof', p, bp])
            d = pd.DataFrame(rows, columns=['sample', 'person', 'y', 'rep', 'fold', 'set', 'pred', 'params'])
            d.to_csv(OUT + f'results/s2/pred_{TASK}_{L}_{M}.csv', index=False)
            te = d[d.set == 'test']
            aucs = [roc_auc_score(x.y, x.pred) for _, x in te.groupby('rep')]
            print(f'{TASK} {L} {M} n={len(idx)} AUC {np.mean(aucs):.3f}+-{np.std(aucs):.3f}  {time.time()-t0:.0f}s', flush=True)
