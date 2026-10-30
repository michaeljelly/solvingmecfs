"""Step 5: permutation importance on held-out person-grouped folds (2 repeats x 5 folds), for chosen task/layer/model.
Hyper-parameters: the most frequently selected setting in the nested runs of s2 (no re-tuning on test folds).
Stability: fraction of the 10 held-out folds in which permuting the feature lowers AUC.
usage: python s5_importance.py TASK LAYER MODEL"""
import sys, json; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from prep import *
from s2_classify import task_samples, task_folds
from sklearn.inspection import permutation_importance

TASK, L, M = sys.argv[1:4]
m = task_samples(TASK); Xall = load_layer(L); idx = Xall.index.intersection(m.index)
X = Xall.loc[idx]; y = m.loc[idx, 'y'].values; g = m.loc[idx, 'person'].values
prev = pd.read_csv(OUT + f'results/s2/pred_{TASK}_{L}_{M}.csv')
params = json.loads(prev.drop_duplicates(['rep', 'fold']).params.value_counts().index[0])
print('params', params, flush=True)
imps = []
for rep in range(2):
    pf = task_folds(TASK, rep); fold = np.array([pf[p] for p in g])
    for k in range(5):
        tr = fold != k; te = fold == k
        pipe = make_pipe(L, M).set_params(**params).fit(X.values[tr], y[tr])
        r = permutation_importance(pipe, X.values[te], y[te], scoring='roc_auc', n_repeats=10, random_state=rep * 10 + k, n_jobs=1)
        imps.append(r.importances_mean)
    print('rep', rep, 'done', flush=True)
I = np.array(imps)
res = pd.DataFrame({'feature': X.columns, 'mean_drop_auc': I.mean(0), 'sd': I.std(0), 'frac_folds_positive': (I > 0).mean(0),
                    'median_rank': pd.DataFrame(-I).rank(axis=1).median(0).values})
# direction of association in the whole data (Spearman with label, on CLR scale for microbiome), for interpretation only
from scipy.stats import spearmanr
Xd = X.values
if L in COMP: Xd = CompFilter('clr', 0.0, COMP[L]).fit(Xd).transform(Xd)
res['rho_with_MECFS'] = [spearmanr(Xd[:, j], y)[0] if np.std(Xd[:, j]) > 0 else np.nan for j in range(Xd.shape[1])]
res = res.sort_values('mean_drop_auc', ascending=False)
res.to_csv(OUT + f'results/s5_importance_{TASK}_{L}_{M}.csv', index=False)
print(res.head(20).round(4).to_string())
