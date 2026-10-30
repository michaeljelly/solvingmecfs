"""Step 1: authors' BioMapAI model under person-grouped 5-fold CV vs their visit-level split.
Within each training fold patient visits are randomly down-sampled to the number of control visits (as the authors did
to balance classes). 3 repeats. Saves out-of-fold disease probabilities and predicted 12 scores."""
import sys; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from biomapai_dnn import *
from sklearn.model_selection import StratifiedKFold

LAY = sys.argv[1:] or ['specie', 'immune', 'metabolome', 'quest', 'kegg', 'omics']
REPS = int(os.environ.get('S1_REPS', 3)); VREPS = int(os.environ.get('S1_VREPS', 2))
os.makedirs(OUT + 'results/s1', exist_ok=True)
for L in LAY:
    X = load_raw(L); X = X.loc[X.index.intersection(score.index)]
    S = score.loc[X.index].values; y = meta.loc[X.index, 'y'].values; persons = meta.loc[X.index, 'person'].values
    rows = []
    for split in ['person', 'visit']:
        for rep in range(REPS if split == 'person' else VREPS):
            if split == 'person':
                pf = person_folds(persons, y, 5, seed=rep); fold = np.array([pf[p] for p in persons])
            else:
                fold = np.zeros(len(y), int)
                for k, (_, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=rep).split(X, y)): fold[te] = k
            rng = np.random.RandomState(100 + rep)
            for k in range(5):
                tr = np.where(fold != k)[0]; te = np.where(fold == k)[0]
                pos = tr[y[tr] == 1]; neg = tr[y[tr] == 0]
                tr_bal = np.concatenate([rng.choice(pos, len(neg), replace=False), neg])
                p, s = fit_predict(X.values[tr_bal], S[tr_bal], y[tr_bal], X.values[te], seed=rep * 10 + k)
                for j, i in enumerate(te):
                    rows.append([X.index[i], persons[i], y[i], split, rep, k, p[j]] + list(s[j]))
            d = pd.DataFrame(rows, columns=['sample', 'person', 'y', 'split', 'rep', 'fold', 'prob'] + SCORES)
            dd = d[(d.split == split) & (d.rep == rep)]
            print(L, split, rep, 'AUC %.3f acc %.3f' % (roc_auc_score(dd.y, dd.prob), ((dd.prob > .5) == dd.y).mean()), flush=True)
            d.to_csv(OUT + f'results/s1/oof_{L}.csv', index=False)
