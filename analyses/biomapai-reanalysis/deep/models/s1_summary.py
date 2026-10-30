"""Summarise step 1 (authors' BioMapAI under person vs visit split) incl. its 12-score predictions within patients."""
import sys, glob; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from common import *
from sklearn.metrics import r2_score, balanced_accuracy_score
from scipy.stats import spearmanr

rows, srows = [], []
for f in sorted(glob.glob(OUT + 'results/s1/oof_*.csv')):
    L = f.split('oof_')[1][:-4]; d = pd.read_csv(f)
    for split, x in d.groupby('split'):
        reps = [roc_auc_score(z.y, z.prob) for _, z in x.groupby('rep')]
        acc = [((z.prob > .5) == z.y).mean() for _, z in x.groupby('rep')]
        bacc = [balanced_accuracy_score(z.y, z.prob > .5) for _, z in x.groupby('rep')]
        a, lo, hi = boot_auc_reps(x, pcol='prob')
        pa, plo, phi = boot_auc_reps(x, pcol='prob', person_level_agg=True)
        rows.append(dict(layer=L, split=split, n_reps=len(reps), auc_mean_reps=np.mean(reps), auc_sd=np.std(reps),
                         auc=a, lo=lo, hi=hi, auc_person=pa, p_lo=plo, p_hi=phi, acc=np.mean(acc), bal_acc=np.mean(bacc)))
        if split == 'person':  # score predictions within patients (their 'symptom' claim)
            sp = x[x.y == 1].groupby('sample')[SCORES].mean()
            truth = score.loc[sp.index]
            spa = x.groupby('sample')[SCORES].mean(); ta = score.loc[spa.index]
            for s in SCORES:
                srows.append(dict(layer=L, score=s, n=len(sp), r2=r2_score(truth[s], sp[s]),
                                  rho=spearmanr(truth[s], sp[s])[0], rho_all_subjects=spearmanr(ta[s], spa[s])[0]))
r = pd.DataFrame(rows); r.to_csv(OUT + 'results/s1_summary.csv', index=False)
s = pd.DataFrame(srows); s.to_csv(OUT + 'results/s1_dnn_scores_within_patients.csv', index=False)
print(r.round(3).to_string())
print(s.pivot(index='score', columns='layer', values='r2').round(2).to_string())
print(s.pivot(index='score', columns='layer', values='rho').round(2).to_string())
