"""Shared helpers for the models/ analysis (person-grouped evaluation of BioMapAI data)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1'); os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

B = '/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT = '/mnt/project-files/biomap/deep/models/'
LAYERS = ['immune', 'specie', 'kegg', 'metabolome', 'quest']

meta = pd.read_csv(B + 'metadata.csv', index_col=0)
meta['person'] = meta['sample_id_tp1']
meta['y'] = (meta.study_ptorhc == 'MECFS').astype(int)
meta['male'] = (meta.gender == 'Male').astype(int)
score = pd.read_csv(B + 'score.csv', index_col=0)
SCORES = list(score.columns)


def load_raw(layer):
    d = pd.read_csv(B + layer + '.csv', index_col=0).T
    return d.loc[d.index.intersection(meta.index)].astype(float)


def clr(X, pseudo=None):
    X = np.asarray(X, float)
    if pseudo is None:
        nz = X[X > 0]
        pseudo = nz.min() / 2 if nz.size else 1e-6
    L = np.log(X + pseudo)
    return L - L.mean(1, keepdims=True)


def person_folds(persons, labels, k=5, seed=0):
    """Assign each person to one of k folds, stratified by the person's label. Returns dict person->fold."""
    rng = np.random.RandomState(seed)
    df = pd.DataFrame({'p': persons, 'y': labels}).drop_duplicates('p')
    f = {}
    for yv, g in df.groupby('y'):
        ps = np.array(g.p.values, dtype=object); rng.shuffle(ps)
        for i, p in enumerate(ps):
            f[p] = i % k
    return f


def boot_auc(y, p, groups, n=1000, seed=0):
    """AUC with percentile bootstrap CI resampling persons (groups)."""
    y = np.asarray(y); p = np.asarray(p); groups = np.asarray(groups)
    auc = roc_auc_score(y, p)
    ug = np.unique(groups); idx = {g: np.where(groups == g)[0] for g in ug}
    rng = np.random.RandomState(seed); bs = []
    for _ in range(n):
        s = rng.choice(ug, len(ug), replace=True)
        ii = np.concatenate([idx[g] for g in s])
        if len(np.unique(y[ii])) < 2: continue
        bs.append(roc_auc_score(y[ii], p[ii]))
    return auc, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def person_level(df, pcol='pred'):
    """Average visit predictions per person. df has columns person, y, pcol."""
    return df.groupby('person').agg(y=('y', 'first'), pred=(pcol, 'mean')).reset_index()


def fmt(a, lo, hi):
    return f"{a:.3f} [{lo:.3f}-{hi:.3f}]"


def boot_auc_reps(d, pcol='pred', n=1000, seed=0, person_level_agg=False):
    """d: rows = held-out predictions with columns rep, person, y, pcol (one row per sample per repeat).
    Returns mean-over-repeats AUC and a 95% CI from bootstrapping persons (same persons resampled in every repeat).
    If person_level_agg, predictions are first averaged per person within each repeat."""
    d = d[['rep', 'person', 'y', pcol]].copy()
    if person_level_agg:
        d = d.groupby(['rep', 'person']).agg(y=('y', 'first'), p=(pcol, 'mean')).reset_index().rename(columns={'p': pcol})
    reps = sorted(d.rep.unique())
    per = [d[d.rep == r] for r in reps]
    est = np.mean([roc_auc_score(x.y, x[pcol]) for x in per])
    persons = d.person.unique(); rng = np.random.RandomState(seed)
    pidx = {p: i for i, p in enumerate(persons)}
    codes = [np.array([pidx[p] for p in x.person.values]) for x in per]
    bs = []
    for _ in range(n):
        cnt = np.bincount(rng.randint(0, len(persons), len(persons)), minlength=len(persons))
        vals = []
        for x, c in zip(per, codes):
            w = cnt[c]; m = w > 0; yy = x.y.values[m]
            if len(np.unique(yy)) == 2: vals.append(roc_auc_score(yy, x[pcol].values[m], sample_weight=w[m]))
        if vals: bs.append(np.mean(vals))
    return est, np.percentile(bs, 2.5), np.percentile(bs, 97.5)
