"""Minimal blood panel from replicated metabolites: train in one cohort (family), test in another (people never shared).
Features = pre-specified replicated items (no caffeine/drug/diet markers) available in both cohorts; values z-scored within cohort
(label-free); L2 logistic regression with fixed C=0.3; AUC with 2000-person bootstrap 95% CI.
Also: 5x5 person-level CV AUC inside BioMapAI and Germain 2022. Output: classifier_transfer.csv"""
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from mwlib import *

S, M = load()
PANEL = ['phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'indolelactate', 'imidazole lactate', 'alpha-hydroxyisovalerate', 'indolepropionate (IPA)',
         'trans-4-hydroxyproline', '5-(galactosylhydroxy)-lysine', '4-hydroxyglutamate', 'methionine sulfoxide / methionine', 'PC(14:0/16:0)',
         'plasmalogen PC(P-16:0/18:2)', 'glucuronate', '5-oxoproline / glutamate', 'carnitine', 'tiglylcarnitine (C5:1-DC)', 'kynurenine', 'TMAO / carnitine', 'dihydroceramide (d18:0/16:0)']
CORE = ['phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'indolelactate', 'imidazole lactate', 'trans-4-hydroxyproline', 'methionine sulfoxide / methionine', 'indolepropionate (IPA)']
D = {}
for st in S:
    Y = build_items(st, S, M)
    Y = Y[[c for c in PANEL if c in Y]]
    D[st] = ((Y - Y.mean()) / Y.std(), S[st]['C'].case.loc[Y.index])
rng = np.random.default_rng(0)


def boot_auc(y, p, B=2000):
    y, p = np.asarray(y), np.asarray(p); a = []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() == y[i].max(): continue
        a.append(roc_auc_score(y[i], p[i]))
    return np.percentile(a, [2.5, 97.5])


rows = []
def run(train, test, feats_pool, tag):
    Xs = [D[s][0] for s in train]; ys = [D[s][1] for s in train]
    feats = [f for f in feats_pool if all(f in X for X in Xs) and f in D[test][0]]
    if len(feats) < 2: return
    X = pd.concat([x[feats] for x in Xs]); y = pd.concat(ys)
    m = LogisticRegression(C=0.3, max_iter=2000).fit(X.values, y.values)
    Xt, yt = D[test][0][feats], D[test][1]
    p = m.predict_proba(Xt.values)[:, 1]; auc = roc_auc_score(yt, p); lo, hi = boot_auc(yt, p)
    rows.append(dict(panel=tag, train='+'.join(train), test=test, n_test=len(yt), n_feats=len(feats), AUC=auc, lo=lo, hi=hi,
                     features='; '.join(feats), coefs='; '.join(f'{f}:{c:+.2f}' for f, c in zip(feats, m.coef_[0]))))

MET = ['BioMapAI', 'Germain2022', 'Germain2018', 'Germain2020']
OTH = ['Che2022', 'NagySzakal2018', 'Naviaux2016']
for tag, pool in [('full panel', PANEL), ('core 7', CORE), ('aromatic lactates only', PANEL[:4])]:
    for te in MET[1:] + OTH: run(['BioMapAI'], te, pool, tag)
    for te in OTH: run(['BioMapAI', 'Germain2022'], te, pool, tag)
    for tr in OTH:
        for te in ['BioMapAI', 'Germain2022']: run([tr], te, pool, tag)
    run(['Germain2022'], 'BioMapAI', pool, tag)
R = pd.DataFrame(rows)
# within-cohort CV
for st in ['BioMapAI', 'Germain2022']:
    for tag, pool in [('full panel', PANEL), ('core 7', CORE)]:
        X, y = D[st][0][[f for f in pool if f in D[st][0]]], D[st][1]
        ps = np.zeros(len(y)); aucs = []
        for tr, te in RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=1).split(X, y):
            m = LogisticRegression(C=0.3, max_iter=2000).fit(X.values[tr], y.values[tr]); aucs.append(roc_auc_score(y.values[te], m.predict_proba(X.values[te])[:, 1]))
        R = pd.concat([R, pd.DataFrame([dict(panel=tag, train=st + ' (5x5 CV)', test=st, n_test=len(y), n_feats=X.shape[1], AUC=np.mean(aucs), lo=np.percentile(aucs, 2.5), hi=np.percentile(aucs, 97.5), features='; '.join(X.columns))])])
R.to_csv(OUT + 'classifier_transfer.csv', index=False)
pd.set_option('display.width', 200)
print(R[['panel', 'train', 'test', 'n_test', 'n_feats', 'AUC', 'lo', 'hi']].round(3).to_string())
print(R[(R.panel == 'full panel') & (R.train == 'BioMapAI')].coefs.iloc[0])
