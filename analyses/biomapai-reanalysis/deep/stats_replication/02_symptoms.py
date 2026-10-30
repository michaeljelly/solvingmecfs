"""Symptom associations in patients (all visits): symptom ~ feature_between + feature_within + age + sex + visit + (1|person).
between = person mean of the feature, within = visit deviation from the person mean (person-mean centring)."""
import numpy as np, pandas as pd
from data import *
from lmmfast import fit_many, bh

m = load_meta()
P0 = m[(m.case == 1) & m.severity.notna()]


def composites(X, mm):
    sp = mm.SUB_PATHWAY.reindex(X.columns)
    z = (X - X.mean()) / X.std()
    grp = {
        'FFA_long_chain': sp.isin(['Long Chain Polyunsaturated Fatty Acid (n3 and n6)', 'Long Chain Saturated Fatty Acid', 'Long Chain Monounsaturated Fatty Acid']),
        'hydroxy_FA': sp == 'Fatty Acid, Monohydroxy',
        'acylcarnitines': sp.fillna('').str.startswith('Fatty Acid Metabolism (Acyl Carnitine'),
        'acylcarnitines_long_chain': sp.isin(['Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)']),
        'ketone_bodies': pd.Series(X.columns.isin(['3-hydroxybutyrate (BHBA)', 'acetoacetate']), index=X.columns),
        'ceramides': sp == 'Ceramides', 'sphingomyelins': sp == 'Sphingomyelins',
        'primary_bile_acids': sp == 'Primary Bile Acid Metabolism', 'secondary_bile_acids': sp.fillna('').str.startswith('Secondary Bile Acid'),
    }
    out = {k: z.loc[:, v.values].mean(axis=1) for k, v in grp.items() if v.sum() > 0}
    n = {k: int(v.sum()) for k, v in grp.items()}
    if 'kynurenine' in X: out['kynurenine'] = z['kynurenine']
    if 'kynurenine' in X and 'tryptophan' in X: out['kyn_trp_ratio_log2'] = X['kynurenine'] - X['tryptophan']
    if '3-hydroxybutyrate (BHBA)' in X: out['BHB'] = z['3-hydroxybutyrate (BHBA)']
    if 'vanillylmandelate (VMA)' in X: out['VMA'] = z['vanillylmandelate (VMA)']
    return pd.DataFrame(out), n


rows = []
for L in LAYERS + ['metab_composites']:
    if L == 'metab_composites':
        Xm = load_layer('metabolome'); X, ncomp = composites(Xm, metab_meta()); print('composite sizes', ncomp)
        rel = pd.read_csv(B + 'data/metagenomics/metaphlan4_species_full.csv', index_col=0).select_dtypes('number').T
        bp = [c for c in rel.columns if any(k in c for k in ['Faecalibacterium', 'Roseburia', 'Coprococcus', 'Anaerostipes', 'Agathobacter'])]
        bpc = np.log10(rel[bp].sum(axis=1) / rel.sum(axis=1) + 1e-4).rename('butyrate_producers_log10')
        X = X.join(bpc, how='outer')
        X = X.loc[:, :]
    else:
        X = load_layer(L)
    idx = X.index.intersection(P0.index); P = P0.loc[idx]; X = X.loc[idx]
    X = (X - X.mean()) / X.std()
    feats = list(X.columns)
    ann = annotate(L, feats) if L != 'metab_composites' else [''] * len(feats)
    for f in feats:  # composites may have NaN where one source missing
        pass
    ok_rows = X.notna().all(axis=1) if L != 'metab_composites' else None
    for s, sname in SYMPTOMS.items():
        if L == 'metab_composites':
            # fit each composite separately on its available rows
            for f in feats:
                d = P.join(X[[f]]).dropna(subset=[f, s])
                bt = d[f].groupby(d.person).transform('mean'); wi = d[f] - bt
                C = np.column_stack([np.ones(len(d)), d.age_c, d.male, d.tp2, d.tp3])
                r = fit_many(d[s].values, C, d.person.values, F=np.stack([bt.values, wi.values], 1)[None])
                rows.append((L, f, sname, len(d), d.person.nunique(), r['est'][0, 0], r['se'][0, 0], r['p'][0, 0],
                             r['est'][0, 1], r['se'][0, 1], r['p'][0, 1], ''))
            continue
        y = P[s].values
        bt = X.groupby(P.person.values).transform('mean'); wi = X - bt
        F = np.stack([bt.values.T, wi.values.T], axis=2)  # (K,n,2)
        C = np.column_stack([np.ones(len(P)), P.age_c, P.male, P.tp2, P.tp3])
        r = fit_many(y, C, P.person.values, F=F)
        for k, f in enumerate(feats):
            rows.append((L, f, sname, len(P), P.person.nunique(), r['est'][k, 0], r['se'][k, 0], r['p'][k, 0],
                         r['est'][k, 1], r['se'][k, 1], r['p'][k, 1], ann[k]))
    print(L, 'done', flush=True)

R = pd.DataFrame(rows, columns=['layer', 'feature', 'symptom', 'n_visits', 'n_persons', 'between_est', 'between_se', 'between_p',
                                'within_est', 'within_se', 'within_p', 'annotation'])
for comp in ['between', 'within']:
    R[comp + '_q'] = R.groupby(['layer', 'symptom'])[comp + '_p'].transform(lambda p: pd.Series(bh(p.values), index=p.index))
R.to_csv(OUT + 'symptom_mixed_models.csv.gz', index=False)
S = R[R.layer != 'metab_composites']
print(S.assign(b=S.between_q < 0.05, w=S.within_q < 0.05).groupby(['symptom', 'layer'])[['b', 'w']].sum().unstack().to_string())
pd.set_option('display.width', 250)
C = R[R.layer == 'metab_composites']
print(C[['feature', 'symptom', 'n_persons', 'between_est', 'between_p', 'within_est', 'within_p']].round(4).to_string())
C.to_csv(OUT + 'symptom_targeted_composites.csv', index=False)
