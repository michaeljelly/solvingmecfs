"""Case-control with all visits (random-intercept LMM), sex-stratified and sex-interaction models,
plus tp1-only Mann-Whitney and tp1-only adjusted OLS for comparison."""
import numpy as np, pandas as pd
from scipy import stats
from data import *
from lmmfast import fit_many, bh

m = load_meta()
rows = []; counts = []


GUARD = {}


def add(L, feats, model, contrast, r, j, ann, extra=None):
    bad = GUARD.get(model)
    if bad is not None:
        for k in ('est', 'se', 'p'): r[k][bad, :] = np.nan
    d = pd.DataFrame({'layer': L, 'feature': feats, 'model': model, 'contrast': contrast,
                      'est': r['est'][:, j], 'se': r['se'][:, j], 'p': r['p'][:, j],
                      'icc': r.get('icc', np.nan), 'annotation': ann})
    rows.append(d)


for L in LAYERS:
    X = load_layer(L)
    idx = X.index.intersection(m.index); M = m.loc[idx]; X = X.loc[idx]
    X = (X - X.mean()) / X.std()
    feats = list(X.columns); ann = annotate(L, feats)
    counts.append((L, len(M), M.person.nunique(), M[M.tp == 1].groupby('illness_duration').size().to_dict(),
                   M.drop_duplicates('person').groupby(['illness_duration', 'gender']).size().to_dict()))
    cov = lambda D, sex=True: [v for v in [D.age_c.values] + ([D.male.values] if sex else []) + [D.bmi_c.values, D.bmi_missing.values, D.tp2.values, D.tp3.values] if np.std(v) > 0]
    one = lambda D: np.ones(len(D))
    # A: patient vs control, all visits
    C = np.column_stack([one(M), M.case] + cov(M))
    add(L, feats, 'LMM_all', 'patient_vs_control', fit_many(X.values, C, M.person.values, L=np.eye(C.shape[1])[[1]]), 0, ann)
    # B: short and long vs control (+ short vs long)
    C = np.column_stack([one(M), M.short, M.long] + cov(M))
    Lc = np.zeros((3, C.shape[1])); Lc[0, 1] = 1; Lc[1, 2] = 1; Lc[2, 1] = 1; Lc[2, 2] = -1
    r = fit_many(X.values, C, M.person.values, L=Lc)
    for j, cn in enumerate(['short_vs_control', 'long_vs_control', 'short_vs_long']): add(L, feats, 'LMM_all', cn, r, j, ann)
    # sensitivity: complete-case BMI
    Mc = M[M.bmi_missing == 0]
    C = np.column_stack([one(Mc), Mc.short, Mc.long, Mc.age_c, Mc.male, Mc.bmi_c, Mc.tp2, Mc.tp3])
    r = fit_many(X.loc[Mc.index].values, C, Mc.person.values, L=np.eye(C.shape[1])[[1, 2]])
    add(L, feats, 'LMM_all_completeBMI', 'short_vs_control', r, 0, ann); add(L, feats, 'LMM_all_completeBMI', 'long_vs_control', r, 1, ann)
    # C/D: sex-stratified
    for sx, nm in [(1, 'LMM_males'), (0, 'LMM_females')]:
        Ms = M[M.male == sx]; Xs = X.loc[Ms.index]
        # guard: features (near-)constant within this sex, e.g. a drug taken by nobody of that sex
        GUARD[nm] = ((Xs.std() < 0.1) | (Xs.apply(lambda v: v.value_counts().iloc[0]) / len(Xs) > 0.9)).values
        C = np.column_stack([one(Ms), Ms.case] + cov(Ms, sex=False))
        add(L, feats, nm, 'patient_vs_control', fit_many(Xs.values, C, Ms.person.values, L=np.eye(C.shape[1])[[1]]), 0, ann)
        C = np.column_stack([one(Ms), Ms.short, Ms.long] + cov(Ms, sex=False))
        r = fit_many(Xs.values, C, Ms.person.values, L=np.eye(C.shape[1])[[1, 2]])
        add(L, feats, nm, 'short_vs_control', r, 0, ann); add(L, feats, nm, 'long_vs_control', r, 1, ann)
    # E: sex interaction (positive = case effect larger/more positive in males)
    C = np.column_stack([one(M), M.case, M.case * M.male] + cov(M))
    add(L, feats, 'LMM_sex_interaction', 'patient_x_male', fit_many(X.values, C, M.person.values, L=np.eye(C.shape[1])[[2]]), 0, ann)
    C = np.column_stack([one(M), M.short, M.long, M.short * M.male, M.long * M.male] + cov(M))
    r = fit_many(X.values, C, M.person.values, L=np.eye(C.shape[1])[[3, 4]])
    add(L, feats, 'LMM_sex_interaction', 'short_x_male', r, 0, ann); add(L, feats, 'LMM_sex_interaction', 'long_x_male', r, 1, ann)
    # tp1-only: adjusted OLS (LMM with a single visit == OLS; grid restricted to ICC 0) and Mann-Whitney
    T = M[M.tp == 1]; XT = X.loc[T.index]
    C = np.column_stack([one(T), T.short, T.long, T.age_c, T.male, T.bmi_c, T.bmi_missing])
    r = fit_many(XT.values, C, T.person.values, L=np.eye(C.shape[1])[[1, 2]], grid=[0.0])
    add(L, feats, 'OLS_tp1', 'short_vs_control', r, 0, ann); add(L, feats, 'OLS_tp1', 'long_vs_control', r, 1, ann)
    C = np.column_stack([one(T), T.case, T.age_c, T.male, T.bmi_c, T.bmi_missing])
    add(L, feats, 'OLS_tp1', 'patient_vs_control', fit_many(XT.values, C, T.person.values, L=np.eye(C.shape[1])[[1]], grid=[0.0]), 0, ann)
    for cn, grp in [('patient_vs_control', T.case == 1), ('short_vs_control', T.short == 1), ('long_vs_control', T.long == 1)]:
        A = XT[grp.values]; Bc = XT[(T.case == 0).values]
        u, pv = stats.mannwhitneyu(A.values, Bc.values, axis=0)
        cliff = 2 * u / (len(A) * len(Bc)) - 1
        rows.append(pd.DataFrame({'layer': L, 'feature': feats, 'model': 'MWU_tp1', 'contrast': cn, 'est': cliff, 'se': np.nan,
                                  'p': pv, 'icc': np.nan, 'annotation': ann}))
    print(L, X.shape, 'done', flush=True)

R = pd.concat(rows, ignore_index=True)
R['q'] = R.groupby(['layer', 'model', 'contrast']).p.transform(lambda p: pd.Series(bh(p.values), index=p.index))
R['ci_lo'] = R.est - 1.96 * R.se; R['ci_hi'] = R.est + 1.96 * R.se
R.to_csv(OUT + 'casecontrol_all_models.csv.gz', index=False)
pd.DataFrame(counts, columns=['layer', 'n_visits', 'n_persons', 'tp1_by_group', 'persons_by_group_sex']).to_csv(OUT + 'casecontrol_sample_counts.csv', index=False)
tab = R.assign(sig=R.q < 0.05).pivot_table(index=['model', 'contrast'], columns='layer', values='sig', aggfunc='sum')
tab.to_csv(OUT + 'casecontrol_n_significant.csv')
print(tab.to_string())
