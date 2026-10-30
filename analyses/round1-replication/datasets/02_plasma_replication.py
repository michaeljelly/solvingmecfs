"""Plasma metabolomics replication in cohorts not used before:
  Germain 2020 Metabolites (Hanson lab, Metabolon + complex lipid panel, 26 ME/CFS / 26 controls, copy from GitHub/mapMECFS)
  Armstrong 2015 (MetaboLights MTBLS161, Australian serum NMR, 34 CFS / 25 controls)  -- only a few overlapping compounds
Plus a whole-profile reproducibility check against Germain 2022 (JCI Insight, D1-PRE baseline), and a model-transfer test
(train on Germain 2022 baseline, test on Germain 2020; person = sample, so every split is person-level).
Effect = SMD of log2 abundance (z-scored within cohort), OLS, no covariates (none released for Germain 2020 / MTBLS161)."""
import numpy as np, pandas as pd, statsmodels.api as sm, os
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

D = '/mnt/project-files/mecfs-research/datasets/'
OUT = D + 'results/'; os.makedirs(OUT, exist_ok=True)
IC = '/mnt/project-files/biomap/BioMapAI/data/indenpendent_cohorts/'


def smd(X, case, cov=None):
    X = (X - X.mean()) / X.std(); C = np.column_stack([np.ones(len(case)), case] + ([] if cov is None else [cov])); out = []
    for f in X.columns:
        y = X[f].values
        if not np.isfinite(y).all() or np.std(y) == 0: continue
        r = sm.OLS(y, C).fit(); out.append((f, r.params[1], r.bse[1], r.pvalues[1]))
    return pd.DataFrame(out, columns=['feature', 'est', 'se', 'p']).set_index('feature')


def log_impute(V, minprev=0.5):
    V = V.apply(pd.to_numeric, errors='coerce'); V = V.where(V > 0)
    V = V.loc[:, V.notna().mean() >= minprev]
    V = V.fillna(V.min() / 2); return np.log2(V)


# ------------------------------------------------ Germain 2020
x = pd.read_excel(D + 'repos/Nadjoubiii_Myalgic-Encephalomyelitis-project/datasets/4.study_data_mecfs_metabolomics.xlsx', sheet_name=None)
dd = x['data_dictionary'].set_index('compound_id'); dm = x['data_matrix'].set_index('compound_id')
sm20 = x['sample_metadata'].set_index('sample_id')
G20 = dm.T.loc[sm20.index]; G20.columns = dd.loc[G20.columns, 'BIOCHEMICAL'].astype(str).values
G20 = G20.loc[:, ~G20.columns.duplicated()]
G20 = log_impute(G20)
# LABEL CHECK: the GitHub copy marks internal IDs C1-C26 as 'case' and P1-P26 as 'control'. Prescription-drug traces say the
# reverse: the 'P' group carries metoprolol (9 vs 2), zolpidem (5 vs 1), duloxetine, venlafaxine, mirtazapine, bupropion,
# lamotrigine, topiramate, methylphenidate, H2 blockers, cetirizine; the 'C' group carries almost none. So P = patient, C = control.
# (see 02b_label_check output results/germain2020_drug_label_check.csv). We use internal ID prefix, not the health_status column.
case20 = sm20.sample_id_internal.astype(str).str.startswith('P').astype(int).values
sub20 = dd.drop_duplicates('BIOCHEMICAL').set_index('BIOCHEMICAL').SUB_PATHWAY.astype(str)
sup20 = dd.drop_duplicates('BIOCHEMICAL').set_index('BIOCHEMICAL').SUPER_PATHWAY.astype(str)
print('Germain2020', G20.shape, case20.sum(), (1 - case20).sum())

# ------------------------------------------------ Germain 2022 baseline (D1-PRE)
d = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='ScaledImpDataZeroDrug&Tobacco', header=None)
fm = d.iloc[8:, :15].copy(); fm.columns = d.iloc[7, :15].values
gm = d.iloc[:7, 15:].T; gm.columns = gm.iloc[0]; gm = gm.iloc[1:]
gm.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'AgeBin', 'BMIBin', 'Bell']
G = d.iloc[8:, 16:].astype(float).T; G.columns = fm.BIOCHEMICAL.astype(str).values; G.index = gm.index
G = G.loc[:, ~G.columns.duplicated()]; G = np.log2(G)
G = G.loc[:, (G > G.min() + 1e-9).mean() >= 0.5]
sub22 = pd.Series(fm['SUB PATHWAY'].astype(str).values, index=fm.BIOCHEMICAL.astype(str).values); sub22 = sub22[~sub22.index.duplicated()]
pre = (gm.Timepoint == 'D1-PRE').values
G22 = G[pre]; g22 = gm[pre]; case22 = (g22.Phenotype == 'CFS').astype(int).values
print('Germain2022 D1-PRE', G22.shape, case22.sum(), (1 - case22).sum())


def composites(X, sub):
    sub = sub.reindex(X.columns).fillna(''); z = (X - X.mean()) / X.std()
    g = {'COMPOSITE long-chain free fatty acids': sub.isin(['Long Chain Polyunsaturated Fatty Acid (n3 and n6)', 'Long Chain Saturated Fatty Acid', 'Long Chain Monounsaturated Fatty Acid', 'Long Chain Fatty Acid', 'Polyunsaturated Fatty Acid (n3 and n6)']),
         'COMPOSITE long-chain acylcarnitines': sub.str.contains('Acyl Carnitine', regex=False) & ~sub.str.contains('Short|Medium', regex=True),
         'COMPOSITE ketone bodies': sub.isin(['Ketone Bodies']),
         'COMPOSITE ceramides': sub.isin(['Ceramides']),
         'COMPOSITE sphingomyelins': sub.isin(['Sphingomyelins']),
         'COMPOSITE plasmalogens': sub.isin(['Plasmalogen']),
         'COMPOSITE diacylglycerols (Metabolon)': sub.isin(['Diacylglycerol']),
         'COMPOSITE primary bile acids': sub.isin(['Primary Bile Acid Metabolism']),
         'COMPOSITE secondary bile acids': sub.str.startswith('Secondary Bile Acid'),
         'COMPOSITE androgenic steroids': sub.isin(['Androgenic Steroids']),
         'COMPOSITE pregnenolone steroids': sub.isin(['Pregnenolone Steroids']),
         'COMPOSITE corticosteroids': sub.isin(['Corticosteroids']),
         'COMPOSITE tryptophan metabolism': sub.isin(['Tryptophan Metabolism'])}
    out = {k: z.loc[:, v.values].mean(axis=1) for k, v in g.items() if v.sum() >= 2}
    if 'kynurenine' in X and 'tryptophan' in X: out['kynurenine/tryptophan (log2 ratio)'] = X['kynurenine'] - X['tryptophan']
    return pd.DataFrame(out), {k: int(v.sum()) for k, v in g.items()}


c20, n20 = composites(G20, sub20); c22, n22 = composites(G22, sub22)
print('composite sizes G20', n20); print('composite sizes G22', n22)
# Germain 2020 complex-lipid platform: TAG / DAG / free fatty acid super-pathway sums
z20 = (G20 - G20.mean()) / G20.std()
for sp in ['Triacylglycerol', 'Diacylglycerol', 'Free Fatty Acids', 'Cholesterol Ester', 'Phosphatidylcholine', 'Phosphatidylethanolamine', 'Lysophosphatidylcholine']:
    cols = [c for c in G20.columns if sup20.get(c) == sp]
    if len(cols) >= 2: c20[f'COMPLEX-LIPID {sp} (mean z)'] = z20[cols].mean(axis=1)

TARGETS = ['indolepropionate', 'indolelactate', 'indoleacetate', '3-indoxyl sulfate', 'kynurenine', 'kynurenate', 'tryptophan', 'serotonin',
           'vanillylmandelate (VMA)', '3-hydroxybutyrate (BHBA)', 'acetoacetate', 'palmitate (16:0)', 'oleate/vaccenate (18:1)', 'linoleate (18:2n6)',
           'glycocholate', 'glycochenodeoxycholate', 'cholate', 'chenodeoxycholate', 'deoxycholate', 'glycodeoxycholate', 'lithocholate', 'ursodeoxycholate',
           'caffeine', 'paraxanthine', 'hippurate', 'catechol sulfate', 'choline', 'carnitine', 'lactate', 'pyruvate', 'citrate', 'succinate',
           'cortisol', 'cortisone', 'dehydroepiandrosterone sulfate (DHEA-S)', 'imidazole propionate', 'trimethylamine N-oxide', 'p-cresol sulfate', 'phenylacetylglutamine']
T20 = G20[[t for t in TARGETS if t in G20]].join(c20)
T22 = G22[[t for t in TARGETS if t in G22]].join(c22)
r20 = smd(T20, case20); cov22 = np.column_stack([g22.AgeBin.astype(float), (g22.Sex == 'M').astype(float), g22.BMIBin.astype(float)])
r22 = smd(T22, case22, cov22)  # Germain 2022 adjusted for age bin, sex, BMI bin
tab = r20.add_prefix('G2020_').join(r22.add_prefix('G2022_'), how='outer')
tab['G2020_q_targets'] = np.nan
m = tab.G2020_p.notna(); tab.loc[m, 'G2020_q_targets'] = multipletests(tab.loc[m, 'G2020_p'], method='fdr_bh')[1]

# ------------------------------------------------ prior per-cohort estimates (BioMapAI deep work) for meta-analysis
prior = pd.read_csv('/mnt/project-files/biomap/deep/stats_replication/replication_per_cohort.csv')
prior = prior[prior.cohort.isin(['Che2022_Columbia_plasma'])].set_index('feature')
ext = pd.read_csv('/mnt/project-files/biomap/deep/coupling/external_replication.csv')
che_ipa = ext[(ext.cohort == 'Che2022') & ext.metabolite.str.startswith('indolepropionate [Indole')]
rt = pd.read_csv('/mnt/project-files/biomap/deep/stats_replication/replication_table.csv').set_index('feature')
bm = pd.DataFrame({'est': rt.BioMapAI_est, 'p': rt.BioMapAI_p}).dropna()
bm['se'] = bm.est.abs() / stats.norm.isf(bm.p / 2)
bm.loc['indolepropionate'] = [-0.3548, 0.0040, (-0.1132 + 0.5965) / 3.92]  # BioMapAI deep coupling/levels_mixed.csv (LMM, age/sex/BMI adj)
che = pd.DataFrame({'est': prior.est, 'se': prior.se, 'p': prior.p})
if len(che_ipa):
    r0 = che_ipa.iloc[0]; che.loc['indolepropionate'] = [r0.beta_SD, (r0.hi - r0.lo) / 3.92, r0.p]


def meta(es, ss):
    es, ss = np.array(es), np.array(ss); w = 1 / ss ** 2; fe = (w * es).sum() / w.sum(); k = len(es)
    Q = (w * (es - fe) ** 2).sum(); tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (ss ** 2 + tau2); re_ = (wr * es).sum() / wr.sum(); rse = np.sqrt(1 / wr.sum())
    return re_, rse, 2 * stats.norm.sf(abs(re_ / rse)), k


rows = []
for f in tab.index:
    E, S, src = [], [], []
    for name, df_e, df_s in [('G2020', tab.get('G2020_est'), tab.get('G2020_se')), ('G2022', tab.get('G2022_est'), tab.get('G2022_se'))]:
        if pd.notna(df_e.get(f)): E.append(df_e[f]); S.append(df_s[f]); src.append(name)
    if f in che.index: E.append(che.loc[f, 'est']); S.append(che.loc[f, 'se']); src.append('Che2022')
    if f in bm.index: E.append(bm.loc[f, 'est']); S.append(bm.loc[f, 'se']); src.append('BioMapAI')
    if len(E) >= 2:
        re_, rse, p, k = meta(E, S)
        rows.append(dict(feature=f, cohorts='+'.join(src), k=k, meta_RE=re_, meta_lo=re_ - 1.96 * rse, meta_hi=re_ + 1.96 * rse, meta_p=p,
                         n_neg=int((np.array(E) < 0).sum()), n_pos=int((np.array(E) > 0).sum())))
MT = pd.DataFrame(rows).set_index('feature')
tab = tab.join(MT)
tab.sort_values('G2020_p').to_csv(OUT + 'plasma_targets_G2020_G2022_meta.csv')
pd.set_option('display.width', 250)
print(tab[['G2020_est', 'G2020_p', 'G2020_q_targets', 'G2022_est', 'G2022_p', 'cohorts', 'meta_RE', 'meta_p', 'n_neg', 'n_pos']].sort_values('G2020_p').round(3).to_string())

# ------------------------------------------------ whole-profile reproducibility: all shared compounds
shared = sorted(set(G20.columns) & set(G22.columns))
a20 = smd(G20[shared], case20); a22 = smd(G22[shared], case22, cov22)
J = a20.add_prefix('G2020_').join(a22.add_prefix('G2022_'), how='inner')
J['G2020_q'] = multipletests(J.G2020_p, method='fdr_bh')[1]; J['G2022_q'] = multipletests(J.G2022_p, method='fdr_bh')[1]
J['sub'] = sub22.reindex(J.index)
J.sort_values('G2020_p').to_csv(OUT + 'plasma_allshared_G2020_vs_G2022.csv')
rho, p = stats.spearmanr(J.G2020_est, J.G2022_est)
print(f'\nshared compounds {len(J)}; effect-size correlation G2020 vs G2022 rho={rho:.2f} p={p:.2g}')
print('G2020 q<0.05:', (J.G2020_q < 0.05).sum(), ' G2022 q<0.05:', (J.G2022_q < 0.05).sum(), ' both q<0.1 same dir:',
      ((J.G2020_q < 0.1) & (J.G2022_q < 0.1) & (np.sign(J.G2020_est) == np.sign(J.G2022_est))).sum())
top22 = J[J.G2022_p < 0.01]
print(f'Germain2022 compounds p<0.01: {len(top22)}; same direction in G2020: {(np.sign(top22.G2020_est) == np.sign(top22.G2022_est)).sum()}, and G2020 p<0.05 same dir: {((np.sign(top22.G2020_est) == np.sign(top22.G2022_est)) & (top22.G2020_p < 0.05)).sum()}')
print(J[(J.G2020_p < 0.01) | (J.G2022_p < 0.005)].sort_values('G2020_p')[['sub', 'G2020_est', 'G2020_p', 'G2022_est', 'G2022_p']].round(3).head(40).to_string())

# ------------------------------------------------ model transfer: train G2022 baseline -> test G2020
Xtr = G22[shared]; Xte = G20[shared]
Xtr = (Xtr - Xtr.mean()) / Xtr.std(); Xte = (Xte - Xte.mean()) / Xte.std()  # within-cohort standardisation (platform/batch)
keep = Xtr.columns[Xtr.notna().all() & Xte.notna().all()]
res = []
for C in [0.01, 0.05, 0.1]:
    mdl = LogisticRegression(penalty='l1', C=C, solver='liblinear', max_iter=5000)
    cv = cross_val_predict(mdl, Xtr[keep].values, case22, cv=StratifiedKFold(5, shuffle=True, random_state=0), method='predict_proba')[:, 1]
    mdl.fit(Xtr[keep].values, case22); pr = mdl.predict_proba(Xte[keep].values)[:, 1]
    res.append((C, roc_auc_score(case22, cv), roc_auc_score(case20, pr), int((mdl.coef_ != 0).sum())))
    print(f'L1-logistic C={C}: internal CV AUC (G2022) {res[-1][1]:.2f}; external AUC on G2020 {res[-1][2]:.2f}; features {res[-1][3]}')
# reverse direction
for C in [0.05, 0.1]:
    mdl = LogisticRegression(penalty='l1', C=C, solver='liblinear', max_iter=5000).fit(Xte[keep].values, case20)
    print(f'reverse G2020->G2022 C={C}: external AUC {roc_auc_score(case22, mdl.predict_proba(Xtr[keep].values)[:, 1]):.2f}')
# pre-specified single-compound tests (no fitting): IPA lower -> AUC
for f in ['indolepropionate', '3-indoxyl sulfate']:
    for nm, X, cs in [('G2020', G20, case20), ('G2022', G22, case22)]:
        if f in X: print(f'{f} {nm}: AUC (lower=case) {roc_auc_score(cs, -X[f]):.2f}')
pd.DataFrame(res, columns=['C', 'cvAUC_G2022', 'extAUC_G2020', 'n_features']).to_csv(OUT + 'plasma_transfer_auc.csv', index=False)
