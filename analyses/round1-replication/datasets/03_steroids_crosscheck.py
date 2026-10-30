"""Cross-check of the steroid signal (androgenic / pregnenolone steroids, DHEA-S, cortisol) that came out of
Germain 2020 + Germain 2022, in cohorts from other labs: BioMapAI (Jackson Lab, Metabolon, LMM with person random
intercept, age/sex/BMI/visit adjusted) and Che 2022 (Columbia, testosterone + androsterone glucuronide only).
Also covariate-adjusted Germain 2022 (sex, age bin, BMI bin) and men-only estimates."""
import sys, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
sys.path.insert(0, '/mnt/project-files/biomap/deep/stats_replication/')
from data import load_meta, load_layer, B
from statsmodels.stats.multitest import multipletests

OUT = '/mnt/project-files/mecfs-research/datasets/results/'
IC = B + 'data/indenpendent_cohorts/'
rows = []

# ---------------- BioMapAI
m = load_meta(); X = load_layer('metabolome')
q = pd.read_csv(B + 'data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv'); q = q[q.TYPE == 'NAMED']
sub = q.drop_duplicates('CHEMICAL_NAME').set_index('CHEMICAL_NAME').SUB_PATHWAY.astype(str) if 'SUB_PATHWAY' in q else None
samp_map = None
idx = X.index.intersection(m.index)
if len(idx) == 0:
    raise SystemExit('BioMapAI metabolome sample ids do not match metadata index; inspect')
X = X.loc[idx]; M = m.loc[idx]
Z = (X - X.mean()) / X.std()
feats = {}
for grp in ['Androgenic Steroids', 'Pregnenolone Steroids', 'Corticosteroids', 'Progestin Steroids']:
    cols = [c for c in X.columns if sub.get(c) == grp]
    if len(cols) >= 2: feats[f'COMPOSITE {grp.lower()}'] = Z[cols].mean(axis=1); print('BioMapAI', grp, len(cols))
for c in ['dehydroepiandrosterone sulfate (DHEA-S)', 'cortisol', 'cortisone', 'androsterone sulfate', 'testosterone sulfate', 'pregnenolone sulfate', '11-ketoetiocholanolone glucuronide']:
    if c in Z: feats[c] = Z[c]
F = pd.DataFrame(feats)
for sex_label, keep in [('all', np.ones(len(M), bool)), ('men', (M.male == 1).values)]:
    for f in F.columns:
        d = pd.DataFrame({'y': F[f].values, 'case': M.case.values, 'age': M.age_c.values, 'male': M.male.values, 'bmi': M.bmi_c.values,
                          'bmim': M.bmi_missing.values, 'tp2': M.tp2.values, 'tp3': M.tp3.values, 'person': M.person.values})[keep]
        fml = 'y ~ case + age + bmi + bmim + tp2 + tp3' + (' + male' if sex_label == 'all' else '')
        try:
            r = smf.mixedlm(fml, d, groups=d.person).fit(reml=True)
            rows.append(dict(cohort='BioMapAI', subset=sex_label, feature=f, est=r.params['case'], se=r.bse['case'], p=r.pvalues['case'],
                             n_case=d[d.case == 1].person.nunique(), n_ctrl=d[d.case == 0].person.nunique()))
        except Exception as e:
            print('fail', f, e)

# ---------------- Germain 2022 D1-PRE, adjusted
d = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='ScaledImpDataZeroDrug&Tobacco', header=None)
fm = d.iloc[8:, :15].copy(); fm.columns = d.iloc[7, :15].values
gm = d.iloc[:7, 15:].T; gm.columns = gm.iloc[0]; gm = gm.iloc[1:]
gm.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'AgeBin', 'BMIBin', 'Bell']
G = d.iloc[8:, 16:].astype(float).T; G.columns = fm.BIOCHEMICAL.astype(str).values; G.index = gm.index
G = np.log2(G.loc[:, ~G.columns.duplicated()])
s22 = pd.Series(fm['SUB PATHWAY'].astype(str).values, index=fm.BIOCHEMICAL.astype(str).values); s22 = s22[~s22.index.duplicated()]
pre = (gm.Timepoint == 'D1-PRE').values; G = G[pre]; g = gm[pre]
Zg = (G - G.mean()) / G.std(); fg = {}
for grp in ['Androgenic Steroids', 'Pregnenolone Steroids', 'Corticosteroids', 'Progestin Steroids']:
    cols = [c for c in G.columns if s22.get(c) == grp]
    if len(cols) >= 2: fg[f'COMPOSITE {grp.lower()}'] = Zg[cols].mean(axis=1)
for c in ['dehydroepiandrosterone sulfate (DHEA-S)', 'cortisol', 'cortisone', 'androsterone sulfate', 'pregnenolone sulfate']:
    if c in Zg: fg[c] = Zg[c]
FG = pd.DataFrame(fg)
for sex_label in ['all', 'men']:
    keep = np.ones(len(g), bool) if sex_label == 'all' else (g.Sex == 'M').values
    for f in FG.columns:
        dd = pd.DataFrame({'y': FG[f].values, 'case': (g.Phenotype == 'CFS').astype(int).values, 'age': g.AgeBin.astype(float).values,
                           'male': (g.Sex == 'M').astype(float).values, 'bmi': g.BMIBin.astype(float).values})[keep]
        r = smf.ols('y ~ case + age + bmi' + (' + male' if sex_label == 'all' else ''), dd).fit()
        rows.append(dict(cohort='Germain2022_D1PRE_adj', subset=sex_label, feature=f, est=r.params['case'], se=r.bse['case'], p=r.pvalues['case'],
                         n_case=int(dd.case.sum()), n_ctrl=int((1 - dd.case).sum())))

# ---------------- Che 2022
che = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_metabolomics_yali_v2.csv', index_col=0)
ph = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_pheno_text_binned_v2.txt', sep='\t', index_col=0)
V = che.iloc[:, 3:].apply(pd.to_numeric, errors='coerce').T; V.columns = che.Molecule.astype(str).values
V = V.loc[:, ~V.columns.duplicated()]
agemap = {'18 - 24': 21, '25 - 34': 30, '35 - 44': 40, '45 - 54': 50, '55 - 64': 60, '65 and over': 70}
bmimap = {'underweight': 0, 'healthy weight': 1, 'overweight': 2, 'obese': 3}
ph = ph.loc[V.index]
for f in ['Testosterone', 'Androstan-3-ol-17-one 3-glucuronide']:
    v = V[f]; v = np.log2(v.where(v > 0).fillna(v[v > 0].min() / 2)); v = (v - v.mean()) / v.std()
    for sex_label in ['all', 'men']:
        dd = pd.DataFrame({'y': v.values, 'case': (ph.Phenotype == 'case').astype(int).values, 'age': ph.age_binned.map(agemap).values,
                           'male': (ph.sex == 'male').astype(float).values, 'bmi': ph.bmi_binned.map(bmimap).values, 'site': ph.site.astype(str).values}).dropna()
        if sex_label == 'men': dd = dd[dd.male == 1]
        r = smf.ols('y ~ case + age + bmi + C(site)' + (' + male' if sex_label == 'all' else ''), dd).fit()
        rows.append(dict(cohort='Che2022', subset=sex_label, feature=f, est=r.params['case'], se=r.bse['case'], p=r.pvalues['case'],
                         n_case=int(dd.case.sum()), n_ctrl=int((1 - dd.case).sum())))

R = pd.DataFrame(rows)
R.to_csv(OUT + 'steroids_crosscheck.csv', index=False)
pd.set_option('display.width', 200)
print(R.round(3).to_string())
