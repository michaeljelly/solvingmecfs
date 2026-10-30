"""External replication of the case-control differences in microbe-derived plasma metabolites.
Germain 2022 (JCI Insight, Metabolon, baseline D1-PRE samples only): log2 scaled-imputed values, OLS with sex, age bin, BMI bin.
Che 2022 (Int J Mol Sci, r56): log2 peak areas (zeros -> half-min), OLS with sex, age bin, BMI bin, site, season."""
from cc import *
import statsmodels.formula.api as smf
CDIR = B + 'data/indenpendent_cohorts/'
KEY = ['indolepropionate', 'indolelactate', 'indoleacetate', 'hippurate', 'catechol sulfate', 'cinnamoylglycine',
       '3-phenylpropionate (hydrocinnamate)', 'p-cresol glucuronide*', 'p-cresol sulfate', 'glycodeoxycholate 3-sulfate',
       'taurodeoxycholate', 'deoxycholate', 'kynurenate', 'kynurenine', 'tryptophan', '3-indoxyl sulfate',
       'imidazole propionate', 'trimethylamine N-oxide', 'phenylacetylglutamine']
rows = []
o = pd.read_excel(CDIR + 'jci.insight.2022.xlsx', sheet_name='ScaledImpDataZeroDrug&Tobacco', header=None)
info = o.iloc[:7, 15:].T; info.columns = info.iloc[0]; info = info.iloc[1:]
info.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'agebin', 'bmibin', 'bell']
names = o.iloc[8:, 1].astype(str).values
vals = o.iloc[8:, 16:].astype(float).values
V = pd.DataFrame(vals.T, index=info.index, columns=names)
info.index = V.index
print('Germain phenotypes', info.Phenotype.value_counts().to_dict(), 'timepoints', info.Timepoint.value_counts().to_dict())
sel = info.Timepoint == 'D1-PRE'
G = info[sel].copy(); G['y'] = (G.Phenotype != 'Control').astype(int)
for k in KEY:
    if k not in V.columns: continue
    x = np.log2(V.loc[G.index, k].astype(float))
    G['v'] = (x - x.mean()) / x.std()
    m = smf.ols('v ~ y + C(Sex) + C(agebin) + C(bmibin)', G).fit()
    ci = m.conf_int().loc['y']
    rows.append(dict(cohort='Germain2022_D1PRE', metabolite=k, n_case=int(G.y.sum()), n_ctrl=int((1 - G.y).sum()),
                     beta_SD=m.params['y'], lo=ci[0], hi=ci[1], p=m.pvalues['y']))
# Che 2022
d = pd.read_csv(CDIR + 'Int.J.Mol.Sci.2022/r56_metabolomics_yali_v2.csv', index_col=0)
ph = pd.read_table(CDIR + 'Int.J.Mol.Sci.2022/r56_pheno_text_binned_v2.txt', index_col=0)
Vc = d.iloc[:, 3:].apply(pd.to_numeric, errors='coerce').T
Vc.columns = d.Molecule.astype(str).values
ph = ph.loc[ph.index.intersection(Vc.index)].copy()
print('Che phenotypes', ph.Phenotype.value_counts().to_dict())
ph['y'] = (ph.Phenotype == 'case').astype(int)
che_map = {'indolepropionate': ['Indole-3-propionic acid', 'indole_3_propionic_acid_PM_a1'],
           'indolelactate': ['DL-Indole-3-lactic acid', 'indole_3_lactate_PM_a1'],
           'indoleacetate': ['indole_3_acetate_PM_a1'], 'glycodeoxycholate (unsulfated)': ['Glycodeoxycholic acid'],
           'kynurenine': ['Kynurenine']}
for k, cols in che_map.items():
    for c in cols:
        x = Vc.loc[ph.index, c]
        if isinstance(x, pd.DataFrame): x = x.iloc[:, 0]
        x = x.astype(float); mn = x[x > 0].min(); x = np.log2(x.where(x > 0, mn / 2))
        ph['v'] = (x - x.mean()) / x.std()
        m = smf.ols('v ~ y + C(sex) + C(age_binned) + C(bmi_binned) + C(site) + C(season_a1)', ph).fit()
        ci = m.conf_int().loc['y']
        rows.append(dict(cohort='Che2022', metabolite=f'{k} [{c}]', n_case=int(ph.y.sum()), n_ctrl=int((1 - ph.y).sum()),
                         beta_SD=m.params['y'], lo=ci[0], hi=ci[1], p=m.pvalues['y']))
O = pd.DataFrame(rows); O['q_within_cohort'] = O.groupby('cohort').p.transform(lambda p: bh(p.values))
O.to_csv(OUT + 'external_replication.csv', index=False)
pd.set_option('display.width', 220); print(O.round(3).to_string())
