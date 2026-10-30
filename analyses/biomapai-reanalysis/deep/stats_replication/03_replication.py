"""Replication of BioMapAI case-control signals in bundled independent cohorts + simple meta-analysis,
and severity (Bell score) tests in Germain 2022.
Effect sizes: standardised mean difference (feature z-scored within cohort), covariate-adjusted
(OLS, or random-intercept LMM when a cohort has repeated samples)."""
import numpy as np, pandas as pd
from scipy import stats
from data import *
from lmmfast import fit_many, bh

IC = B + 'data/indenpendent_cohorts/'
BP_GEN = ['Faecalibacterium', 'Roseburia', 'Coprococcus', 'Anaerostipes', 'Agathobacter']
SPECIES = ['Faecalibacterium_prausnitzii', 'Roseburia_faecis', 'Roseburia_intestinalis', 'Roseburia_hominis', 'Roseburia_inulinivorans',
           'Coprococcus_eutactus', 'Coprococcus_comes', 'Coprococcus_catus', 'Anaerostipes_hadrus',
           'Ruminococcus_gnavus', 'Flavonifractor_plautii', 'Dysosmobacter_welbionis', 'Enterocloster_bolteae', 'Clostridium_symbiosum']
METABS = ['kynurenine', 'tryptophan', 'vanillylmandelate (VMA)', 'choline', '3-hydroxybutyrate (BHBA)', 'kynurenate',
          'glycocholate', 'glycochenodeoxycholate', 'glycodeoxycholate', 'deoxycholate', 'cholate', 'chenodeoxycholate', 'lithocholate',
          'caffeine', 'paraxanthine', 'quinate', 'hippurate', 'catechol sulfate', 'palmitate (16:0)', 'oleate/vaccenate (18:1)', 'linoleate (18:2n6)']


# ---------------------------------------------------------------- helpers
def zs(X): return (X - X.mean()) / X.std()


def ols_many(Y, C, col=1):
    """Y (n,K) DataFrame, C (n,p) with target at column col. Returns est, se, p arrays."""
    r = fit_many(Y.values, C, np.arange(len(C)), L=np.eye(C.shape[1])[[col]], grid=[0.0])
    return r['est'][:, 0], r['se'][:, 0], r['p'][:, 0]


def species_composites(rel):
    rel = rel.div(rel.sum(axis=1), axis=0)
    bp = [c for c in rel.columns if any(c.startswith(g) for g in BP_GEN)]
    return pd.DataFrame({'COMPOSITE butyrate producers (log10 rel. abund.)': np.log10(rel[bp].sum(axis=1) + 1e-4)})


def clr_keep(rel, keep_extra):
    rel = rel.div(rel.sum(axis=1), axis=0)
    keep = ((rel > 1e-6).mean() >= 0.1) | rel.columns.isin(keep_extra)
    r = rel.loc[:, keep]; pc = r[r > 0].min().min() / 2
    lg = np.log(r + pc); return lg.sub(lg.mean(axis=1), axis=0)


def metabolon_composites(X, sub):
    """X: log2 samples x compounds, sub: Series compound->SUB_PATHWAY (Metabolon nomenclature)."""
    sub = sub.reindex(X.columns).fillna('')
    z = zs(X)
    g = {'COMPOSITE long-chain free fatty acids': sub.isin(['Long Chain Polyunsaturated Fatty Acid (n3 and n6)', 'Long Chain Saturated Fatty Acid', 'Long Chain Monounsaturated Fatty Acid']),
         'COMPOSITE long-chain acylcarnitines': sub.isin(['Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)']),
         'COMPOSITE ceramides': sub.isin(['Ceramides']),
         'COMPOSITE sphingomyelins': sub.isin(['Sphingomyelins']),
         'COMPOSITE primary bile acids': sub.isin(['Primary Bile Acid Metabolism']),
         'COMPOSITE secondary bile acids': sub.str.startswith('Secondary Bile Acid')}
    out = {k: z.loc[:, v.values].mean(axis=1) for k, v in g.items() if v.sum() >= 2}
    if 'kynurenine' in X and 'tryptophan' in X: out['kynurenine/tryptophan (log2 ratio)'] = X['kynurenine'] - X['tryptophan']
    return pd.DataFrame(out), {k: int(v.sum()) for k, v in g.items()}


# ---------------------------------------------------------------- BioMapAI reference
m = load_meta()
ref = pd.read_csv(OUT + 'casecontrol_all_models.csv.gz')
ref = ref[ref.model == 'LMM_all'].pivot_table(index='feature', columns='contrast', values=['est', 'se', 'p'])
biomap = {}


def biomap_fit(X):
    idx = X.index.intersection(m.index); M = m.loc[idx]; X = zs(X.loc[idx])
    cov = [M.age_c, M.male, M.bmi_c, M.bmi_missing, M.tp2, M.tp3]
    out = {}
    C = np.column_stack([np.ones(len(M)), M.case] + cov)
    r = fit_many(X.values, C, M.person.values, L=np.eye(C.shape[1])[[1]])
    C2 = np.column_stack([np.ones(len(M)), M.short, M.long] + cov)
    r2 = fit_many(X.values, C2, M.person.values, L=np.eye(C2.shape[1])[[1, 2]])
    for k, f in enumerate(X.columns):
        out[f] = dict(est=r['est'][k, 0], se=r['se'][k, 0], p=r['p'][k, 0], short_est=r2['est'][k, 0], short_p=r2['p'][k, 0], long_est=r2['est'][k, 1], long_p=r2['p'][k, 1],
                      n_case=int(M.drop_duplicates('person').case.sum()), n_ctrl=int((1 - M.drop_duplicates('person').case).sum()))
    return out


rel_bm = pd.read_csv(B + 'data/metagenomics/metaphlan4_species_full.csv', index_col=0).select_dtypes('number').T
biomap.update(biomap_fit(species_composites(rel_bm)))
Xm = load_layer('metabolome'); mm = metab_meta()
comp_bm, nbm = metabolon_composites(Xm, mm.SUB_PATHWAY); print('BioMapAI composite sizes', nbm)
biomap.update(biomap_fit(comp_bm))
top_sp = ref['p']['patient_vs_control'][ref.index.isin(clr(rel_bm).columns)].nsmallest(10).index.tolist()
top_mb = ref['p']['patient_vs_control'][ref.index.isin(Xm.columns)].nsmallest(15).index.tolist()
SPECIES_ALL = list(dict.fromkeys(SPECIES + top_sp)); METABS_ALL = list(dict.fromkeys(METABS + top_mb))
for f in SPECIES_ALL + METABS_ALL:
    if f in ref.index:
        n = biomap.get('COMPOSITE long-chain free fatty acids', {})
        biomap[f] = dict(est=ref.loc[f, ('est', 'patient_vs_control')], se=ref.loc[f, ('se', 'patient_vs_control')], p=ref.loc[f, ('p', 'patient_vs_control')],
                         short_est=ref.loc[f, ('est', 'short_vs_control')], short_p=ref.loc[f, ('p', 'short_vs_control')],
                         long_est=ref.loc[f, ('est', 'long_vs_control')], long_p=ref.loc[f, ('p', 'long_vs_control')])

ext = []  # rows: feature, cohort, est, se, p, n_case, n_ctrl


def record(cohort, X, C, case, groups=None):
    X = zs(X.loc[:, X.std() > 0])
    if groups is None:
        e, s, p = ols_many(X, C)
    else:
        r = fit_many(X.values, C, groups, L=np.eye(C.shape[1])[[1]]); e, s, p = r['est'][:, 0], r['se'][:, 0], r['p'][:, 0]
    nc = len(set(np.asarray(groups)[case == 1])) if groups is not None else int(case.sum())
    nh = len(set(np.asarray(groups)[case == 0])) if groups is not None else int((1 - case).sum())
    for k, f in enumerate(X.columns): ext.append((f, cohort, e[k], s[k], p[k], nc, nh))


# ---------------------------------------------------------------- Raijmakers 2020 (NL stool, no covariates bundled)
r = pd.read_csv(IC + 'Raijmakers/metaphlan4_species_full.csv', index_col=0).select_dtypes('number').T
r = r[r.index.str.startswith('CFS') | r.index.str.startswith('HV')]
case = r.index.str.startswith('CFS').astype(int)
X = clr_keep(r, SPECIES_ALL)[[s for s in SPECIES_ALL if s in r.columns]].join(species_composites(r))
record('Raijmakers2020_stool', X, np.column_stack([np.ones(len(X)), case]), case)

# ---------------------------------------------------------------- Columbia / Guo 2023 (US stool): age, sex, BMI, site
c = pd.read_csv(IC + 'columbia/metaphlan4_species_full.csv', index_col=0).select_dtypes('number').T
sra = pd.read_csv(IC + 'columbia/SraRunTable.txt').set_index('Run').loc[c.index]
case = (sra.sampletype == 'Case').astype(int).values
site = pd.get_dummies(sra.Site, drop_first=True).astype(float).values
C = np.column_stack([np.ones(len(c)), case, (sra.age_at_a1 - 45) / 10, (sra.sex == 'male').astype(float), (sra.BMI - 26) / 5, site])
X = clr_keep(c, SPECIES_ALL)[[s for s in SPECIES_ALL if s in c.columns]].join(species_composites(c))
record('Guo2023_Columbia_stool', X, C, case)

# ---------------------------------------------------------------- Che 2022 (US plasma, same participants as Columbia stool)
che = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_metabolomics_yali_v2.csv', index_col=0)
ph = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_pheno_text_binned_v2.txt', sep='\t', index_col=0)
V = che.iloc[:, 3:].apply(pd.to_numeric, errors='coerce').T
V.columns = che.Molecule.astype(str).values
V = V.loc[:, ~V.columns.duplicated()]
V = V.loc[:, (V > 0).mean() >= 0.5]
V = np.log2(V + V[V > 0].min() / 2)
ph = ph.loc[V.index]
# crosswalk BioMapAI name -> Che name: manual for key compounds, else InChIKey first block
manual = {'kynurenine': 'Kynurenine', 'tryptophan': 'Tryptophan', 'choline': 'Choline cation', '3-hydroxybutyrate (BHBA)': '3_hydroxybutyric_acid_PM_a1',
          'caffeine': 'Caffeine', 'glycocholate': 'Glycocholic acid', 'glycodeoxycholate': 'Glycodeoxycholic acid',
          'palmitate (16:0)': 'palmitic_acid_PM_a1', 'oleate/vaccenate (18:1)': 'FA (18:1) (oleic acid)', 'linoleate (18:2n6)': 'FA (18:2) (linoleic acid)'}
ik = che.drop_duplicates('Molecule').set_index(che.drop_duplicates('Molecule').InChiKey.astype(str).str[:14]).Molecule
xw = {}
for f in METABS_ALL:
    if f in manual and manual[f] in V: xw[f] = manual[f]; continue
    k = str(mm.INCHIKEY.get(f, ''))[:14]
    if k and k in ik.index:
        cand = ik.loc[[k]].tolist(); cand = [x for x in cand if x in V]
        if cand: xw[f] = cand[0]
print('Che crosswalk', xw)
nmc = pd.Series(V.columns, index=V.columns)
cheg = {'COMPOSITE long-chain free fatty acids': [x for x in V.columns if (x.startswith('FA (') and int(x[4:6]) >= 14) or x in ['palmitic_acid_PM_a1', 'stearic_acid_PM_a1', 'arachidic_acid_PM_a1', 'behenic_acid_PM_a1', 'myristic_acid_PM_a1', 'oleic_acid_PM_a1', 'linoleic_acid_PM_a1', 'palmitoleic_acid_PM_a1']],
        'COMPOSITE long-chain acylcarnitines': [x for x in V.columns if x in ['Linoleoylcarnitine', 'Oleoyl-L-carnitine', '3-Hydroxyoleylcarnitine', 'Palmitoylcarnitine', 'Stearoylcarnitine', 'Myristoylcarnitine']],
        'COMPOSITE ceramides': [x for x in V.columns if x.startswith('Cer(')],
        'COMPOSITE sphingomyelins': [x for x in V.columns if x.startswith('SM (')],
        'COMPOSITE primary bile acids': [x for x in V.columns if x in ['Glycocholic acid', 'Cholic acid', 'Glycochenodeoxycholic acid', 'Taurocholic acid']],
        'COMPOSITE secondary bile acids': [x for x in V.columns if x in ['Glycodeoxycholic acid', 'Tauroursodeoxycholic acid', 'Deoxycholic acid', 'Lithocholic acid']]}
print('Che composite sizes', {k: len(v) for k, v in cheg.items()})
Z = zs(V)
Xc = pd.DataFrame({f: V[x] for f, x in xw.items()})
for k, v in cheg.items():
    if len(v) >= 1: Xc[k] = Z[v].mean(axis=1)
if 'Kynurenine' in V and 'Tryptophan' in V: Xc['kynurenine/tryptophan (log2 ratio)'] = V['Kynurenine'] - V['Tryptophan']
agemap = {'18 - 24': 21, '25 - 34': 30, '35 - 44': 40, '45 - 54': 50, '55 - 64': 60, '65 and over': 70}
bmimap = {'underweight': 0, 'healthy weight': 1, 'overweight': 2, 'obese': 3}
case = (ph.Phenotype == 'case').astype(int).values
C = np.column_stack([np.ones(len(ph)), case, (ph.age_binned.map(agemap) - 45) / 10, (ph.sex == 'male').astype(float), ph.bmi_binned.map(bmimap),
                     pd.get_dummies(ph.site, drop_first=True).astype(float).values])
record('Che2022_Columbia_plasma', Xc, C, case)

# ---------------------------------------------------------------- Germain 2022 (US plasma, Metabolon; pre-exercise samples D1-PRE and D2-PRE)
d = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='ScaledImpDataZeroDrug&Tobacco', header=None)
fm = d.iloc[8:, :15].copy(); fm.columns = d.iloc[7, :15].values
gm = d.iloc[:7, 15:].T; gm.columns = gm.iloc[0]; gm = gm.iloc[1:]
gm.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'AgeBin', 'BMIBin', 'Bell']
G = d.iloc[8:, 16:].astype(float).T; G.columns = fm.BIOCHEMICAL.astype(str).values; G.index = gm.index
G = G.loc[:, ~G.columns.duplicated()]
G = np.log2(G)
G = G.loc[:, (G > G.min() + 1e-9).mean() >= 0.2]
gsub = pd.Series(fm['SUB PATHWAY'].astype(str).values, index=fm.BIOCHEMICAL.astype(str).values)
gsub = gsub[~gsub.index.duplicated()]
pre = gm.Timepoint.str.endswith('PRE').values
Gp = G[pre]; gp = gm[pre]
comp_g, ng = metabolon_composites(Gp, gsub); print('Germain composite sizes', ng)
Xg = Gp[[f for f in METABS_ALL if f in Gp.columns]].join(comp_g)
case = (gp.Phenotype == 'CFS').astype(int).values
C = np.column_stack([np.ones(len(gp)), case, gp.AgeBin.astype(float), (gp.Sex == 'M').astype(float), gp.BMIBin.astype(float), (gp.Timepoint == 'D2-PRE').astype(float)])
record('Germain2022_plasma', Xg, C, case, groups=gp.ID.values)

# ---------------------------------------------------------------- assemble + meta-analysis
E = pd.DataFrame(ext, columns=['feature', 'cohort', 'est', 'se', 'p', 'n_case', 'n_ctrl'])
E.to_csv(OUT + 'replication_per_cohort.csv', index=False)


def meta(e, s):
    w = 1 / s ** 2; fe = np.sum(w * e) / w.sum(); fse = np.sqrt(1 / w.sum())
    k = len(e); Q = np.sum(w * (e - fe) ** 2)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (s ** 2 + tau2); re = np.sum(wr * e) / wr.sum(); rse = np.sqrt(1 / wr.sum())
    I2 = max(0, (Q - (k - 1)) / Q) if (k > 1 and Q > 0) else 0
    return fe, fse, 2 * stats.norm.sf(abs(fe / fse)), re, rse, 2 * stats.norm.sf(abs(re / rse)), I2


rows = []
for f in list(dict.fromkeys([*SPECIES_ALL, 'COMPOSITE butyrate producers (log10 rel. abund.)', *METABS_ALL, *comp_bm.columns])):
    b = biomap.get(f); e = E[E.feature == f]
    if b is None and len(e) == 0: continue
    row = dict(feature=f, layer='stool' if (f in SPECIES_ALL or 'butyrate' in f) else 'plasma',
               BioMapAI_est=b['est'] if b else np.nan, BioMapAI_p=b['p'] if b else np.nan,
               BioMapAI_short_est=b['short_est'] if b else np.nan, BioMapAI_short_p=b['short_p'] if b else np.nan,
               BioMapAI_long_est=b['long_est'] if b else np.nan, BioMapAI_long_p=b['long_p'] if b else np.nan)
    for _, x in e.iterrows():
        row[x.cohort + '_est'] = x.est; row[x.cohort + '_p'] = x.p
    if len(e):
        fe, fse, fp, re, rse, rp, I2 = meta(e.est.values, e.se.values)
        row.update(ext_k=len(e), ext_meta_FE=fe, ext_meta_FE_lo=fe - 1.96 * fse, ext_meta_FE_hi=fe + 1.96 * fse, ext_meta_FE_p=fp,
                   ext_meta_RE=re, ext_meta_RE_p=rp, ext_I2=I2)
        if b:
            ee = np.r_[e.est.values, b['est']]; ss = np.r_[e.se.values, b['se']]
            fe2, _, fp2, re2, _, rp2, I22 = meta(ee, ss)
            row.update(all_meta_RE=re2, all_meta_RE_p=rp2, all_I2=I22)
            # discovery signal = BioMapAI patient-vs-control, or the short/long contrast if only that one was significant
            if b['p'] < 0.05: ref_sign, disc = np.sign(b['est']), 'patient'
            elif b['short_p'] < 0.05: ref_sign, disc = np.sign(b['short_est']), 'short'
            elif b['long_p'] < 0.05: ref_sign, disc = np.sign(b['long_est']), 'long'
            else: ref_sign, disc = np.sign(fe), 'none'
            row['BioMapAI_discovery_contrast'] = disc
            row['n_ext_same_dir'] = int((np.sign(e.est) == ref_sign).sum())
            row['n_ext_same_dir_p05'] = int(((np.sign(e.est) == ref_sign) & (e.p < 0.05)).sum())
            if disc == 'none':
                row['replicates'] = 'not a BioMapAI signal' + (' (external meta p<0.05)' if fp < 0.05 else '')
            else:
                row['replicates'] = ('yes' if (np.sign(fe) == ref_sign and fp < 0.05) else ('direction only' if np.sign(fe) == ref_sign else 'no (opposite)'))
    rows.append(row)
T = pd.DataFrame(rows)
T['in_top_list'] = T.feature.isin(top_sp + top_mb)
T.to_csv(OUT + 'replication_table.csv', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 30)
show = ['feature', 'BioMapAI_est', 'BioMapAI_p', 'BioMapAI_short_est', 'BioMapAI_long_est', 'BioMapAI_long_p', 'BioMapAI_discovery_contrast', 'Raijmakers2020_stool_est', 'Raijmakers2020_stool_p', 'Guo2023_Columbia_stool_est', 'Guo2023_Columbia_stool_p',
        'Che2022_Columbia_plasma_est', 'Che2022_Columbia_plasma_p', 'Germain2022_plasma_est', 'Germain2022_plasma_p', 'ext_meta_FE', 'ext_meta_FE_p', 'ext_I2', 'replicates']
print(T[[c for c in show if c in T]].round(3).to_string(index=False))
print('N per cohort', E.groupby('cohort')[['n_case', 'n_ctrl']].first().to_string())

# ---------------------------------------------------------------- Germain: severity (Bell score, higher = better function) in CFS cases
cf = gp[gp.Phenotype == 'CFS']; Xs = zs(Xg.loc[cf.index])
bell = cf.Bell.astype(float)
res = []
C = np.column_stack([np.ones(len(cf)), (bell - bell.mean()) / bell.std(), cf.AgeBin.astype(float), (cf.Sex == 'M').astype(float), cf.BMIBin.astype(float), (cf.Timepoint == 'D2-PRE').astype(float)])
r = fit_many(Xs.values, C, cf.ID.values, L=np.eye(C.shape[1])[[1]])
pm = Xs.groupby(cf.ID.values).mean(); bm = bell.groupby(cf.ID.values).mean()
for k, f in enumerate(Xs.columns):
    rho, pr = stats.spearmanr(pm[f], bm.loc[pm.index])
    res.append((f, len(bm), bell.min(), bell.max(), r['est'][k, 0], r['se'][k, 0], r['p'][k, 0], rho, pr))
S = pd.DataFrame(res, columns=['feature', 'n_cases', 'bell_min', 'bell_max', 'est_SD_per_SD_Bell', 'se', 'p', 'spearman_rho_personmean', 'spearman_p'])
S.to_csv(OUT + 'germain_bell_severity.csv', index=False)
print(S[S.feature.str.contains('kynuren|COMPOSITE|BHB|tryptophan')].round(3).to_string(index=False))
