"""Build FULL (all-metabolite) per-study matrices for ME/CFS blood metabolomics cohorts.
Output /tmp/work_mw/studies_full.pkl : {study: dict(L, C, ann, D, family, person_level)}
  L   : samples(persons) x features, log2, half-min imputed, features detected in >=50% of samples
  D   : samples x ALL features, detected (bool) -- used for drug/xenobiotic exposure calls
  C   : covariates (case, age, male, bmi, ...)
  ann : per-feature annotation (name, comp_id, inchikey, hmdb, kegg, pubchem, sub, super, detect)
Loader logic follows downloads/mwb/01_build.py (previous round), extended to keep every metabolite.
"""
import re, sys, pickle, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, '/mnt/project-files/mecfs-research/downloads/mwb/')
from mwparse import ssf, data, RAW
sys.path.insert(0, '/mnt/project-files/biomap/deep/stats_replication/')
from data import load_meta

BM = '/mnt/project-files/biomap/BioMapAI/'
IC = BM + 'data/indenpendent_cohorts/'
S = {}


def clean_id(x):
    x = str(x).strip()
    return '' if x.lower() in ('nan', '', '0', 'none', '-', 'n/a') else x


def mk_ann(names, **cols):
    a = pd.DataFrame({'name': list(names)})
    for k in ['comp_id', 'inchikey', 'hmdb', 'kegg', 'pubchem', 'sub', 'super', 'chebi']:
        v = cols.get(k)
        a[k] = [clean_id(x) for x in v] if v is not None else ''
    return a


def finish(study, V, ann, C, family, platform, cov, note='', imputed_min=False, det_thr=0.5):
    """V raw-scale samples x features (column names = ann.name, unique)."""
    ann = ann.drop_duplicates('name').set_index('name')
    V = V.loc[:, ~V.columns.duplicated()]
    V = V.loc[:, [c for c in V.columns if 'istd' not in c.lower() and c in ann.index]]
    V = V.where(V > 0)
    if imputed_min:   # data already imputed with the minimum: detected = above column minimum
        D = V.gt(V.min() * (1 + 1e-6))
    else:
        D = V.notna()
    C = C.loc[V.index]
    det = D.mean()
    keep = (det >= det_thr) & (V.notna().sum() > 5)
    Vk = V.loc[:, keep]
    L = np.log2(Vk.fillna(Vk.min() / 2))
    L = L.loc[:, L.std() > 1e-9]
    ann = ann.copy(); ann['detect'] = det.reindex(ann.index)
    S[study] = dict(L=L, D=D, C=C, ann=ann, family=family, platform=platform, cov=cov, note=note)
    print(f'{study:28s} n={len(L)} (case {int(C.case.sum())}/ctrl {int((1-C.case).sum())}) features all={V.shape[1]} kept={L.shape[1]}', flush=True)


# ============ 1. BioMapAI (Xiong 2025) QC-normalised Metabolon; averaged per PERSON over visits (geometric mean)
q = pd.read_csv(BM + 'data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv')
q = q[q.TYPE == 'NAMED']
Sm = [c for c in q.columns if c.startswith('SAM')]
V = q[Sm].astype(float).T; V.columns = q.CHEMICAL_NAME.values
ann = mk_ann(q.CHEMICAL_NAME, comp_id=q.COMP_ID, inchikey=q.INCHIKEY, hmdb=q.HMDB, kegg=q.KEGG, pubchem=q.PUBCHEM, sub=q.SUB_PATHWAY, super=q.SUPER_PATHWAY)
m = load_meta(); idx = V.index.intersection(m.index); V = V.loc[idx]; M = m.loc[idx]
# detection per visit (imputed minimum) then per-person: detected in any visit
Dv = V.gt(V.min() * (1 + 1e-6))
Lv = np.log2(V)
Lp = Lv.groupby(M.person.values).mean()
Dp = Dv.groupby(M.person.values).mean()   # fraction of visits detected
Mp = M.drop_duplicates('person').set_index('person')
num = M.groupby('person')[['severity', 'NLT.Orthostatic Intolerance', 'RAND36.Physical Functioning', 'RAND36.Energy/Fatigue',
                           'DANA.Cognitive Efficiency', 'PSQI.sleeping_problem', 'GSRS.gastrointestinal', 'diet_meat', 'diet_sugar', 'diet_veg', 'diet_grains', 'diet_fruit']].mean()
yes = lambda c: M.assign(v=(M[c] == 'Yes').astype(float)).groupby('person').v.max()
C = pd.DataFrame({'case': Mp.case, 'age': Mp.age, 'male': Mp.male, 'bmi': Mp.bmi_imp, 'bmi_missing': Mp.bmi_missing,
                  'short': Mp.short, 'nvis': M.groupby('person').size()}).join(num)
for c in ['antibiotics', 'probiotics', 'antivirals', 'antifungals']: C[c] = yes(c)
C = C.loc[Lp.index]
Vp = np.exp2(Lp)
Dbool = Dp > 0
# keep: detected in >=50% of visits overall (same rule)
det_thr_cols = Dv.mean() >= 0.5
L = Lp.loc[:, det_thr_cols.values]; L = L.loc[:, L.std() > 1e-9]
annb = ann.drop_duplicates('name').set_index('name'); annb['detect'] = Dv.mean().reindex(annb.index)
S['BioMapAI'] = dict(L=L, D=Dbool, C=C, ann=annb, family='Metabolon', platform='Metabolon HD4 (QC-normalised)', cov='age, sex, BMI (person mean of visits)', note='person-level geometric mean over 1-4 visits',
                     Lvisit=Lv, Mvisit=M)
print(f"BioMapAI n={len(L)} case {int(C.case.sum())} feats kept {L.shape[1]}")

# ============ 2. Germain 2022 JCI Insight: pre-exercise D1/D2 averaged per person (Metabolon)
d = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='Original Dataset', header=None)
fm = d.iloc[8:, :15].copy(); fm.columns = d.iloc[7, :15].values
gm = d.iloc[:7, 16:].T.copy(); gm.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'AgeBin', 'BMIBin', 'Bell']; gm = gm.reset_index(drop=True)
V = d.iloc[8:, 16:].apply(pd.to_numeric, errors='coerce').T.reset_index(drop=True); V.columns = fm.BIOCHEMICAL.astype(str).values
pre = gm.Timepoint.astype(str).str.endswith('PRE').values
Vp_ = V[pre].copy(); gp = gm[pre]; Vp_.index = gp.ID.astype(str).values
Dg = Vp_.notna().groupby(level=0).max()
Vg = np.exp2(np.log2(Vp_.where(Vp_ > 0)).groupby(level=0).mean())
g1 = gp.drop_duplicates('ID'); g1 = g1.set_index(g1.ID.astype(str)).loc[Vg.index]
C = pd.DataFrame({'case': (g1.Phenotype == 'CFS').astype(int), 'age': g1.AgeBin.astype(float), 'male': (g1.Sex == 'M').astype(int),
                  'bmi': g1.BMIBin.astype(float), 'bell': pd.to_numeric(g1.Bell, errors='coerce')}, index=Vg.index)
ann = mk_ann(fm.BIOCHEMICAL.astype(str), comp_id=fm['COMP ID'], hmdb=fm.HMDB, kegg=fm.KEGG, pubchem=fm.PUBCHEM, sub=fm['SUB PATHWAY'], super=fm['SUPER PATHWAY'])
# also keep post-exercise matrices for possible use
finish('Germain2022', Vg, ann, C, 'Metabolon', 'Metabolon HD4', 'age bin, sex, BMI bin', 'mean of D1-PRE and D2-PRE')
S['Germain2022']['D'] = Dg.loc[S['Germain2022']['L'].index]

# ============ 3/4. Germain 2018, 2020 Metabolites (Metabolon; Europe PMC supplements; no covariates)
for study, f in [('Germain2018', RAW + 'Germain2018_Metabolites/Supplementary_File_1.xlsx'), ('Germain2020', RAW + 'Germain2020_Metabolites/Supplementary File 1.xlsx')]:
    g = pd.read_excel(f, header=0); g.columns = [str(c).strip() for c in g.columns]
    samp = [c for c in g.columns if re.fullmatch(r'[CP]\d+', c)]
    V = g[samp].apply(pd.to_numeric, errors='coerce').T; V.columns = g.BIOCHEMICAL.astype(str).values
    hm = [c for c in g.columns if 'HMDB' in c][0]
    gc = lambda k: g[[c for c in g.columns if c.replace('_', ' ') == k][0]]
    ann = mk_ann(g.BIOCHEMICAL.astype(str), comp_id=gc('COMP ID'), hmdb=g[hm], kegg=g.KEGG, pubchem=g.PUBCHEM, sub=gc('SUB PATHWAY'), super=gc('SUPER PATHWAY'))
    C = pd.DataFrame({'case': [int(s.startswith('P')) for s in V.index]}, index=V.index)
    finish(study, V, ann, C, 'Metabolon', 'Metabolon', 'none', 'Hanson lab; may overlap other Hanson cohorts')

# ============ 5. Che 2022 (MCAM, Fiehn WCMC)
che = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_metabolomics_yali_v2.csv', index_col=0)
ph = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_pheno_text_binned_v2.txt', sep='\t', index_col=0)
che = che[che.Molecule.notna()]
che['Molecule'] = che.Molecule.astype(str).str.strip()
che = che.drop_duplicates('Molecule')
V = che.iloc[:, 3:].apply(pd.to_numeric, errors='coerce').T; V.columns = che.Molecule.values
dbi = che.database_identifier.astype(str)
ann = mk_ann(che.Molecule, inchikey=che.InChiKey, hmdb=[x if x.startswith('HMDB') else '' for x in dbi],
             kegg=[x if re.fullmatch(r'C\d{5}', x) else '' for x in dbi], pubchem=[x.replace('CID', '') if 'CID' in x else '' for x in dbi], chebi=[x if 'CHEBI' in x else '' for x in dbi])
ph = ph.loc[V.index]
agemap = {'18 - 24': 21, '25 - 34': 30, '35 - 44': 40, '45 - 54': 50, '55 - 64': 60, '65 and over': 70}
bmimap = {'underweight': 18, 'healthy weight': 22, 'overweight': 27, 'obese': 33}
C = pd.DataFrame({'case': (ph.Phenotype == 'case').astype(int), 'age': ph.age_binned.map(agemap), 'male': (ph.sex == 'male').astype(int), 'bmi': ph.bmi_binned.map(bmimap),
                  'ibs': (ph.ibs_dx.astype(str).str.lower().str.startswith('y')).astype(int)}, index=ph.index)
C['age'] = C.age.fillna(C.age.median()); C['bmi'] = C.bmi.fillna(C.bmi.median())
for s in sorted(ph.site.unique())[1:]: C['site_' + str(s)] = (ph.site == s).astype(int)
finish('Che2022', V, ann, C, 'Fiehn', 'Fiehn GC-TOF/HILIC/CSH', 'age bin, sex, BMI bin, site')

# ============ 6. Nagy-Szakal 2018 (Fiehn; MW ST000800/801/813)
s = ssf(RAW + 'AN001271.mwtab.txt'); s = s[s.CaseStatus.isin(['Case', 'Control'])]
num2esp = {int(re.search(r'ESP-\d+_(\d+)', r.subject).group(1)): r['sample'] for _, r in s.iterrows()}
blocks, anns = [], []
for f in ['ST000800_AN001271_Results.txt', 'ST000800_AN001272_Results.txt']:
    t = pd.read_csv(RAW + f, sep='\t', low_memory=False); t = t[t.Annotation.notna()]
    t['nm'] = t.Annotation.astype(str).str.strip().str.replace(r'^1_\s*', '', regex=True)
    t = t[~t.nm.str.lower().str.contains('istd')].drop_duplicates('nm')
    cols = [c for c in t.columns if re.fullmatch(r'\d+', str(c))]
    V = t[cols].apply(pd.to_numeric, errors='coerce').T; V.columns = t.nm.values
    V.index = [num2esp.get(int(c)) for c in cols]; V = V[V.index.notna()]; V = V.groupby(level=0).mean()
    blocks.append(V); anns.append(mk_ann(t.nm, inchikey=t['InChI Key']))
t = pd.read_csv(RAW + 'ST000801_AN001273_Results.txt', sep='\t', low_memory=False); t = t[t.Annotation.notna()]
t['nm'] = t.Annotation.astype(str).str.strip(); t = t.drop_duplicates('nm')
cols = [c for c in t.columns if 'ESP-' in c]
V = t[cols].apply(pd.to_numeric, errors='coerce').T; V.columns = t.nm.values
V.index = [re.search(r'(ESP-\d+)', c).group(1) for c in cols]; V = V.groupby(level=0).mean(); blocks.append(V)
anns.append(mk_ann(t.nm, pubchem=[str(int(x)) if pd.notna(x) and x > 1 else '' for x in t['PubChem CID']]))
t = pd.read_csv(RAW + 'ST000813_AN001286_Results.txt', sep='\t', low_memory=False)
t = t[t[' Metabolite Name'].astype(str).str.contains('[a-z]', regex=True) & (t[' Metabolite Name'] != 'Factors')]
t['nm'] = t[' Metabolite Name'].astype(str).str.strip(); t = t.drop_duplicates('nm')
cols = [c for c in t.columns if c.startswith('ESP-')]
V = t[cols].apply(pd.to_numeric, errors='coerce').T; V.columns = t.nm.values
V.index = [re.search(r'(ESP-\d+)', c).group(1) for c in cols]; V = V.groupby(level=0).mean(); blocks.append(V)
anns.append(mk_ann(t.nm, inchikey=t['InChI Key'], kegg=t.KEGG, pubchem=t.PubChem))
order = [3, 2, 1, 0]   # GC-TOF first, then HILIC, then CSH
V = pd.concat([blocks[i] for i in order], axis=1, join='inner'); ann = pd.concat([anns[i] for i in order])
V = V.loc[:, ~V.columns.duplicated()]
s2 = s.set_index('sample').loc[V.index]
C = pd.DataFrame({'case': (s2.CaseStatus == 'Case').astype(int), 'age_over50': (s2.AGE_DICH == 'over50').astype(int), 'male': (s2.SEX == 'Male').astype(int),
                  'bmi': pd.to_numeric(s2.x_bmi, errors='coerce'), 'ibs': (s2.x_IBS.astype(str).str.lower().str.contains('yes|ibs\\+|^1')).astype(int)}, index=V.index)
C['bmi'] = C.bmi.fillna(C.bmi.median())
finish('NagySzakal2018', V, ann, C, 'Fiehn', 'Fiehn GC-TOF + HILIC + CSH', 'age >50, sex, BMI')

# ============ 7. Naviaux 2016 (targeted LC-MS/MS)
Vs, ans = [], []
for a in ['AN000705', 'AN000706']:
    V, an = data(RAW + a + '.mwtab.txt'); Vs.append(V); ans.append(an)
V = pd.concat(Vs, axis=1); an = pd.concat(ans).drop_duplicates('metabolite_name')
ann = mk_ann(an.metabolite_name.str.strip(), kegg=an['KEGG ID'], pubchem=an['PubChem ID'])
V.columns = [c.strip() for c in V.columns]
s = ssf(RAW + 'AN000705.mwtab.txt').set_index('sample').loc[V.index]
C = pd.DataFrame({'case': (s.Disease == 'CFS').astype(int), 'male': (s.gender == 'Male').astype(int)}, index=V.index)
finish('Naviaux2016', V, ann, C, 'Other', 'targeted LC-MS/MS (Naviaux)', 'sex')

# ============ 8. NCNP Japan serum tryptophan panel
V, an = data(RAW + 'AN008373.mwtab.txt')
ann = mk_ann(an.metabolite_name.str.strip(), kegg=an['KEGG ID'], pubchem=an['PubChem ID']); V.columns = [c.strip() for c in V.columns]
s = ssf(RAW + 'AN008373.mwtab.txt').set_index('sample').loc[V.index]
C = pd.DataFrame({'case': (s.Condition == 'ME/CFS').astype(int)}, index=V.index)
finish('NCNP2026', V, ann, C, 'Other', 'targeted LC-MS/MS tryptophan panel', 'none')

# ============ 9. Armstrong 2015 NMR (MTBLS161)
mf = pd.read_csv(RAW + 'MTBLS161/m_MTBLS161_mecfs_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv', sep='\t')
sc = [c for c in mf.columns if c.endswith('serum')]
V = mf[sc].apply(lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')).T; V.columns = mf.metabolite_identification.astype(str).str.strip().values
ann = mk_ann(V.columns, chebi=mf.database_identifier.astype(str))
C = pd.DataFrame({'case': [int(i.startswith('CFS')) for i in V.index]}, index=V.index)
finish('Armstrong2015', V, ann, C, 'Other', '1H-NMR serum', 'none')

pickle.dump(S, open('/tmp/work_mw/studies_full.pkl', 'wb'))
print('saved')
