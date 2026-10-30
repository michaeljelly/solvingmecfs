"""Build per-study lead-metabolite matrices (log2, samples x canonical features) + covariates.
Output: work/studies.pkl  {study: dict(X, C, tissue, platform, n_case, n_ctrl, cov_used, note)}
        name_mapping.csv (canonical feature -> exact platform name(s) used, with IDs where the source gives them)
"""
import re, sys, pickle, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, '/mnt/project-files/mecfs-research/downloads/mwb/')
from mwparse import ssf, data, RAW

OUT = '/mnt/project-files/mecfs-research/downloads/mwb/'
BM = '/mnt/project-files/biomap/BioMapAI/'
IC = BM + 'data/indenpendent_cohorts/'

# ---------------------------------------------------------------- canonical lead features: candidate exact names (case/space-insensitive)
LEADS = {
    'IPA (indolepropionate)': ['indolepropionate', 'Indole-3-propionic acid', 'indole_3_propionic_acid_PM_a1'],
    'indolelactate': ['indolelactate', 'indole-3-lactate', 'Indole-3-lactic acid', 'DL-Indole-3-lactic acid', 'indole_3_lactate_PM_a1'],
    'indoleacetate': ['indoleacetate', 'indole-3-acetate', 'Indole-3-acetic acid', 'indole_3_acetate_PM_a1'],
    'indoxyl sulfate': ['3-indoxyl sulfate', 'Indoxyl sulfate'],
    'tryptophan': ['tryptophan', 'L-Tryptophan', 'Tryptophan', 'tryptophan_PM_a1'],
    'kynurenine': ['kynurenine', 'L-Kynurenine', 'Kynurenine'],
    'kynurenate': ['kynurenate', 'Kynurenic acid'],
    'quinolinate': ['quinolinate', 'Quinolinic Acid'],
    'serotonin': ['serotonin', 'Serotonin'],
    '3-hydroxybutyrate': ['3-hydroxybutyrate (BHBA)', '3-hydroxybutyric acid', '3_hydroxybutyric_acid_PM_a1', '3-hydroxybutyrate'],
    'acetoacetate': ['acetoacetate', 'Acetoacetic acid'],
    'palmitate (16:0)': ['palmitate (16:0)', 'palmitic acid', 'palmitic_acid_PM_a1', 'FA (16:0) (palmitic acid)', 'palmitate 160'],
    'oleate (18:1)': ['oleate/vaccenate (18:1)', 'oleic acid', 'FA (18:1) (oleic acid)', 'oleate/vaccenate 181'],
    'linoleate (18:2)': ['linoleate (18:2n6)', 'linoleic acid', 'FA (18:2) (linoleic acid)', 'linoleate 182n6'],
    'carnitine': ['carnitine', 'L-Carnitine', 'L-carnitine', 'Carnitine'],
    'acetylcarnitine': ['acetylcarnitine (C2)', 'L-Acetylcarnitine', 'Acetylcarnitine', 'Acetyl-DL-carnitine', 'acetylcarnitine C2'],
    'palmitoylcarnitine (C16)': ['palmitoylcarnitine (C16)', 'L-Palmitoylcarnitine', 'Palmityl-L-carnitine', 'Acylcarnitine C16:0', 'palmitoylcarnitine C16'],
    'VMA (vanillylmandelate)': ['vanillylmandelate (VMA)', 'vanillylmandelate VMA'],
    'HVA (homovanillate)': ['homovanillate (HVA)', 'Homovanillic acid', 'homovanillate HVA'],
    'glycocholate': ['glycocholate', 'Glycocholic acid'],
    'cholate': ['cholate', 'Cholic acid'],
    'chenodeoxycholate': ['chenodeoxycholate', 'Chenodeoxycholic acid'],
    'deoxycholate': ['deoxycholate', 'Deoxycholic acid'],
    'glycodeoxycholate': ['glycodeoxycholate', 'Glycodeoxycholic acid'],
    'caffeine': ['caffeine', 'Caffeine'],
    'paraxanthine': ['paraxanthine'],
    'hippurate': ['hippurate', 'Hippuric acid', 'hippuric acid'],
}
norm = lambda s: re.sub(r'\s+', ' ', str(s)).strip().lower()

# Metabolon sub-pathway composites (same definitions as the BioMapAI deep analysis)
MB_COMP = {
    'COMP long-chain free fatty acids': lambda s: s.isin(['Long Chain Polyunsaturated Fatty Acid (n3 and n6)', 'Long Chain Saturated Fatty Acid', 'Long Chain Monounsaturated Fatty Acid']),
    'COMP long-chain acylcarnitines': lambda s: s.isin(['Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)']),
    'COMP ketone bodies': lambda s: s.isin(['Ketone Bodies']),
    'COMP ceramides': lambda s: s.isin(['Ceramides']),
    'COMP sphingomyelins': lambda s: s.isin(['Sphingomyelins']),
    'COMP phosphatidylcholines': lambda s: s.isin(['Phosphatidylcholine (PC)']),
    'COMP primary bile acids': lambda s: s.isin(['Primary Bile Acid Metabolism']),
    'COMP secondary bile acids': lambda s: s.str.startswith('Secondary Bile Acid'),
}
# name-regex composites for non-Metabolon platforms
FFA_RE = re.compile(r'^(myristic|palmitic|palmitoleic|stearic|oleic|linoleic|linolenic|arachidonic|eicosapentaenoic|docosahexaenoic)[ _]acid|^FA \((1[4-9]|2[0-2]):\d\)|^(Oleic|Linoleic|Arachidonic|Palmitic|Stearic) acid', re.I)
LCAC_RE = re.compile(r'(myristoyl|palmitoyl|palmityl|hexadecenoyl|octadecanoyl|stearoyl|oleoyl|linoleoyl|linoleyl|arachidonyl|tetradecadien|hydroxyhexadecano|hydroxyhexadeceno|hydroxyoleoyl|hydroxylinoleyl|hydroxyoleyl|tetradecenoyl).{0,6}carnitine|Acylcarnitine C1[4-8]|carnitine.*C1[4-8]', re.I)
CER_RE = re.compile(r'^(Cer\(|Ceramide ?\()', re.I)
SM_RE = re.compile(r'^SM ?\(', re.I)
PC_RE = re.compile(r'^PC ?\(', re.I)
PBA_RE = re.compile(r'^(glyco|tauro)?(cholic|chenodeoxycholic|chenodesoxycholic)[ _]acid$|^Chenodeoxyglycocholic acid$|^Taurochenodesoxycholic acid$', re.I)
SBA_RE = re.compile(r'^(glyco|tauro)?(deoxycholic|lithocholic|ursodeoxycholic)[ _]acid$', re.I)
KET_RE = re.compile(r'^(3-hydroxybutyric acid|3_hydroxybutyric_acid_PM_a1|Acetoacetic acid|3-hydroxybutyrate)$', re.I)
RE_COMP = {'COMP long-chain free fatty acids': FFA_RE, 'COMP long-chain acylcarnitines': LCAC_RE, 'COMP ketone bodies': KET_RE, 'COMP ceramides': CER_RE,
           'COMP sphingomyelins': SM_RE, 'COMP phosphatidylcholines': PC_RE, 'COMP primary bile acids': PBA_RE, 'COMP secondary bile acids': SBA_RE}

MAP = []  # mapping rows


def zs(X): return (X - X.mean()) / X.std()


def build(study, V, ann_ids=None, sub=None, drop_istd=True):
    """V: raw-scale samples x metabolites (positive values; NaN = missing). Returns log2 X with canonical leads + composites."""
    V = V.loc[:, ~V.columns.duplicated()]
    V = V.loc[:, [c for c in V.columns if not (drop_istd and 'istd' in c.lower())]]
    V = V.where(V > 0)
    keep = V.notna().mean() >= 0.5
    V = V.loc[:, keep]
    L = np.log2(V.fillna(V.min() / 2))       # half-minimum imputation
    L = L.loc[:, L.std() > 0]
    lut = {}
    for c in L.columns: lut.setdefault(norm(c), c)
    out = {}
    for f, cands in LEADS.items():
        for c in cands:
            if norm(c) in lut:
                col = lut[norm(c)]; out[f] = L[col]
                ids = ann_ids.get(col, '') if ann_ids is not None else ''
                MAP.append(dict(study=study, feature=f, platform_name=col, ids=ids)); break
    X = pd.DataFrame(out, index=L.index)
    if 'kynurenine' in X and 'tryptophan' in X: X['kynurenine/tryptophan ratio'] = X['kynurenine'] - X['tryptophan']
    Z = zs(L)
    for k in RE_COMP:
        if sub is not None:
            s = sub.reindex(L.columns).fillna('')
            m = MB_COMP[k](s).values
        else:
            m = np.array([bool(RE_COMP[k].search(str(c).strip())) for c in L.columns])
        if m.sum() >= 2 or (k == 'COMP ketone bodies' and m.sum() >= 1):
            X[k] = Z.loc[:, m].mean(axis=1)
            MAP.append(dict(study=study, feature=k, platform_name='; '.join(map(str, L.columns[m])), ids=f'n={int(m.sum())}'))
    return X


S = {}


def add(study, X, C, tissue, platform, cov_used, note='', person=None):
    C = C.loc[X.index]
    S[study] = dict(X=X, C=C, tissue=tissue, platform=platform, cov_used=cov_used, note=note,
                    n_case=int(C.case.sum()), n_ctrl=int((1 - C.case).sum()))
    print(f'{study:32s} {tissue:6s} X{X.shape} case {int(C.case.sum())} ctrl {int((1-C.case).sum())} feats: {len(X.columns)}')


# ================================================================ 1. BioMapAI (Xiong 2025): QC-normalised Metabolon, all visits
sys.path.insert(0, '/mnt/project-files/biomap/deep/stats_replication/')
from data import load_meta
q = pd.read_csv(BM + 'data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv')
q = q[q.TYPE == 'NAMED']; q.index = q.CHEMICAL_NAME
Sm = [c for c in q.columns if c.startswith('SAM')]
V = q[Sm].astype(float).T
ids = {r.CHEMICAL_NAME: f"HMDB:{r.get('HMDB','')} KEGG:{r.get('KEGG','')} PUBCHEM:{r.get('PUBCHEM','')}" for _, r in q.iterrows()}
X = build('BioMapAI_Xiong2025', V, ids, q.SUB_PATHWAY)
m = load_meta(); idx = X.index.intersection(m.index); X = X.loc[idx]; M = m.loc[idx]
C = pd.DataFrame({'case': M.case, 'age': M.age, 'male': M.male, 'bmi': M.bmi_imp, 'bmi_missing': M.bmi_missing, 'tp2': M.tp2, 'tp3': M.tp3, 'person': M.person}, index=idx)
add('BioMapAI_Xiong2025', X, C, 'plasma', 'Metabolon HD4', 'age, sex, BMI, visit; random intercept per person', 'discovery cohort')

# ================================================================ 2. Che 2022 (Columbia/MCAM plasma, Fiehn West Coast Metabolomics Center)
che = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_metabolomics_yali_v2.csv', index_col=0)
ph = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_pheno_text_binned_v2.txt', sep='\t', index_col=0)
V = che.iloc[:, 3:].apply(pd.to_numeric, errors='coerce').T; V.columns = che.Molecule.astype(str).values
ids = {r.Molecule: f"InChIKey:{r.InChiKey} {r.database_identifier}" for _, r in che.iterrows()}
X = build('Che2022_MCAM', V, ids)
ph = ph.loc[X.index]
agemap = {'18 - 24': 21, '25 - 34': 30, '35 - 44': 40, '45 - 54': 50, '55 - 64': 60, '65 and over': 70}
bmimap = {'underweight': 18, 'healthy weight': 22, 'overweight': 27, 'obese': 33}
C = pd.DataFrame({'case': (ph.Phenotype == 'case').astype(int), 'age': ph.age_binned.map(agemap), 'male': (ph.sex == 'male').astype(int), 'bmi': ph.bmi_binned.map(bmimap)}, index=ph.index)
for s in sorted(ph.site.unique())[1:]: C['site_' + s] = (ph.site == s).astype(int)
add('Che2022_MCAM', X, C, 'plasma', 'Fiehn GC-TOF/HILIC/CSH/oxylipin (MW ST002000-3)', 'age bin, sex, BMI bin, site')

# ================================================================ 3. Germain 2022 JCI Insight (Hanson, 2-day CPET), pre-exercise samples averaged per person
d = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='Original Dataset', header=None)
fm = d.iloc[8:, :15].copy(); fm.columns = d.iloc[7, :15].values
gm = d.iloc[:7, 16:].T.copy(); gm.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'AgeBin', 'BMIBin', 'Bell']; gm = gm.reset_index(drop=True)
V = d.iloc[8:, 16:].apply(pd.to_numeric, errors='coerce').T.reset_index(drop=True); V.columns = fm.BIOCHEMICAL.astype(str).values
pre = gm.Timepoint.astype(str).str.endswith('PRE').values
Vp = V[pre].copy(); gp = gm[pre]
Vp.index = gp.ID.astype(str).values
Vp = np.exp2(np.log2(Vp.where(Vp > 0)).groupby(level=0).mean())   # geometric mean of D1-PRE and D2-PRE
sub = pd.Series(fm['SUB PATHWAY'].astype(str).values, index=fm.BIOCHEMICAL.astype(str).values); sub = sub[~sub.index.duplicated()]
ids = {r.BIOCHEMICAL: f"HMDB:{r.get('HMDB', '')} KEGG:{r.get('KEGG', '')} PUBCHEM:{r.get('PUBCHEM', '')}" for _, r in fm.iterrows()}
X = build('Germain2022_JCI', Vp, ids, sub)
g1 = gp.drop_duplicates('ID').set_index(gp.drop_duplicates('ID').ID.astype(str)).loc[X.index]
C = pd.DataFrame({'case': (g1.Phenotype == 'CFS').astype(int), 'age': g1.AgeBin.astype(float), 'male': (g1.Sex == 'M').astype(int), 'bmi': g1.BMIBin.astype(float)}, index=X.index)
add('Germain2022_JCI', X, C, 'plasma', 'Metabolon HD4', 'age bin, sex, BMI bin', 'mean of D1-PRE and D2-PRE per person')

# ================================================================ 4/5. Germain 2018 & 2020 Metabolites (Hanson; Metabolon; per-sample supplementary)
for study, f, sheetcol in [('Germain2018_Metabolites', RAW + 'Germain2018_Metabolites/Supplementary_File_1.xlsx', 'SUB PATHWAY'),
                           ('Germain2020_Metabolites', RAW + 'Germain2020_Metabolites/Supplementary File 1.xlsx', 'SUB_PATHWAY')]:
    g = pd.read_excel(f, header=0)
    g.columns = [str(c).strip() for c in g.columns]
    samp = [c for c in g.columns if re.fullmatch(r'[CP]\d+', c)]
    V = g[samp].apply(pd.to_numeric, errors='coerce').T; V.columns = g.BIOCHEMICAL.astype(str).values
    hm = [c for c in g.columns if 'HMDB' in c][0]
    ids = {r['BIOCHEMICAL']: f"HMDB:{r[hm]} KEGG:{r.get('KEGG', '')} PUBCHEM:{r.get('PUBCHEM', '')}" for _, r in g.iterrows()}
    sub = pd.Series(g[sheetcol].astype(str).values, index=g.BIOCHEMICAL.astype(str).values); sub = sub[~sub.index.duplicated()]
    X = build(study, V, ids, sub)
    C = pd.DataFrame({'case': [int(s.startswith('P')) for s in X.index]}, index=X.index)
    add(study, X, C, 'plasma', 'Metabolon', 'none (no per-sample covariates published)', 'may share participants with other Hanson-lab cohorts')

# ================================================================ 6. Naviaux 2016 PNAS (MW ST000450), targeted LC-MS/MS
Xs = []
for a in ['AN000705', 'AN000706']:
    V, an = data(RAW + a + '.mwtab.txt'); Xs.append(V)
    ids = {r.metabolite_name: f"KEGG:{r['KEGG ID']} PUBCHEM:{r['PubChem ID']}" for _, r in an.iterrows()}
V = pd.concat(Xs, axis=1)
ids = {}
for a in ['AN000705', 'AN000706']:
    _, an = data(RAW + a + '.mwtab.txt'); ids.update({r.metabolite_name: f"KEGG:{r['KEGG ID']} PUBCHEM:{r['PubChem ID']}" for _, r in an.iterrows()})
X = build('Naviaux2016_PNAS', V, ids)
s = ssf(RAW + 'AN000705.mwtab.txt').set_index('sample').loc[X.index]
C = pd.DataFrame({'case': (s.Disease == 'CFS').astype(int), 'male': (s.gender == 'Male').astype(int)}, index=X.index)
add('Naviaux2016_PNAS', X, C, 'plasma', 'targeted LC-MS/MS (Naviaux lab)', 'sex')

# ================================================================ 7. Nagy-Szakal 2018 Sci Rep (MW ST000800/801/802/813; Fiehn)
s = ssf(RAW + 'AN001271.mwtab.txt'); s = s[s.CaseStatus.isin(['Case', 'Control'])]
num2esp = {int(re.search(r'ESP-\d+_(\d+)', r.subject).group(1)): r['sample'] for _, r in s.iterrows()}
blocks = []
for f, namecol in [('ST000800_AN001271_Results.txt', 'Annotation'), ('ST000800_AN001272_Results.txt', 'Annotation')]:
    t = pd.read_csv(RAW + f, sep='\t', low_memory=False)
    t = t[t[namecol].notna()]
    cols = [c for c in t.columns if re.fullmatch(r'\d+', str(c))]
    V = t[cols].apply(pd.to_numeric, errors='coerce').T; V.columns = t[namecol].astype(str).str.strip().values
    V.index = [num2esp.get(int(c)) for c in cols]; V = V[V.index.notna()]
    blocks.append(V)
t = pd.read_csv(RAW + 'ST000801_AN001273_Results.txt', sep='\t', low_memory=False); t = t[t.Annotation.notna()]
cols = [c for c in t.columns if 'ESP-' in c]
V = t[cols].apply(pd.to_numeric, errors='coerce').T; V.columns = t.Annotation.astype(str).str.strip().values
V.index = [re.search(r'(ESP-\d+)', c).group(1) for c in cols]; V = V.groupby(level=0).mean(); blocks.append(V)
t = pd.read_csv(RAW + 'ST000813_AN001286_Results.txt', sep='\t', low_memory=False)
t = t[t[' Metabolite Name'].astype(str).str.contains('[a-z]', regex=True) & (t[' Metabolite Name'] != 'Factors')]
cols = [c for c in t.columns if c.startswith('ESP-')]
V = t[cols].apply(pd.to_numeric, errors='coerce').T; V.columns = t[' Metabolite Name'].astype(str).str.strip().values
V.index = [re.search(r'(ESP-\d+)', c).group(1) for c in cols]; blocks.append(V)
# keep GC-TOF names first for FA/3-HB/indoles (GC is the validated platform for these), then HILIC, then CSH
V = pd.concat([blocks[3], blocks[2], blocks[1], blocks[0]], axis=1, join='inner')
X = build('NagySzakal2018_SciRep', V)
s = s.set_index('sample').loc[X.index]
C = pd.DataFrame({'case': (s.CaseStatus == 'Case').astype(int), 'age_over50': (s.AGE_DICH == 'over50').astype(int), 'male': (s.SEX == 'Male').astype(int),
                  'bmi': pd.to_numeric(s.x_bmi, errors='coerce')}, index=X.index)
C['bmi'] = C.bmi.fillna(C.bmi.median())
add('NagySzakal2018_SciRep', X, C, 'plasma', 'Fiehn GC-TOF + HILIC + CSH lipidomics', 'age (<=50 vs >50), sex, BMI')

# ================================================================ 8. NCNP Japan (MW ST004940) serum tryptophan/indole panel
V, an = data(RAW + 'AN008373.mwtab.txt')
ids = {r.metabolite_name: f"KEGG:{r['KEGG ID']} PUBCHEM:{r['PubChem ID']}" for _, r in an.iterrows()}
X = build('NCNP_Japan_serum_ST004940', V, ids)
s = ssf(RAW + 'AN008373.mwtab.txt').set_index('sample').loc[X.index]
C = pd.DataFrame({'case': (s.Condition == 'ME/CFS').astype(int)}, index=X.index)
add('NCNP_Japan_serum_ST004940', X, C, 'serum', 'targeted LC-MS/MS tryptophan panel', 'none published per sample')

# ================================================================ 9. Walitt 2024 NIH plasma lipid classes (Supplementary Data 20)
w = pd.read_excel(RAW + 'Walitt2024_NatCommun/Supplementary Data 20.xlsx', header=None)
grp = w.iloc[2, 1:-1].values; V = w.iloc[3:, 1:-1].astype(float).T; V.columns = w.iloc[3:, 0].astype(str).values
V.index = [f'W{i}' for i in range(len(V))]
L = np.log2(V); Z = zs(L)
X = pd.DataFrame({'COMP long-chain free fatty acids': Z['FFA'], 'COMP ceramides': Z['CER'], 'COMP sphingomyelins': Z['SM'], 'COMP phosphatidylcholines': Z['PC']}, index=V.index)
for k, v in [('COMP long-chain free fatty acids', 'FFA (total, Lipidyzer)'), ('COMP ceramides', 'CER (total)'), ('COMP sphingomyelins', 'SM (total)'), ('COMP phosphatidylcholines', 'PC (total)')]:
    MAP.append(dict(study='Walitt2024_NIH_plasma_lipids', feature=k, platform_name=v, ids='lipid class sum'))
C = pd.DataFrame({'case': (grp == 'PI-ME/CFS').astype(int)}, index=V.index)
add('Walitt2024_NIH_plasma_lipids', X, C, 'plasma', 'Lipidyzer lipid-class sums', 'none (per-sample covariates not linked)')

# ================================================================ 10. Armstrong 2015 (MetaboLights MTBLS161) serum NMR
mf = pd.read_csv(RAW + 'MTBLS161/m_MTBLS161_mecfs_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv', sep='\t')
sc = [c for c in mf.columns if c.endswith('serum')]
V = mf[sc].apply(lambda s: pd.to_numeric(s.astype(str).str.replace(',', ''), errors='coerce')).T; V.columns = mf.metabolite_identification.astype(str).values
ids = dict(zip(mf.metabolite_identification.astype(str), mf.database_identifier.astype(str)))
V.columns = [c.replace('L-Tryptophan', 'tryptophan').replace('Hippurate', 'hippurate') for c in V.columns]
X = build('Armstrong2015_MTBLS161', V, {k.replace('L-Tryptophan', 'tryptophan'): v for k, v in ids.items()})
C = pd.DataFrame({'case': [int(i.startswith('CFS')) for i in X.index]}, index=X.index)
add('Armstrong2015_MTBLS161', X, C, 'serum', '1H-NMR', 'none')

# ================================================================ CSF: NIH Walitt 2024 (MW ST003178), Metabolon
Vs, anns = [], []
for a in ['AN005217', 'AN005218', 'AN005219', 'AN005220']:
    V, an = data(RAW + a + '.mwtab.txt'); Vs.append(V); anns.append(an)
V = pd.concat(Vs, axis=1); an = pd.concat(anns)
ids = {r.metabolite_name: f"HMDB:{r.HMDB} KEGG:{r.KEGG} PUBCHEM:{r.PUBCHEM}" for _, r in an.iterrows()}
sub = pd.Series(an.SUB_PATHWAY.values, index=an.metabolite_name.values); sub = sub[~sub.index.duplicated()]
X = build('Walitt2024_NIH_CSF', V, ids, sub)
s = ssf(RAW + 'AN005217.mwtab.txt').set_index('sample').loc[X.index]
C = pd.DataFrame({'case': (s.MECFS_STATUS == 'MECFS').astype(int), 'age': pd.to_numeric(s.x_AGE, errors='coerce'), 'male': (s.SEX == 'M').astype(int)}, index=X.index)
add('Walitt2024_NIH_CSF', X, C, 'CSF', 'Metabolon', 'age, sex')

# ================================================================ CSF: ST000910-913 (Columbia/Fiehn CSF), ME/CFS vs no-neurological-disease controls
Vs, ids = [], {}
for a in ['AN001480', 'AN001481', 'AN001482']:
    V, an = data(RAW + a + '.mwtab.txt')
    sm = ssf(RAW + a + '.mwtab.txt').set_index('sample').subject
    V.index = [sm.get(i, i) for i in V.index]; V = V.groupby(level=0).mean(); Vs.append(V)
    ids.update({r.metabolite_name: f"KEGG:{r['KEGG ID']} PUBCHEM:{r['PubChem ID']} InChIKey:{r['InChI Key']}" for _, r in an.iterrows()})
t = pd.read_csv(RAW + 'ST000913_AN001483_Results.txt', sep='\t', low_memory=False)
V4 = None
V = pd.concat(Vs, axis=1, join='inner')
s = ssf(RAW + 'AN001480.mwtab.txt').drop_duplicates('subject').set_index('subject')
s = s.loc[V.index]; keep = s.Diagnosis.isin(['MECFS', 'ND']).values
X = build('CSF2018_ST000910', V[keep], ids)
s = s.loc[X.index]
C = pd.DataFrame({'case': (s.Diagnosis == 'MECFS').astype(int), 'age': pd.to_numeric(s.x_Age, errors='coerce'), 'male': (s.Sex.str.upper() == 'MALE').astype(int)}, index=X.index)
C['age'] = C.age.fillna(C.age.median())
add('CSF2018_ST000910', X, C, 'CSF', 'Fiehn GC-TOF + CSH', 'age, sex', 'controls = no neurological disease')

# ================================================================ CSF: Baraniuk 2025 Sci Rep (per-sample supplementary, z-scored values, with/without prior exercise)
b = pd.read_excel(RAW + 'Baraniuk2025_CSF/41598_2025_91324_MOESM2_ESM.xlsx', header=0)
b = b[b.Disease.isin(['cfs', 'sc'])].reset_index(drop=True)
bm = {'tryptophan': 'TRYPTOPHAN_pos_1', '3-hydroxybutyrate': '3-HYDROXYBUTANOATE_pos_1', 'carnitine': 'L-CARNITINE_pos_1', 'acetylcarnitine': 'L-ACETYLCARNITINE_pos_1'}
X = pd.DataFrame({k: pd.to_numeric(b[v], errors='coerce') for k, v in bm.items()})
for k, v in bm.items(): MAP.append(dict(study='Baraniuk2025_CSF', feature=k, platform_name=v, ids='supplied as standardised values'))
X.index = [f'B{i}' for i in range(len(b))]
C = pd.DataFrame({'case': (b.Disease == 'cfs').astype(int).values, 'age': pd.to_numeric(b.Age_jnb, errors='coerce').values,
                  'male': (b.Male.astype(str).str.lower() == 'male').astype(int).values, 'postexercise': (b.Exercise == 'postexercise').astype(int).values}, index=X.index)
C['age'] = C.age.fillna(C.age.median())
add('Baraniuk2025_CSF', X, C, 'CSF', 'untargeted LC-MS (values already standardised)', 'age, sex, exercise group', 'values supplied already scaled; not log-transformed again')

# ================================================================ Stool: ST003346 (UConn fast targeted indoles), ST004941 (NCNP fecal)
V, an = data(RAW + 'AN005481.mwtab.txt'); s = ssf(RAW + 'AN005481.mwtab.txt').set_index('sample').loc[V.index]
V = np.exp2(np.log2(V.where(V > 0)).groupby(s.subject.values).mean())
X = build('UConn_stool_ST003346', V)
cs = s.drop_duplicates('subject').set_index('subject').loc[X.index]
C = pd.DataFrame({'case': (cs.Condition == 'MECFS').astype(int), 'male': (cs.x_Sex == 'Male').astype(int)}, index=X.index)
add('UConn_stool_ST003346', X, C, 'stool', 'targeted LC-MS indoles', 'sex', 'triplicates and repeat aliquots averaged per person')

V, an = data(RAW + 'AN008374.mwtab.txt'); ids = {r.metabolite_name: f"KEGG:{r['KEGG ID']} PUBCHEM:{r['PubChem ID']}" for _, r in an.iterrows()}
X = build('NCNP_Japan_stool_ST004941', V, ids)
s = ssf(RAW + 'AN008374.mwtab.txt').set_index('sample').loc[X.index]
C = pd.DataFrame({'case': (s.Condition == 'ME/CFS').astype(int)}, index=X.index)
add('NCNP_Japan_stool_ST004941', X, C, 'stool', 'targeted LC-MS/MS tryptophan panel', 'none')

pickle.dump(S, open(OUT + 'work_studies.pkl', 'wb'))
M = pd.DataFrame(MAP).drop_duplicates()
M.to_csv(OUT + 'name_mapping.csv', index=False)
print(M.pivot_table(index='feature', columns='study', values='platform_name', aggfunc='first').notna().sum(axis=1))
