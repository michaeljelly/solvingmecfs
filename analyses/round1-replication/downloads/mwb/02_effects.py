"""Per-study patient-vs-control effects (SD units, covariate-adjusted) and random-effects meta-analysis.
Outputs: per_study_effects.csv, meta_blood.csv, meta_csf_stool.csv, fig_forest_leads.png
"""
import sys, pickle, re, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
warnings.filterwarnings('ignore')
OUT = '/mnt/project-files/mecfs-research/downloads/mwb/'
RAW = OUT + 'raw/'
sys.path.insert(0, '/mnt/project-files/biomap/deep/stats_replication/')
from lmmfast import fit_many

S = pickle.load(open(OUT + 'work_studies.pkl', 'rb'))


def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = ~np.isnan(p)
    if ok.sum():
        pp = p[ok]; o = np.argsort(pp); r = pp[o] * ok.sum() / (np.arange(ok.sum()) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]; qq = np.empty_like(pp); qq[o] = np.minimum(r, 1); q[ok] = qq
    return q


rows = []
for st, d in S.items():
    X, C = d['X'], d['C'].copy()
    X = (X - X.mean()) / X.std()
    covs = [c for c in C.columns if c not in ('case', 'person')]
    for c in covs:
        if C[c].nunique() > 2: C[c] = (C[c] - C[c].mean()) / C[c].std()
    if st == 'BioMapAI_Xiong2025':
        D = np.column_stack([np.ones(len(C)), C.case] + [C[c] for c in covs]).astype(float)
        Y = X.fillna(X.mean())
        r = fit_many(Y.values, D, C.person.values, L=np.eye(D.shape[1])[[1]])
        npc = C.drop_duplicates('person'); nc, nh = int(npc.case.sum()), int((1 - npc.case).sum())
        for k, f in enumerate(X.columns):
            rows.append(dict(study=st, tissue=d['tissue'], feature=f, est=r['est'][k, 0], se=r['se'][k, 0], p=r['p'][k, 0], n_case=nc, n_ctrl=nh, adjusted=d['cov_used']))
        continue
    for f in X.columns:
        y = X[f]; ok = y.notna()
        D = sm.add_constant(C.loc[ok, ['case'] + covs].astype(float), has_constant='add')
        fit = sm.OLS(y[ok].astype(float), D).fit()
        rows.append(dict(study=st, tissue=d['tissue'], feature=f, est=fit.params['case'], se=fit.bse['case'], p=fit.pvalues['case'],
                         n_case=int(C.loc[ok, 'case'].sum()), n_ctrl=int((1 - C.loc[ok, 'case']).sum()), adjusted=d['cov_used']))

# ---------------------------------------------------------------- Hoel 2021 (summary level: Welch p + mean fold change, 83 ME/CFS vs 35 HC)
h = pd.read_excel(RAW + 'Hoel2021/sdd1.xlsx', '3.FC & stats', header=1)
h = h.rename(columns={'MEall': 'FC', 'MEall.1': 'p'})[['Met', 'SUB_PATHWAY', 'HMDB_ID', 'FC', 'p']]
x = pd.read_excel(RAW + 'Hoel2021/sdd1.xlsx', '11.HD4 Xenobiotics', header=2).rename(columns={'Biochemical': 'Met', 'Sub_pathway': 'SUB_PATHWAY', 'MEall': 'FC', 'MEall.1': 'p'})[['Met', 'SUB_PATHWAY', 'FC', 'p']]
h = pd.concat([h, x]).dropna(subset=['FC', 'p'])
n1, n0 = 83, 35
vd = 1 / n1 + 1 / n0


def d_from(fc, p):
    z = stats.norm.isf(np.clip(p, 1e-300, 1) / 2) * np.sign(np.log(fc))
    dd = z * np.sqrt(vd); return dd, np.sqrt(vd + dd ** 2 / (2 * (n1 + n0)))


HOEL = {'IPA (indolepropionate)': 'indolepropionate', 'indolelactate': 'indolelactate', 'indoleacetate': 'indoleacetate', 'indoxyl sulfate': '3-indoxyl sulfate',
        'tryptophan': 'tryptophan', 'kynurenine': 'kynurenine', 'kynurenate': 'kynurenate', 'quinolinate': 'quinolinate', 'serotonin': 'serotonin',
        '3-hydroxybutyrate': '3-hydroxybutyrate', 'acetoacetate': 'acetoacetate', 'palmitate (16:0)': 'palmitate 160', 'oleate (18:1)': 'oleate/vaccenate 181',
        'linoleate (18:2)': 'linoleate 182n6', 'carnitine': 'carnitine', 'acetylcarnitine': 'acetylcarnitine C2', 'palmitoylcarnitine (C16)': 'palmitoylcarnitine C16',
        'VMA (vanillylmandelate)': 'vanillylmandelate VMA', 'HVA (homovanillate)': 'homovanillate HVA', 'glycocholate': 'glycocholate', 'cholate': 'cholate',
        'chenodeoxycholate': 'chenodeoxycholate', 'deoxycholate': 'deoxycholate', 'glycodeoxycholate': 'glycodeoxycholate', 'caffeine': 'caffeine',
        'paraxanthine': 'paraxanthine', 'hippurate': 'hippurate'}
maprows = []
for f, nm in HOEL.items():
    r = h[h.Met.astype(str).str.strip() == nm]
    if len(r):
        r = r.iloc[0]; dd, se = d_from(r.FC, r.p)
        rows.append(dict(study='Hoel2021_JCI_summary', tissue='serum', feature=f, est=dd, se=se, p=r.p, n_case=n1, n_ctrl=n0, adjusted='none (from published Welch p and fold change)'))
        maprows.append(dict(study='Hoel2021_JCI_summary', feature=f, platform_name=nm, ids=f"HMDB:{r.get('HMDB_ID', '')}"))
# composites for Hoel: mean d of members (SE = mean member SE; conservative, assumes full correlation)
HSUB = {'COMP long-chain acylcarnitines': ['Fatty Acid Metabolism (Acyl Carnitine, Long Chain Saturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Monounsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Polyunsaturated)', 'Fatty Acid Metabolism (Acyl Carnitine, Hydroxy)'],
        'COMP ketone bodies': ['Ketone Bodies'], 'COMP primary bile acids': ['Primary Bile Acid Metabolism'], 'COMP secondary bile acids': ['Secondary Bile Acid Metabolism']}
for f, subs in HSUB.items():
    r = h[h.SUB_PATHWAY.isin(subs)]
    dd, se = d_from(r.FC.values.astype(float), r.p.values.astype(float))
    m, s = dd.mean(), se.mean()
    rows.append(dict(study='Hoel2021_JCI_summary', tissue='serum', feature=f, est=m, se=s, p=2 * stats.norm.sf(abs(m / s)), n_case=n1, n_ctrl=n0, adjusted='none; composite = mean of member effects'))
    maprows.append(dict(study='Hoel2021_JCI_summary', feature=f, platform_name='; '.join(r.Met.astype(str)), ids=f'n={len(r)}'))
cl = pd.read_excel(RAW + 'Hoel2021/sdd2.xlsx', '3.FC CLP classes', header=None).iloc[4:, [2, 3, 8]]
cl.columns = ['cls', 'FC', 'p']
for f, c in [('COMP long-chain free fatty acids', 'FFA'), ('COMP ceramides', 'CER'), ('COMP sphingomyelins', 'SM'), ('COMP phosphatidylcholines', 'PC')]:
    r = cl[cl.cls == c].iloc[0]; dd, se = d_from(float(r.FC), float(r.p))
    rows.append(dict(study='Hoel2021_JCI_summary', tissue='serum', feature=f, est=dd, se=se, p=float(r.p), n_case=n1, n_ctrl=n0, adjusted='none; Lipidyzer class sum'))
    maprows.append(dict(study='Hoel2021_JCI_summary', feature=f, platform_name=f'{c} (Lipidyzer class sum)', ids='class'))
mp = pd.read_csv(OUT + 'name_mapping.csv'); pd.concat([mp, pd.DataFrame(maprows)]).drop_duplicates().to_csv(OUT + 'name_mapping.csv', index=False)

E = pd.DataFrame(rows)
E['q_within_study'] = E.groupby('study').p.transform(bh)
E.to_csv(OUT + 'per_study_effects.csv', index=False)


# ---------------------------------------------------------------- meta-analysis
def meta(e, s):
    e, s = np.asarray(e, float), np.asarray(s, float)
    w = 1 / s ** 2; fe = np.sum(w * e) / w.sum(); k = len(e); Q = np.sum(w * (e - fe) ** 2)
    tau2 = max(0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0
    wr = 1 / (s ** 2 + tau2); re_ = np.sum(wr * e) / wr.sum(); rse = np.sqrt(1 / wr.sum())
    I2 = max(0, (Q - (k - 1)) / Q) if (k > 1 and Q > 0) else 0
    return re_, rse, 2 * stats.norm.sf(abs(re_ / rse)), I2, tau2


BLOOD_ALL = [s for s in E.study.unique() if E[E.study == s].tissue.iloc[0] in ('plasma', 'serum')]
out = []
for f in E.feature.unique():
    b = E[(E.feature == f) & E.study.isin(BLOOD_ALL)]
    r = dict(feature=f)
    bmrow = b[b.study == 'BioMapAI_Xiong2025']
    r['BioMapAI_est'] = bmrow.est.iloc[0] if len(bmrow) else np.nan
    r['BioMapAI_p'] = bmrow.p.iloc[0] if len(bmrow) else np.nan
    for lab, sub in [('ext', b[b.study != 'BioMapAI_Xiong2025']), ('all', b),
                     ('ext_noG1820', b[~b.study.isin(['BioMapAI_Xiong2025', 'Germain2018_Metabolites', 'Germain2020_Metabolites'])])]:
        if len(sub) >= 2:
            m, s, p, I2, t2 = meta(sub.est, sub.se)
            r.update({f'{lab}_k': len(sub), f'{lab}_est': m, f'{lab}_lo': m - 1.96 * s, f'{lab}_hi': m + 1.96 * s, f'{lab}_p': p, f'{lab}_I2': I2,
                      f'{lab}_n_neg': int((sub.est < 0).sum()), f'{lab}_n_sig': int((sub.p < 0.05).sum()), f'{lab}_N': int(sub.n_case.sum() + sub.n_ctrl.sum())})
        elif len(sub) == 1:
            r.update({f'{lab}_k': 1, f'{lab}_est': sub.est.iloc[0], f'{lab}_p': sub.p.iloc[0]})
    out.append(r)
Mb = pd.DataFrame(out)
Mb['ext_q'] = bh(Mb.ext_p.where(Mb.ext_k >= 2))
Mb['all_q'] = bh(Mb.all_p.where(Mb.all_k >= 2))


def verdict(r):
    if not (r.get('ext_k', 0) >= 2): return 'too few external cohorts'
    disc = (not np.isnan(r.BioMapAI_p)) and r.BioMapAI_p < 0.05
    same = (not np.isnan(r.BioMapAI_est)) and np.sign(r.ext_est) == np.sign(r.BioMapAI_est)
    if r.ext_q < 0.05 and disc and same: return 'REPLICATES: BioMapAI p<0.05 and external q<0.05, same direction'
    if r.ext_q < 0.05 and disc and not same: return 'external q<0.05 but OPPOSITE to BioMapAI'
    if r.ext_q < 0.05: return 'external signal (q<0.05) not seen in BioMapAI'
    if r.ext_p < 0.05 and disc and same: return 'partial: same direction, external p<0.05 but q>=0.05'
    if r.ext_p < 0.05: return 'nominal external signal only (p<0.05, q>=0.05)'
    if disc: return 'BioMapAI signal does NOT replicate'
    return 'null in all'


Mb['verdict'] = Mb.apply(verdict, axis=1)
Mb = Mb.sort_values('ext_p')
Mb.to_csv(OUT + 'meta_blood.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100)
print(Mb[['feature', 'BioMapAI_est', 'BioMapAI_p', 'ext_k', 'ext_N', 'ext_est', 'ext_lo', 'ext_hi', 'ext_p', 'ext_q', 'ext_I2', 'ext_n_neg', 'all_est', 'all_p', 'ext_noG1820_est', 'ext_noG1820_p', 'verdict']].round(3).to_string())

# CSF and stool: per-study only + small meta within tissue
o2 = []
for t in ['CSF', 'stool']:
    for f in E.feature.unique():
        b = E[(E.feature == f) & (E.tissue == t)]
        if len(b) == 0: continue
        r = dict(tissue=t, feature=f, k=len(b), studies='; '.join(f'{s}:{e:+.2f} (p={p:.3g})' for s, e, p in zip(b.study, b.est, b.p)))
        if len(b) >= 2:
            m, s, p, I2, _ = meta(b.est, b.se); r.update(est=m, lo=m - 1.96 * s, hi=m + 1.96 * s, p=p, I2=I2)
        else: r.update(est=b.est.iloc[0], p=b.p.iloc[0])
        o2.append(r)
M2 = pd.DataFrame(o2); M2['q'] = M2.groupby('tissue').p.transform(bh)
M2.to_csv(OUT + 'meta_csf_stool.csv', index=False)
print(M2.round(3).to_string())
