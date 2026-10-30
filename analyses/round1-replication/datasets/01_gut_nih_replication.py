"""Gut replication in a third independent stool cohort: NIH intramural PI-ME/CFS (Walitt 2024, 17 ME/CFS / 15 HV,
JAMS species relative abundance from the authors' GitHub). Tests:
 (a) the 14 pre-specified BioMapAI gut species + butyrate-producer composite (direction fixed from earlier work);
 (b) a pre-specified signature score (mean z of 'up' species minus mean z of 'down' species) -> AUC, no fitting;
 (c) all species shared by NIH, Raijmakers 2020 (NL) and Guo 2023 (Columbia): per-cohort SMD, direction concordance,
     random-effects meta, BH-FDR.
Effect = standardised mean difference of CLR abundance (OLS, NIH adjusted for age+sex, Columbia for age/sex/BMI/site)."""
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests

OUT = '/mnt/project-files/mecfs-research/datasets/results/'
import os; os.makedirs(OUT, exist_ok=True)
NIH = '/mnt/project-files/mecfs-research/datasets/repos/docwalitt_National-Institutes-of-Health-Myalgic-Encephalomyelitis-Chronic-Fatigue-Syndrome-Code-Repository/data/'
IC = '/mnt/project-files/biomap/BioMapAI/data/indenpendent_cohorts/'
REP = '/mnt/project-files/biomap/deep/stats_replication/replication_table.csv'

DOWN = ['Faecalibacterium_prausnitzii', 'Roseburia_faecis', 'Roseburia_intestinalis', 'Roseburia_hominis', 'Roseburia_inulinivorans',
        'Coprococcus_eutactus', 'Coprococcus_comes', 'Coprococcus_catus', 'Anaerostipes_hadrus']
UP = ['Ruminococcus_gnavus', 'Flavonifractor_plautii', 'Dysosmobacter_welbionis', 'Enterocloster_bolteae', 'Clostridium_symbiosum']
BP_GEN = ['Faecalibacterium', 'Roseburia', 'Coprococcus', 'Anaerostipes', 'Agathobacter']


def clr(rel, prev=0.2, keep=()):
    rel = rel.div(rel.sum(axis=1), axis=0)
    k = ((rel > 0).mean() >= prev) | rel.columns.isin(keep)
    r = rel.loc[:, k]; pc = r[r > 0].min().min() / 2
    lg = np.log(r + pc); return lg.sub(lg.mean(axis=1), axis=0)


def bp(rel):
    rel = rel.div(rel.sum(axis=1), axis=0)
    cols = [c for c in rel.columns if any(c.startswith(g) for g in BP_GEN)]
    return np.log10(rel[cols].sum(axis=1) + 1e-4)


def smd(X, case, cov=None):
    X = (X - X.mean()) / X.std()
    C = np.column_stack([np.ones(len(case)), case] + ([] if cov is None else [cov]))
    out = []
    for f in X.columns:
        y = X[f].values
        if np.nanstd(y) == 0: continue
        r = sm.OLS(y, C).fit(); out.append((f, r.params[1], r.bse[1], r.pvalues[1]))
    return pd.DataFrame(out, columns=['feature', 'est', 'se', 'p']).set_index('feature')


# ---------------------------------------------------------------- NIH
x = pd.read_excel(NIH + 'Walitt2023_Microbiome_metadata_and_PPM_relative_abundance.xlsx', sheet_name=None)
md = x['Metadata'].set_index('Sample')
ab = x['LKT_PPM'].set_index('row.names')
ab = ab[ab.index.str.startswith('LKT__s__')]
ab.index = ab.index.str.replace('LKT__s__', '', regex=False)
nih = ab.T.loc[md.index]
case_n = (md.Group == 'ME/CFS').astype(int).values
cov_n = np.column_stack([(md.Age - 40) / 15, (md.Sex == 'Male').astype(float)])
print('NIH n', case_n.sum(), (1 - case_n).sum())

# ---------------------------------------------------------------- other stool cohorts (MetaPhlAn4)
r = pd.read_csv(IC + 'Raijmakers/metaphlan4_species_full.csv', index_col=0).select_dtypes('number').T
r = r[r.index.str.startswith('CFS') | r.index.str.startswith('HV')]
case_r = r.index.str.startswith('CFS').astype(int)
c = pd.read_csv(IC + 'columbia/metaphlan4_species_full.csv', index_col=0).select_dtypes('number').T
sra = pd.read_csv(IC + 'columbia/SraRunTable.txt').set_index('Run').loc[c.index]
case_c = (sra.sampletype == 'Case').astype(int).values
cov_c = np.column_stack([(sra.age_at_a1 - 45) / 10, (sra.sex == 'male').astype(float), ((sra.BMI - 26) / 5).fillna(0),
                         pd.get_dummies(sra.Site, drop_first=True).astype(float).values])

cohorts = {'NIH_Walitt2024': (nih, case_n, cov_n), 'Raijmakers2020': (r, case_r, None), 'Guo2023_Columbia': (c, case_c, cov_c)}
res = {}
for k, (rel, cs, cv) in cohorts.items():
    X = clr(rel, keep=DOWN + UP)
    X['COMPOSITE_butyrate_producers'] = bp(rel)
    res[k] = smd(X, cs, cv)

# (a) pre-specified species in NIH
rep = pd.read_csv(REP).set_index('feature')
rows = []
for f in DOWN + UP + ['COMPOSITE_butyrate_producers']:
    exp = -1 if (f in DOWN or f.startswith('COMPOSITE')) else 1
    e = res['NIH_Walitt2024'].loc[f] if f in res['NIH_Walitt2024'].index else None
    rows.append(dict(feature=f, expected_dir=exp, NIH_est=None if e is None else e.est, NIH_p=None if e is None else e.p,
                     same_dir=None if e is None else np.sign(e.est) == exp,
                     Raij_est=res['Raijmakers2020'].est.get(f), Columbia_est=res['Guo2023_Columbia'].est.get(f)))
A = pd.DataFrame(rows)
nsame = A.same_dir.sum(); ntot = A.same_dir.notna().sum()
print(A.round(3).to_string())
print(f'pre-specified species same direction in NIH: {nsame}/{ntot}, binomial p={stats.binomtest(int(nsame), int(ntot), 0.5, alternative="greater").pvalue:.3g}')
A.to_csv(OUT + 'gut_nih_prespecified.csv', index=False)

# (b) pre-specified signature score, no fitting
for k, (rel, cs, cv) in cohorts.items():
    X = clr(rel, keep=DOWN + UP); Z = (X - X.mean()) / X.std()
    s = Z[[u for u in UP if u in Z]].mean(axis=1) - Z[[d for d in DOWN if d in Z]].mean(axis=1)
    auc = roc_auc_score(cs, s)
    # bootstrap CI by person (one sample per person here)
    rng = np.random.default_rng(0); bs = []
    for _ in range(2000):
        i = rng.integers(0, len(cs), len(cs))
        if len(set(cs[i])) == 2: bs.append(roc_auc_score(cs[i], s.values[i]))
    perm = np.mean([roc_auc_score(rng.permutation(cs), s) >= auc for _ in range(5000)])
    print(f'{k}: signature AUC {auc:.3f} (95% CI {np.percentile(bs, 2.5):.2f}-{np.percentile(bs, 97.5):.2f}), perm p={perm:.4f}, n={len(cs)}')
    with open(OUT + 'gut_signature_auc.txt', 'a') as fh:
        fh.write(f'{k}\tAUC={auc:.3f}\tCI={np.percentile(bs, 2.5):.2f}-{np.percentile(bs, 97.5):.2f}\tperm_p={perm:.4f}\tn_case={cs.sum()}\tn_ctrl={(1 - cs).sum()}\n')

# (c) all shared species
shared = set(res['NIH_Walitt2024'].index) & set(res['Raijmakers2020'].index) & set(res['Guo2023_Columbia'].index)
shared = sorted(shared)
M = []
for f in shared:
    e = np.array([res[k].loc[f, 'est'] for k in cohorts]); s = np.array([res[k].loc[f, 'se'] for k in cohorts])
    w = 1 / s ** 2; fe = (w * e).sum() / w.sum(); Q = (w * (e - fe) ** 2).sum()
    tau2 = max(0, (Q - 2) / (w.sum() - (w ** 2).sum() / w.sum())); wr = 1 / (s ** 2 + tau2)
    re_ = (wr * e).sum() / wr.sum(); rse = np.sqrt(1 / wr.sum())
    M.append(dict(feature=f, **{k + '_est': res[k].loc[f, 'est'] for k in cohorts}, **{k + '_p': res[k].loc[f, 'p'] for k in cohorts},
                  meta_RE=re_, meta_se=rse, meta_p=2 * stats.norm.sf(abs(re_ / rse)), n_same_dir=int(max((e > 0).sum(), (e < 0).sum()))))
M = pd.DataFrame(M); M['meta_q'] = multipletests(M.meta_p, method='fdr_bh')[1]
M['all3_same_dir_and_each_p<0.1'] = (M.n_same_dir == 3) & (M[[k + '_p' for k in cohorts]] < 0.1).all(axis=1)
M.sort_values('meta_p').to_csv(OUT + 'gut_3cohort_meta.csv', index=False)
print('shared species', len(M), 'meta q<0.05:', (M.meta_q < 0.05).sum(), '; all 3 same dir & each p<0.1:', M['all3_same_dir_and_each_p<0.1'].sum())
print(M.sort_values('meta_p').head(25)[['feature', 'NIH_Walitt2024_est', 'Raijmakers2020_est', 'Guo2023_Columbia_est', 'meta_RE', 'meta_q', 'n_same_dir']].round(3).to_string())
# cross-cohort agreement of effect sizes (all shared species)
for a, b in [('NIH_Walitt2024', 'Raijmakers2020'), ('NIH_Walitt2024', 'Guo2023_Columbia'), ('Raijmakers2020', 'Guo2023_Columbia')]:
    rho, p = stats.spearmanr(M[a + '_est'], M[b + '_est'])
    print(f'effect-size correlation {a} vs {b}: rho={rho:.2f} p={p:.2g} (n species={len(M)})')
