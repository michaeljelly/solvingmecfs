"""BioMapAI paired stool metagenomes: which gut species track the replicated blood metabolites?
Person level: CLR species and log2 metabolites averaged over the visits that have both.
Partial Spearman (ranks residualised on case, age, sex, BMI); BH within each target across species.
Genus/taxon sums tested as pre-specified candidates (Bifidobacterium, Lactobacillaceae, C. sporogenes, Peptostreptococcus, Akkermansia...).
Cross-validated (person-level 5x5 fold) elastic-net R2 from all species to each target.
Does the case-control difference in the aromatic-lactate module shrink after adjusting for the implicated taxa?
Also: correlation among gut-derived co-metabolites (BioMapAI, Germain 2022, Che) - do they move together?
Outputs: species_metabolite_partial.csv, taxa_candidates.csv, species_cv_r2.csv, cometabolite_correlations.csv
"""
import sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')
from mwlib import *
sys.path.insert(0, '/mnt/project-files/biomap/deep/stats_replication/')
from data import clr

S, M = load()
B = S['BioMapAI']; Lv, Mv = B['Lvisit'], B['Mvisit']
sp = pd.read_csv('/mnt/project-files/biomap/BioMapAI/data/metagenomics/metaphlan4_species_full.csv', index_col=0)
tax = sp[['phylum', 'family', 'genus']]
rel = sp.select_dtypes('number').T
common = rel.index.intersection(Lv.index)
print('visits with both', len(common))
relc = rel.loc[common]; relc = relc.div(relc.sum(axis=1), axis=0)
Xs = clr(relc, prev=0.1)
person = Mv.loc[common, 'person']
Xp = Xs.groupby(person.values).mean()
Lp = Lv.loc[common].groupby(person.values).mean()
C = B['C'].loc[Xp.index]
print('persons', len(Xp), 'case', int(C.case.sum()), 'species', Xs.shape[1])

lab2 = dict(zip(M[M.study == 'BioMapAI'].label, M[M.study == 'BioMapAI'].feature))
def col(l): return lab2.get(l, l)
T = {}
for l in ['indolelactate', 'indolepropionate', 'phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'imidazole lactate', 'alpha-hydroxyisocaproate', 'alpha-hydroxyisovalerate',
          'hippurate', 'p-cresol sulfate', 'phenylacetylglutamine', '3-indoxyl sulfate', 'trimethylamine N-oxide', 'imidazole propionate', 'trans-4-hydroxyproline',
          '5-(galactosylhydroxy)-lysine', 'indoleacetate', 'tryptophan', 'kynurenine']:
    if col(l) in Lp: T[l] = Lp[col(l)]
z = lambda s: (s - s.mean()) / s.std()
arom = ['phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'indolelactate', 'imidazole lactate']
T['MODULE aromatic lactates'] = pd.concat([z(T[a]) for a in arom if a in T], axis=1).mean(axis=1)
coll = [c for c in ['trans-4-hydroxyproline', '5-(galactosylhydroxy)-lysine', '4-hydroxyglutamate'] if col(c) in Lp]
T['MODULE collagen'] = pd.concat([z(Lp[col(c)]) for c in coll], axis=1).mean(axis=1)
if col('methionine sulfoxide') in Lp and col('methionine') in Lp: T['methionine sulfoxide/methionine'] = Lp[col('methionine sulfoxide')] - Lp[col('methionine')]
T = pd.DataFrame(T)

Z = pd.DataFrame({'case': C.case, 'age': C.age, 'male': C.male, 'bmi': C.bmi}).astype(float)
Zc = np.column_stack([np.ones(len(Z)), Z.values])
def resid_rank(D):
    R = D.rank().values
    b = np.linalg.lstsq(Zc, R, rcond=None)[0]
    return pd.DataFrame(R - Zc @ b, index=D.index, columns=D.columns)
RX = resid_rank(Xp); RT = resid_rank(T)
rows = []
for t in T.columns:
    y = RT[t].values
    r = (RX.values * (y[:, None] - y.mean())).sum(0) / np.sqrt(((RX.values - RX.values.mean(0)) ** 2).sum(0) * ((y - y.mean()) ** 2).sum())
    n = len(y); tt = r * np.sqrt((n - 6) / (1 - r ** 2)); p = 2 * stats.t.sf(np.abs(tt), n - 6)
    d = pd.DataFrame({'target': t, 'species': Xp.columns, 'rho_partial': r, 'p': p}); d['q'] = bh(p)
    rows.append(d)
P = pd.concat(rows)
P['genus'] = P.species.map(tax.genus); P['family'] = P.species.map(tax.family)
P.to_csv(OUT + 'species_metabolite_partial.csv', index=False)
pd.set_option('display.width', 220)
for t in T.columns:
    d = P[P.target == t].sort_values('p')
    print(f'\n## {t}: q<0.05 n={int((d.q < 0.05).sum())}; top:', '; '.join(f"{r.species} {r.rho_partial:+.2f} (q={r.q:.3f})" for r in d.head(8).itertuples()))

# candidate taxa (pre-specified), relative abundance sums -> log10(x + 1e-6), prevalence reported
rel_p = relc.groupby(person.values).mean().loc[Xp.index]
def taxsum(mask): return np.log10(rel_p.loc[:, mask].sum(axis=1) + 1e-6)
sname = pd.Series(rel_p.columns, index=rel_p.columns)
gen = tax.genus.reindex(rel_p.columns).astype(str); fam = tax.family.reindex(rel_p.columns).astype(str)
CAND = {'Bifidobacterium (genus)': gen.str.contains('Bifidobacterium').values,
        'Lactobacillaceae (family)': fam.str.contains('Lactobacillaceae').values,
        'Clostridium sporogenes': sname.str.contains('sporogenes').values,
        'Peptostreptococcus (genus)': gen.str.contains('Peptostreptococcus').values,
        'Akkermansia muciniphila': sname.str.contains('Akkermansia').values,
        'Faecalibacterium prausnitzii': sname.str.contains('prausnitzii').values,
        'Roseburia (genus)': gen.str.contains('Roseburia').values,
        'Streptococcus (genus)': gen.str.contains('Streptococcus').values,
        'Enterococcus (genus)': gen.str.contains('Enterococcus').values,
        'Bacteroides+Phocaeicola': gen.str.contains('Bacteroides|Phocaeicola').values,
        'Escherichia (genus)': gen.str.contains('Escherichia').values,
        'Ruminococcus gnavus': sname.str.contains('gnavus').values,
        'Anaerostipes (genus)': gen.str.contains('Anaerostipes').values}
TX = pd.DataFrame({k: taxsum(m) for k, m in CAND.items() if m.sum() > 0})
prev = pd.Series({k: float((rel_p.loc[:, m].sum(axis=1) > 0).mean()) for k, m in CAND.items() if m.sum() > 0})
RTX = resid_rank(TX)
rows = []
for t in T.columns:
    for k in TX.columns:
        r, p = stats.pearsonr(RTX[k], RT[t]); rows.append(dict(target=t, taxon=k, prevalence=prev[k], rho_partial=r, p=p))
TC = pd.DataFrame(rows); TC['q'] = bh(TC.p.values)
# case-control difference of each candidate taxon
cc = {}
for k in TX.columns:
    r = ols_effect(TX[[k]], Z.assign(case=Z.case), 'case', ['age', 'male', 'bmi']).iloc[0]; cc[k] = (r.est, r.p)
TC['taxon_case_SD'] = TC.taxon.map(lambda k: cc[k][0]); TC['taxon_case_p'] = TC.taxon.map(lambda k: cc[k][1])
TC.to_csv(OUT + 'taxa_candidates.csv', index=False)
print('\nCandidate taxa (partial Spearman, BH over all target x taxon):')
print(TC[TC.target.isin(['indolelactate', 'indolepropionate', 'phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'imidazole lactate', 'MODULE aromatic lactates', 'MODULE collagen', 'hippurate'])]
      .pivot(index='taxon', columns='target', values='rho_partial').round(2).to_string())
print(TC.drop_duplicates('taxon')[['taxon', 'prevalence', 'taxon_case_SD', 'taxon_case_p']].round(3).to_string())
print(TC[TC.q < 0.05].sort_values('p')[['target', 'taxon', 'rho_partial', 'q']].round(3).to_string())

# mediation-style: does adjusting for taxa shrink case effect on aromatic lactate module / ILA / IPA ?
med = []
for t in ['MODULE aromatic lactates', 'indolelactate', 'indolepropionate', 'phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'imidazole lactate']:
    base = ols_effect(T[[t]], Z, 'case', ['age', 'male', 'bmi']).iloc[0]
    top = P[(P.target == t)].sort_values('p').head(10).species.tolist()
    Zt = Z.join(Xp[top]); Zt.columns = [str(c) for c in Zt.columns]
    adj = ols_effect(T[[t]], Zt, 'case', ['age', 'male', 'bmi'] + [str(c) for c in top]).iloc[0]
    Zb = Z.join(TX[['Bifidobacterium (genus)', 'Lactobacillaceae (family)']].rename(columns=lambda c: c.split(' ')[0]))
    adjb = ols_effect(T[[t]], Zb, 'case', ['age', 'male', 'bmi', 'Bifidobacterium', 'Lactobacillaceae']).iloc[0]
    med.append(dict(target=t, case_SD=base.est, p=base.p, case_SD_adj_Bifido_Lacto=adjb.est, p_adj_BL=adjb.p, case_SD_adj_top10species=adj.est, p_adj10=adj.p))
print(pd.DataFrame(med).round(3).to_string())
pd.DataFrame(med).to_csv(OUT + 'taxa_mediation.csv', index=False)

# cross-validated elastic net species -> targets (person-level folds)
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import RepeatedKFold
cvr = []
Xa = Xp.values
for t in T.columns:
    y = T[t].values; pr = np.zeros(len(y)); r2s = []
    for rep, (tr, te) in enumerate(RepeatedKFold(n_splits=5, n_repeats=3, random_state=1).split(Xa)):
        mu, sd = Xa[tr].mean(0), Xa[tr].std(0) + 1e-9
        m = ElasticNetCV(l1_ratio=0.5, cv=4, alphas=30, max_iter=5000, random_state=0).fit((Xa[tr] - mu) / sd, y[tr])
        p_ = m.predict((Xa[te] - mu) / sd)
        r2s.append(1 - ((y[te] - p_) ** 2).sum() / ((y[te] - y[tr].mean()) ** 2).sum())
    cvr.append(dict(target=t, cv_R2_mean=np.mean(r2s), cv_R2_sd=np.std(r2s)))
CV = pd.DataFrame(cvr); CV.to_csv(OUT + 'species_cv_r2.csv', index=False); print(CV.round(3).to_string())

# co-metabolite correlation structure (partial on case) across cohorts
CM = ['indolelactate', 'indolepropionate', 'phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'imidazole lactate', 'hippurate', 'p-cresol sulfate',
      'phenylacetylglutamine', '3-indoxyl sulfate', 'trimethylamine N-oxide', 'imidazole propionate', 'indoleacetate', 'lactate', 'phenylpyruvate', 'alpha-hydroxyisocaproate']
out = []
for st in ['BioMapAI', 'Germain2022', 'Germain2018', 'Germain2020', 'Che2022', 'NagySzakal2018']:
    L = S[st]['L']; Cs = S[st]['C']; fm = dict(zip(M[M.study == st].label, M[M.study == st].feature))
    cols = {l: fm[l] for l in CM if l in fm and fm[l] in L}
    D = L[list(cols.values())].rank(); D.columns = list(cols.keys())
    X_ = np.column_stack([np.ones(len(D)), Cs.case.values])
    Rr = D.values - X_ @ np.linalg.lstsq(X_, D.values, rcond=None)[0]
    cr = pd.DataFrame(np.corrcoef(Rr.T), index=D.columns, columns=D.columns)
    for a in ['indolelactate', 'indolepropionate']:
        if a in cr:
            for b in cr.columns:
                if b != a: out.append(dict(study=st, a=a, b=b, r=cr.loc[a, b]))
CO = pd.DataFrame(out); CO.to_csv(OUT + 'cometabolite_correlations.csv', index=False)
print(CO.pivot_table(index=['a', 'b'], columns='study', values='r').round(2).to_string())
