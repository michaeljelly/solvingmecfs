"""Product/precursor ratios and pathway modules across all per-sample blood cohorts.
Ratio = log2(sum of numerators) - log2(sum of denominators) (raw-scale sums). Module = mean z of available members (>= half present, min 2).
Each tested per cohort with base covariates (M0) and + caffeine (M1); random-effects meta; BH across all ratios+modules.
Also: automatic sulfate-vs-glucuronide pairs (same parent compound) as a phase-II balance test.
Output: ratios_modules_meta.csv, ratios_per_study.csv
"""
import re
import numpy as np, pandas as pd
from mwlib import *

S, M = load()
feat = {st: dict(zip(g.label, g.feature)) for st, g in M.groupby('study')}   # label -> study feature
RATIOS = {
    # tryptophan
    'indolelactate / tryptophan': (['indolelactate'], ['tryptophan']),
    'IPA / tryptophan': (['indolepropionate'], ['tryptophan']),
    'indoleacetate / tryptophan': (['indoleacetate'], ['tryptophan']),
    '3-indoxyl sulfate / tryptophan': (['3-indoxyl sulfate'], ['tryptophan']),
    'IPA / indolelactate (reductive step)': (['indolepropionate'], ['indolelactate']),
    'kynurenine / tryptophan (IDO/TDO)': (['kynurenine'], ['tryptophan']),
    'kynurenate / kynurenine (KAT)': (['kynurenate'], ['kynurenine']),
    'quinolinate / kynurenine': (['quinolinate'], ['kynurenine']),
    'xanthurenate / kynurenine': (['xanthurenate'], ['kynurenine']),
    'anthranilate / kynurenine': (['anthranilate'], ['kynurenine']),
    'serotonin / tryptophan': (['serotonin'], ['tryptophan']),
    'C-glycosyltryptophan / tryptophan': (['C-glycosyltryptophan'], ['tryptophan']),
    # aromatic / branched 2-hydroxy acids vs keto acids (NADH-dependent reduction) and vs amino acids
    'phenyllactate / phenylalanine': (['phenyllactate (PLA)'], ['phenylalanine']),
    '4-OH-phenyllactate / tyrosine': (['3-(4-hydroxyphenyl)lactate'], ['tyrosine']),
    'imidazole lactate / histidine': (['imidazole lactate'], ['histidine']),
    'phenyllactate / phenylpyruvate': (['phenyllactate (PLA)'], ['phenylpyruvate']),
    '4-OH-phenyllactate / 4-OH-phenylpyruvate': (['3-(4-hydroxyphenyl)lactate'], ['4-hydroxyphenylpyruvate']),
    'phenylpyruvate / phenylalanine (transamination)': (['phenylpyruvate'], ['phenylalanine']),
    'alpha-hydroxyisocaproate / 4-methyl-2-oxopentanoate': (['alpha-hydroxyisocaproate'], ['4-methyl-2-oxopentanoate']),
    '2-hydroxy-3-methylvalerate / 3-methyl-2-oxovalerate': (['2-hydroxy-3-methylvalerate'], ['3-methyl-2-oxovalerate']),
    'alpha-hydroxyisovalerate / 3-methyl-2-oxobutyrate': (['alpha-hydroxyisovalerate'], ['3-methyl-2-oxobutyrate']),
    'lactate / pyruvate (cytosolic NADH/NAD+)': (['lactate'], ['pyruvate']),
    '3-hydroxybutyrate / acetoacetate (mito NADH/NAD+)': (['3-hydroxybutyrate (BHBA)'], ['acetoacetate']),
    # NO / urea cycle
    'citrulline / arginine (NOS proxy)': (['citrulline'], ['arginine']),
    'ornithine / arginine (arginase)': (['ornithine'], ['arginine']),
    'arginine / (ornithine + citrulline) (GABR)': (['arginine'], ['ornithine', 'citrulline']),
    'ADMA+SDMA / arginine': (['dimethylarginine (SDMA + ADMA)'], ['arginine']),
    'homoarginine / arginine': (['homoarginine'], ['arginine']),
    # energy
    'creatine / creatinine': (['creatine'], ['creatinine']),
    'guanidinoacetate / arginine (AGAT)': (['guanidinoacetate'], ['arginine']),
    'creatine / guanidinoacetate (GAMT)': (['creatine'], ['guanidinoacetate']),
    'acetylcarnitine / carnitine': (['acetylcarnitine (C2)'], ['carnitine']),
    'propionylcarnitine / carnitine': (['propionylcarnitine (C3)'], ['carnitine']),
    '(C16+C18:1) / carnitine (CPT1)': (['palmitoylcarnitine (C16)', 'oleoylcarnitine (C18:1)'], ['carnitine']),
    'deoxycarnitine / carnitine (BBOX)': (['deoxycarnitine'], ['carnitine']),
    'glutamine / glutamate': (['glutamine'], ['glutamate']),
    # glycine conjugation / microbial
    'hippurate / benzoate (glycine conjugation)': (['hippurate'], ['benzoate']),
    'phenylacetylglutamine / phenylalanine': (['phenylacetylglutamine'], ['phenylalanine']),
    'p-cresol sulfate / tyrosine': (['p-cresol sulfate'], ['tyrosine']),
    'TMAO / choline': (['trimethylamine N-oxide'], ['choline']),
    'TMAO / carnitine': (['trimethylamine N-oxide'], ['carnitine']),
    'imidazole propionate / histidine': (['imidazole propionate'], ['histidine']),
    'p-cresol glucuronide / p-cresol sulfate': (['p-cresol glucuronide*'], ['p-cresol sulfate']),
    # polyamines, methylation, purines
    'N-acetylputrescine / ornithine': (['N-acetylputrescine'], ['ornithine']),
    'MTA / spermidine': (['5-methylthioadenosine (MTA)'], ['spermidine']),
    'SAH / methionine': (['S-adenosylhomocysteine (SAH)'], ['methionine']),
    'hypoxanthine / xanthine': (['hypoxanthine'], ['xanthine']),
    'urate / xanthine (XO)': (['urate'], ['xanthine']),
    'inosine / hypoxanthine': (['inosine'], ['hypoxanthine']),
    # taurine, NAD, redox
    'taurine / hypotaurine': (['taurine'], ['hypotaurine']),
    'taurine / cysteine': (['taurine'], ['cysteine']),
    '1-methylnicotinamide / nicotinamide (NNMT)': (['1-methylnicotinamide'], ['nicotinamide']),
    '2PY / nicotinamide (NAD turnover)': (['N1-methyl-2-pyridone-5-carboxamide'], ['nicotinamide']),
    'quinolinate / tryptophan (de novo NAD)': (['quinolinate'], ['tryptophan']),
    'methionine sulfoxide / methionine (oxidative stress)': (['methionine sulfoxide'], ['methionine']),
    '5-oxoproline / glutamate': (['5-oxoproline'], ['glutamate']),
    'Fischer ratio BCAA / (Phe+Tyr)': (['leucine', 'isoleucine', 'valine'], ['phenylalanine', 'tyrosine']),
}
MODULES = {
    'MODULE 2-hydroxy acids (aromatic + branched lactates)': ['phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'indolelactate', 'imidazole lactate', 'alpha-hydroxyisocaproate', 'alpha-hydroxyisovalerate', '2-hydroxy-3-methylvalerate'],
    'MODULE aromatic lactates (PLA, OH-PLA, ILA, imidazole lactate)': ['phenyllactate (PLA)', '3-(4-hydroxyphenyl)lactate', 'indolelactate', 'imidazole lactate'],
    'MODULE branched-chain 2-hydroxy acids': ['alpha-hydroxyisocaproate', 'alpha-hydroxyisovalerate', '2-hydroxy-3-methylvalerate'],
    'MODULE collagen turnover (Hyp, Gal-Hyl, 4-OH-Glu, HOG)': ['trans-4-hydroxyproline', '5-(galactosylhydroxy)-lysine', '4-hydroxyglutamate', '4-hydroxy-2-oxoglutaric acid', '5-hydroxylysine'],
    'MODULE microbial co-metabolites (all)': ['p-cresol sulfate', 'phenylacetylglutamine', 'trimethylamine N-oxide', 'hippurate', 'imidazole propionate', 'indolepropionate', 'indolelactate', '3-indoxyl sulfate', '4-ethylphenylsulfate', 'phenol sulfate'],
    'MODULE proteolytic fermentation (pCS, PAG, IS, phenol sulfate, 4-EPS, ImP)': ['p-cresol sulfate', 'phenylacetylglutamine', '3-indoxyl sulfate', 'phenol sulfate', '4-ethylphenylsulfate', 'imidazole propionate'],
    'MODULE IPA + indolelactate + hippurate': ['indolepropionate', 'indolelactate', 'hippurate'],
    'MODULE ketone bodies': ['3-hydroxybutyrate (BHBA)', 'acetoacetate'],
    'MODULE caffeine metabolites': ['caffeine', 'paraxanthine', 'theophylline', '1-methylxanthine', '1,7-dimethylurate', '1-methylurate'],
    'MODULE androgen/pregnenolone sulfates': ['dehydroepiandrosterone sulfate (DHEA-S)', 'androsterone sulfate', 'epiandrosterone sulfate', 'androstenediol (3beta,17beta) disulfate (1)', '5alpha-androstan-3beta,17beta-diol disulfate', 'pregnenolone sulfate'],
    'MODULE plasmalogens (P-16:0)': ['1-(1-enyl-palmitoyl)-2-linoleoyl-GPC (P-16:0/18:2)*', '1-(1-enyl-palmitoyl)-2-arachidonoyl-GPC (P-16:0/20:4)*', '1-(1-enyl-palmitoyl)-2-oleoyl-GPC (P-16:0/18:1)*', '1-(1-enyl-palmitoyl)-2-palmitoyl-GPC (P-16:0/16:0)*'],
    'MODULE long-chain acylcarnitines': ['palmitoylcarnitine (C16)', 'stearoylcarnitine (C18)', 'oleoylcarnitine (C18:1)', 'linoleoylcarnitine (C18:2)*', 'myristoylcarnitine (C14)', 'arachidonoylcarnitine (C20:4)'],
    'MODULE glutathione/oxidative (5-oxoproline, 2-HB, cys-GSH, Met-SO)': ['5-oxoproline', '2-hydroxybutyrate/2-hydroxyisobutyrate', 'cysteine-glutathione disulfide', 'methionine sulfoxide'],
}
# automatic sulfate vs glucuronide pairs (same parent name) in Metabolon labels
labels = set(M.label)
PAIRS = {}
for l in labels:
    m = re.fullmatch(r'(.+?)\s*glucuronide\*?', l)
    if m:
        base = m.group(1).strip()
        for cand in [base + ' sulfate', base + ' sulfate*', base + 'sulfate']:
            if cand in labels: PAIRS[f'PAIR {l} / {cand}'] = ([l], [cand])
RATIOS.update(PAIRS)


def lab2col(st, lab):
    f = feat[st].get(lab)
    return f if (f is not None and f in S[st]['L'].columns) else None


import pickle
CADJ = pickle.load(open(W + 'C_adj.pkl', 'rb'))
rows = []
for st, d in S.items():
    L = d['L']; C, MOD = CADJ[st]; C = C.copy()
    caf = lab2col(st, 'caffeine')
    if caf: C['caf'] = (L[caf] - L[caf].mean()) / L[caf].std()
    Y = {}
    for name, (num, den) in RATIOS.items():
        nc = [lab2col(st, x) for x in num]; dc = [lab2col(st, x) for x in den]
        if any(c is None for c in nc) or any(c is None for c in dc): continue
        Y[name] = np.log2(np.exp2(L[nc]).sum(axis=1)) - np.log2(np.exp2(L[dc]).sum(axis=1))
    for name, mem in MODULES.items():
        cs = [c for c in (lab2col(st, x) for x in mem) if c]
        if len(cs) >= max(2, int(np.ceil(len(mem) / 2))) or (len(cs) >= 2 and 'ketone' in name):
            Z = (L[cs] - L[cs].mean()) / L[cs].std(); Y[name] = Z.mean(axis=1)
    if not Y: continue
    Y = pd.DataFrame(Y)
    for mn in ['M0', 'M1', 'M2', 'M3', 'M4']:
        cv, keep = MOD[mn]
        idx = C.index if keep is None else C.index[keep.values]
        cv = [c for c in cv if C.loc[idx, c].std() > 0]
        r = ols_effect(Y.loc[idx], C.loc[idx], 'case', cv)
        r['study'] = st; r['model'] = mn; r['ratio'] = r.index; r['n_case'] = int(C.loc[idx].case.sum()); r['n_ctrl'] = int((1 - C.loc[idx].case).sum())
        rows.append(r)
A = pd.concat(rows, ignore_index=True)
A.to_csv(OUT + 'ratios_per_study.csv', index=False)
out = []
for (rt, mn), g in A.groupby(['ratio', 'model']):
    r = meta(g.est, g.se); r.update(ratio=rt, model=mn, N=int((g.n_case + g.n_ctrl).sum()), studies=','.join(g.study),
                              n_neg=int((g.est < 0).sum()), n_pos=int((g.est > 0).sum()))
    # leave-BioMapAI-out
    g2 = g[g.study != 'BioMapAI']; r2 = meta(g2.est, g2.se); r['exBM_est'] = r2['est']; r['exBM_p'] = r2['p']; r['exBM_k'] = r2['k']
    out.append(r)
R = pd.DataFrame(out)
for mn in ['M0', 'M1', 'M2', 'M3', 'M4']:
    m = (R.model == mn) & (R.k >= 2)
    R.loc[m, 'q'] = bh(R.loc[m, 'p'])
R = R.sort_values(['model', 'p'])
R.to_csv(OUT + 'ratios_modules_meta.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
print(R[(R.model == 'M0') & (R.k >= 2)][['ratio', 'k', 'N', 'est', 'se', 'p', 'q', 'I2', 'n_neg', 'n_pos', 'exBM_est', 'exBM_p']].round(3).head(60).to_string())
P = R[R.k >= 2].pivot(index='ratio', columns='model', values='est').join(R[R.k >= 2].pivot(index='ratio', columns='model', values='q'), rsuffix='_q')
top = R[(R.model == 'M0') & (R.q < 0.1)].ratio
print(P.loc[top].round(3).to_string())
