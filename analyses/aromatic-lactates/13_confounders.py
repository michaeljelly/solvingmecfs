"""Confounder handling for the metabolome-wide results.
Drug exposure is called from detected drug metabolites (Metabolon cohorts); coffee from caffeine; diet from plant/fish markers;
fasting from ketone bodies. Every harmonised metabolite is re-tested under several models and re-meta-analysed.
Models (per-sample cohorts only; Hoel summary data cannot be adjusted and is left out of all models here for a fair comparison):
  M0 base covariates (age/sex/BMI/site where available)
  M1 + caffeine
  M2 + drug classes detected (Metabolon cohorts; others get M1)
  M3 + caffeine + drugs + diet index + ketone index
  M4 base covariates, people with no prescription-type drug detected (gabapentinoid, antihistamine, beta-blocker, psychoactive, acid suppressant, other Rx); Metabolon cohorts; others full sample
Outputs: drug_exposure_by_group.csv, adjusted_meta.csv
"""
import re, pickle
import numpy as np, pandas as pd
from mwlib import *

S, M = load()
DRUG = {'acetaminophen': r'acetamin|acetamidophen', 'nsaid': r'ibuprofen|naproxen|meloxicam|celecoxib|diclofenac',
        'gabapentinoid': r'gabapentin|pregabalin', 'antihistamine': r'cetirizine|fexofenadine|diphenhydramine|loratadine|hydroxyzine',
        'betablocker': r'metoprolol|atenolol|propranolol', 'psychoactive': r'trazadone|trazodone|quetiapine|bupropion|venlafaxine|zolpidem|mirtazapine|diazepam|sertraline|citalopram|fluoxetine|duloxetine|amitriptyline|nortriptyline|lamotrigine|topiramate|levetiracetam|carbamazepine',
        'acid_suppr': r'omeprazole|pantoprazole|lansoprazole|ranitidine|famotidine', 'other_rx': r'montelukast|hydrochlorothiazide|valsartan|losartan|sitagliptin|metformin|warfarin|pseudoephedrine|guaifenesin|lidocaine|glycinexylidide|prednisol|azithromycin|lamivudine|cyclobenzaprine|chlorzoxazone'}
DIET = ['stachydrine', 'ergothioneine', '3-carboxy-4-methyl-5-propyl-2-furanpropanoate (CMPF)', 'proline betaine', 'homostachydrine*', 'beta-cryptoxanthin', 'docosahexaenoate (DHA; 22:6n3)',
        'Ergothioneine', 'Stachydrine', 'Docosahexaenoic acid', 'FA (22:6) (docosahexaenoic acid)', 'docosahexaenoic acid']
KET = ['3-hydroxybutyrate (BHBA)', 'acetoacetate', '3-hydroxybutyric acid', '3-Hydroxybutyric acid', '3_hydroxybutyric_acid_PM_a1', 'Acetoacetic acid', '3-hydroxybutyrate']
CAF = ['caffeine', 'Caffeine', 'caffeine_PM_a1']


def zmean(L, names):
    c = [n for n in names if n in L.columns]
    if not c: return None, []
    Z = (L[c] - L[c].mean()) / L[c].std()
    return Z.mean(axis=1), c


expo = []
MODELS = {}
for st, d in S.items():
    L, C, D, a = d['L'], d['C'].copy(), d['D'], d['ann']
    caf, _ = zmean(L, CAF); diet, dc = zmean(L, DIET); ket, kc = zmean(L, KET)
    if caf is not None: C['caf'] = caf
    if diet is not None: C['diet'] = diet
    if ket is not None: C['ket'] = ket
    dcols = []
    if d['family'] == 'Metabolon':
        drugf = [f for f in D.columns if f in a.index and str(a.loc[f, 'sub']).startswith('Drug')]
        for k, rx in DRUG.items():
            fs = [f for f in drugf if re.search(rx, f, re.I)]
            if not fs: continue
            C['drug_' + k] = D[fs].astype(bool).any(axis=1).astype(int).reindex(C.index).fillna(0)
            if C['drug_' + k].mean() >= 0.03: dcols.append('drug_' + k)
            expo.append(dict(study=st, drug_class=k, compounds='; '.join(fs[:6]), prev_case=C.loc[C.case == 1, 'drug_' + k].mean(), prev_ctrl=C.loc[C.case == 0, 'drug_' + k].mean()))
        RX = ['drug_' + k for k in ['gabapentinoid', 'antihistamine', 'betablocker', 'psychoactive', 'acid_suppr', 'other_rx'] if 'drug_' + k in C]
        C['drug_any'] = C[RX].max(axis=1)   # prescription-type drugs (paracetamol/NSAIDs not counted)
        expo.append(dict(study=st, drug_class='ANY', prev_case=C.loc[C.case == 1, 'drug_any'].mean(), prev_ctrl=C.loc[C.case == 0, 'drug_any'].mean(),
                         compounds=f'n drug-free: case {int(((C.drug_any == 0) & (C.case == 1)).sum())}, ctrl {int(((C.drug_any == 0) & (C.case == 0)).sum())}'))
    base = covs(st, C)
    extra_caf = ['caf'] if 'caf' in C else []
    m = {'M0': (base, None), 'M1': (base + extra_caf, None), 'M2': (base + dcols, None),
         'M3': (base + extra_caf + dcols + [c for c in ['diet', 'ket'] if c in C], None)}
    if 'drug_any' in C: m['M4'] = (base, C.drug_any == 0)
    else: m['M4'] = (base, None)
    MODELS[st] = (C, m, dc, kc)
    S[st]['C'] = C
pickle.dump({st: (S[st]['C'], MODELS[st][1]) for st in S}, open(W + 'C_adj.pkl', 'wb'))
E = pd.DataFrame(expo); E.to_csv(OUT + 'drug_exposure_by_group.csv', index=False)
print(E.round(2).to_string())

rows = []
for st, (C, m, dc, kc) in MODELS.items():
    L = S[st]['L']
    for mn, (cv, keep) in m.items():
        idx = C.index if keep is None else C.index[keep.values]
        Cc = C.loc[idx]
        if Cc.case.sum() < 5 or (1 - Cc.case).sum() < 5: continue
        cv = [c for c in cv if Cc[c].std() > 0]
        r = ols_effect(L.loc[idx], Cc, 'case', cv)
        r['study'] = st; r['model'] = mn; r['feature'] = r.index; r['n_case'] = int(Cc.case.sum()); r['n_ctrl'] = int((1 - Cc.case).sum())
        rows.append(r)
A = pd.concat(rows, ignore_index=True).merge(M[['hkey', 'study', 'feature']], on=['study', 'feature'])
A.to_csv(W + 'adjusted_per_study.csv', index=False)
out = []
for (mn, hk), g in A.groupby(['model', 'hkey']):
    r = meta(g.est, g.se); r.update(model=mn, hkey=hk, N=int((g.n_case + g.n_ctrl).sum())); out.append(r)
R = pd.DataFrame(out)
R = R.pivot(index='hkey', columns='model', values=['est', 'p', 'k', 'N'])
R.columns = [f'{b}_{a}' for a, b in R.columns]
R0 = pd.read_csv(OUT + 'meta_all.csv', index_col=0)
R = R0[['label', 'sub', 'all_k', 'all_est', 'all_p', 'all_q', 'rep_A', 'rep_B', 'rep_C']].join(R, how='inner')
for mn in ['M0', 'M1', 'M2', 'M3', 'M4']:
    ok = R[f'M0_k'] >= 3
    R[f'{mn}_q'] = np.nan; R.loc[ok, f'{mn}_q'] = bh(R.loc[ok, f'{mn}_p'])
R.to_csv(OUT + 'adjusted_meta.csv')
cand = R[((R.all_q < 0.05) & (R.all_k >= 3)) | R.rep_A | R.rep_B | R.rep_C]
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 40)
print(cand.sort_values('all_p')[['label', 'all_est', 'M0_est', 'M0_q', 'M1_est', 'M1_q', 'M2_est', 'M2_q', 'M3_est', 'M3_q', 'M4_est', 'M4_q', 'M4_N']].round(3).to_string())
