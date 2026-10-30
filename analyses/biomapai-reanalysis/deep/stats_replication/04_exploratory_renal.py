"""Exploratory follow-up: the metabolites most associated with between-person severity/sleep/cognition are
classic kidney-filtration markers (C-glycosyltryptophan, pseudouridine, N6-carbamoylthreonyladenosine, ...).
Build a 'filtration-marker' composite, test vs symptoms (between/within), vs clinical eGFR/creatinine,
with BMI adjustment, and in Germain 2022 vs Bell score and case status."""
import numpy as np, pandas as pd
from scipy import stats
from data import *
from lmmfast import fit_many

MARK = ['C-glycosyltryptophan', 'pseudouridine', 'N6-carbamoylthreonyladenosine', 'N2,N2-dimethylguanosine', 'N-acetylserine',
        'erythronate*', 'N1-methylinosine', 'N-acetylalanine', '7-methylguanine', 'myo-inositol', 'N-acetylthreonine']
m = load_meta()
X = load_layer('metabolome'); mk = [c for c in MARK if c in X]; print('markers used', mk)
z = (X[mk] - X[mk].mean()) / X[mk].std()
comp = z.mean(axis=1).rename('filtration_markers')
Q = load_layer('quest')[['eGFR non-African American', 'Creatinine']]
D = m.join(comp, how='inner').join(Q, how='left')
out = []
T = D[D.tp == 1]
r, p = stats.spearmanr(T.filtration_markers, T['eGFR non-African American'], nan_policy='omit'); out.append(('BioMapAI tp1 all', 'spearman vs eGFR (quest, residualised)', r, p, len(T)))
# case-control
C = np.column_stack([np.ones(len(D)), D.case, D.age_c, D.male, D.bmi_c, D.bmi_missing, D.tp2, D.tp3])
rr = fit_many(((D.filtration_markers - D.filtration_markers.mean()) / D.filtration_markers.std()).values[:, None], C, D.person.values, L=np.eye(8)[[1]])
out.append(('BioMapAI all visits', 'patient vs control (SD)', rr['est'][0, 0], rr['p'][0, 0], D.person.nunique()))
P = D[(D.case == 1) & D.severity.notna()]
f = (P.filtration_markers - D.filtration_markers.mean()) / D.filtration_markers.std()
bt = f.groupby(P.person).transform('mean'); wi = f - bt
for s, nm in SYMPTOMS.items():
    for adj in ['age+sex', 'age+sex+BMI+eGFRlab']:
        PP = P.copy(); keep = np.ones(len(PP), bool)
        cols = [np.ones(len(PP)), PP.age_c, PP.male, PP.tp2, PP.tp3]
        if adj != 'age+sex':
            keep = PP['eGFR non-African American'].notna().values
            cols += [PP.bmi_c, PP.bmi_missing, PP['eGFR non-African American'].fillna(0)]
        Cc = np.column_stack(cols)[keep]
        r2 = fit_many(PP[s].values[keep], Cc, PP.person.values[keep], F=np.stack([bt.values[keep], wi.values[keep]], 1)[None])
        out.append(('BioMapAI patients', f'{nm} ~ between ({adj})', r2['est'][0, 0] / P[s].std(), r2['p'][0, 0], PP.person[keep].nunique()))
        out.append(('BioMapAI patients', f'{nm} ~ within ({adj})', r2['est'][0, 1] / P[s].std(), r2['p'][0, 1], PP.person[keep].nunique()))
# Germain
IC = B + 'data/indenpendent_cohorts/'
d = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='ScaledImpDataZeroDrug&Tobacco', header=None)
fm = d.iloc[8:, :15].copy(); fm.columns = d.iloc[7, :15].values
gm = d.iloc[:7, 15:].T; gm.columns = gm.iloc[0]; gm = gm.iloc[1:]
gm.columns = ['ID', 'Phenotype', 'Timepoint', 'Sex', 'AgeBin', 'BMIBin', 'Bell']
G = d.iloc[8:, 16:].astype(float).T; G.columns = fm.BIOCHEMICAL.astype(str).values; G.index = gm.index
G = np.log2(G.loc[:, ~G.columns.duplicated()])
gk = [c for c in mk if c in G]; print('Germain markers', gk)
pre = gm.Timepoint.str.endswith('PRE').values; Gp = G[pre]; gp = gm[pre]
gc = ((Gp[gk] - Gp[gk].mean()) / Gp[gk].std()).mean(axis=1); gc = (gc - gc.mean()) / gc.std()
case = (gp.Phenotype == 'CFS').astype(int).values
C = np.column_stack([np.ones(len(gp)), case, gp.AgeBin.astype(float), (gp.Sex == 'M').astype(float), gp.BMIBin.astype(float), (gp.Timepoint == 'D2-PRE').astype(float)])
rr = fit_many(gc.values[:, None], C, gp.ID.values, L=np.eye(6)[[1]])
out.append(('Germain2022', 'patient vs control (SD)', rr['est'][0, 0], rr['p'][0, 0], gp.ID.nunique()))
cf = gp.Phenotype == 'CFS'
pmn = gc[cf].groupby(gp.ID[cf]).mean(); bell = gp.Bell[cf].astype(float).groupby(gp.ID[cf]).mean()
r, p = stats.spearmanr(pmn, bell.loc[pmn.index]); out.append(('Germain2022 cases', 'spearman vs Bell (higher Bell = better)', r, p, len(pmn)))
O = pd.DataFrame(out, columns=['cohort', 'test', 'estimate', 'p', 'n_persons'])
O.to_csv(OUT + 'exploratory_filtration_markers.csv', index=False)
print(O.round(4).to_string(index=False))
