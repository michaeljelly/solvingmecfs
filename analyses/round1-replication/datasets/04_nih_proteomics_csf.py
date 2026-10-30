"""NIH Walitt 2024 source data: (1) SomaScan plasma and CSF proteomics (published split by sex, recombined here,
OLS with sex (+age for plasma), BH-FDR); (2) CSF catecholamines and their relation to effort/grip measures
(re-test of the paper's main CSF claim with exact permutation-free Mann-Whitney and effect sizes)."""
import numpy as np, pandas as pd, rdata, statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

R = '/mnt/project-files/mecfs-research/datasets/repos/docwalitt_National-Institutes-of-Health-Myalgic-Encephalomyelitis-Chronic-Fatigue-Syndrome-Code-Repository/'
OUT = '/mnt/project-files/mecfs-research/datasets/results/'
d = rdata.read_rds(R + 'Nature Comms MECFS Source Data.RData')


def fix(df):
    df = df.copy(); df.columns = [str(c) for c in df.columns]
    for c in df.columns:
        if c not in ('SID', 'GRP', 'Gender', 'Group', 'Study.code'): df[c] = pd.to_numeric(df[c], errors='coerce').astype(float)
    return df


out = []
for fluid, a, b, va in [('plasma', 'Figure S15B_plasma_M', 'Figure S16B_plasma_F', 'Figure S15B_VIP_plasma_M'), ('CSF', 'Figure S15D_CSF_M', 'Figure S16D_CSF_F', 'Figure S15D_VIP_CSF_M')]:
    M = fix(d[a]); F = fix(d[b]); ann = d[va].copy(); ann.columns = [str(c) for c in ann.columns]; ann['SeqId'] = ann['SeqId'].astype(str)
    common = [c for c in M.columns if c in F.columns and c not in ('SID', 'GRP', 'Gender', 'Age')]
    X = pd.concat([M[['SID', 'GRP', 'Gender'] + common], F[['SID', 'GRP', 'Gender'] + common]], ignore_index=True)
    print(fluid, 'n', len(X), X.GRP.value_counts().to_dict(), X.Gender.value_counts().to_dict(), 'proteins', len(common))
    grp = X.GRP.astype(str); case = (~grp.str.upper().str.contains('HV|CONTROL|HEALTHY')).astype(int).values
    Y = np.log2(X[common].astype(float).clip(lower=1e-3))
    C = np.column_stack([np.ones(len(X)), case, (X.Gender.astype(str).str.upper().str.startswith('M')).astype(float)])
    rows = []
    for p in common:
        y = Y[p].to_numpy(dtype=float)
        if np.std(y) == 0 or not np.isfinite(y).all(): continue
        y = (y - y.mean()) / y.std(); r = sm.OLS(y, C).fit(); rows.append((p, r.params[1], r.pvalues[1]))
    T = pd.DataFrame(rows, columns=['SeqId', 'est_SD', 'p']); T['q'] = multipletests(T.p, method='fdr_bh')[1]

    T = T.merge(ann[['SeqId', 'EntrezGeneSymbol', 'TargetFullName']].drop_duplicates('SeqId'), on='SeqId', how='left')
    T['fluid'] = fluid; out.append(T)
    print(fluid, 'q<0.05:', (T.q < 0.05).sum(), 'q<0.2:', (T.q < 0.2).sum(), 'min q', T.q.min().round(3))
    print(T.sort_values('p').head(12)[['EntrezGeneSymbol', 'est_SD', 'p', 'q']].round(4).to_string())
P = pd.concat(out); P.to_csv(OUT + 'nih_somascan_plasma_csf.csv', index=False)

# ---------------- CSF catecholamines
cat = []
for k in ['Figure 6A', 'Figure 6B', 'Figure 6C', 'Figure 6D', 'Figure 6F']:
    df = fix(d[k]); col = [c for c in df.columns if c.startswith('CSF')][0]
    g = df.Group.astype(str)
    a = df[col][g == 'PI-ME/CFS'].dropna(); b = df[col][g == 'HV'].dropna()
    u = stats.mannwhitneyu(a, b); hl = np.median(np.subtract.outer(a.values, b.values))
    cat.append(dict(measure=col, n_case=len(a), n_hv=len(b), median_case=a.median(), median_hv=b.median(),
                    pct_diff=100 * (a.median() / b.median() - 1), MWU_p=u.pvalue, AUC_lower_in_case=1 - u.statistic / (len(a) * len(b))))
C = pd.DataFrame(cat); C['q'] = multipletests(C.MWU_p, method='fdr_bh')[1]
print(C.round(3).to_string()); C.to_csv(OUT + 'nih_csf_catecholamines.csv', index=False)
# within-group correlations reported in the paper (CSF NE vs grip time-to-failure etc.)
for k in ['Figure 6D', 'Figure 6E', 'Figure 6F', 'Figure 6G']:
    df = fix(d[k]); cols = [c for c in df.columns if c != 'Group']
    for gname in ['PI-ME/CFS', 'HV']:
        s = df[df.Group.astype(str) == gname][cols].dropna()
        if len(s) > 4:
            rho, p = stats.spearmanr(s[cols[0]], s[cols[1]]); print(f'{k} {gname}: {cols[0]} vs {cols[1]} rho={rho:.2f} p={p:.3f} n={len(s)}')
