"""Armstrong 2015 (MetaboLights MTBLS161) serum + urine NMR, 34 CFS / 25 controls (Australia). Copy from GitHub.
Tests ketones, lactate, acetate and amino acids; BH-FDR over all measured compounds per fluid. No covariates available."""
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
f = 'repos/markothaurelio_cfs-metabolomics-biomarker-discovery/m_MTBLS161_mecfs_metabolite_profiling.tsv'
d = pd.read_csv(f, sep='\t').set_index('metabolite_identification')
out = []
for fluid in ['serum', 'urine']:
    cols = [c for c in d.columns if c.endswith(fluid)]
    X = d[cols].apply(pd.to_numeric, errors='coerce').T
    X = X.loc[:, X.notna().mean() >= 0.6]
    case = X.index.str.startswith('CFS').astype(int)
    for m in X.columns:
        v = np.log2(X[m].where(X[m] > 0)).dropna(); c = case[X.index.get_indexer(v.index)]
        a, b = v[c == 1], v[c == 0]
        if len(a) < 5 or len(b) < 5: continue
        t = stats.ttest_ind(a, b, equal_var=False); sd = v.std()
        out.append(dict(fluid=fluid, metabolite=m, n_case=len(a), n_ctrl=len(b), smd=(a.mean() - b.mean()) / sd, p=t.pvalue,
                        mwu_p=stats.mannwhitneyu(a, b).pvalue))
R = pd.DataFrame(out); R['q'] = R.groupby('fluid').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
R.sort_values(['fluid', 'p']).to_csv('results/mtbls161_nmr.csv', index=False)
pd.set_option('display.width', 200); print(R.sort_values(['fluid', 'p']).round(3).to_string())
