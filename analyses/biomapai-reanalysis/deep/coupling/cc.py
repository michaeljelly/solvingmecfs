"""Shared loaders for the gut-blood coupling analysis."""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
B = '/mnt/project-files/biomap/BioMapAI/'
I = B + 'codes/AI/input/'
OUT = '/mnt/project-files/biomap/deep/coupling/'

meta = pd.read_csv(I + 'metadata.csv', index_col=0)
meta['subject'] = meta.sample_id_tp1
meta['y'] = (meta.study_ptorhc == 'MECFS').astype(int)
meta['male'] = (meta.gender == 'Male').astype(int)
meta['tp'] = meta.timepoints.str[2:].astype(int)
meta['bmi_i'] = meta.bmi.fillna(meta.bmi.median())
meta['bmi_na'] = meta.bmi.isna().astype(int)
COV = ['age', 'male', 'bmi_i', 'bmi_na']

score = pd.read_csv(I + 'score.csv', index_col=0)
_s = score.copy()
_s['RAND36.Role Limitations due to Physical Health'] /= 4
_s['RAND36.Role Limitations due to Emotional Problems'] /= 3
score['severity'] = _s.apply(lambda c: (c - c.mean()) / c.std()).mean(axis=1)

mmeta = pd.read_csv(I + 'feature_meta/metabolome_feature_meta.csv', index_col=0)


def metab(max_imp=0.2):
    """QC-normalised Metabolon data, log2, drop metabolites with >max_imp fraction at minimum (imputed)."""
    q = pd.read_csv(B + 'data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv', index_col=0)
    ann = q.iloc[:, :17].copy()
    v = q.iloc[:, 17:].astype(float)
    frac = v.eq(v.min(axis=1), axis=0).mean(axis=1)
    keep = frac <= max_imp
    X = np.log2(v[keep]).T
    X.columns = q.CHEMICAL_NAME[keep].values
    ann = ann[keep].set_index('CHEMICAL_NAME')
    return X, ann


def species_clr(path=I + 'specie.csv'):
    s = pd.read_csv(path, index_col=0).T
    s = s.div(s.sum(axis=1), axis=0)
    pc = s[s > 0].min().min() / 2
    l = np.log(s + pc)
    return l.sub(l.mean(axis=1), axis=0), s


def species_full():
    f = pd.read_csv(B + 'data/metagenomics/metaphlan4_species_full.csv', index_col=0)
    tax = f.iloc[:, :7]
    a = f.iloc[:, 7:].astype(float).T
    a = a.div(a.sum(axis=1), axis=0)
    return a, tax


def kegg_logcpm():
    k = pd.read_csv(I + 'kegg.csv', index_col=0).T
    cpm = k.div(k.sum(axis=1), axis=0) * 1e6
    return np.log2(cpm + 1), k


def immune(raw=False):
    if raw:
        return pd.read_csv(B + 'data/immuno/immune_percentage.csv', index_col=0).T
    return pd.read_csv(I + 'immune.csv', index_col=0).T


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    o = np.argsort(pv)
    r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    q = np.empty(n); q[o] = np.minimum(r, 1)
    out[ok] = q
    return out


def resid(Y, C):
    """OLS residualise columns of Y (DataFrame) on covariate DataFrame C (in-sample)."""
    Cm = np.column_stack([np.ones(len(C)), C.values.astype(float)])
    beta, *_ = np.linalg.lstsq(Cm, Y.values, rcond=None)
    return pd.DataFrame(Y.values - Cm @ beta, index=Y.index, columns=Y.columns)
