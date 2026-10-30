"""Shared loaders/transforms for the subtype analysis."""
import pandas as pd, numpy as np
from scipy import stats
B='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT='/mnt/project-files/biomap/deep/subtypes/'
meta=pd.read_csv(B+'metadata.csv',index_col=0)
meta['subject']=meta.sample_id_tp1
meta['y']=(meta.study_ptorhc=='MECFS').astype(int)
meta['tp']=meta.timepoints.str[2:].astype(int)
meta['male']=(meta.gender=='Male').astype(int)
meta['bmi_miss']=meta.bmi.isna().astype(int)
meta['bmi_i']=meta.bmi.fillna(meta.bmi.median())
score=pd.read_csv(B+'score.csv',index_col=0)
score['RAND36.Role Limitations due to Physical Health']/=4
score['RAND36.Role Limitations due to Emotional Problems']/=3
SC=list(score.columns)
SHORT={c:c.split('.')[-1].replace('Role Limitations due to ','RoleLim ').replace(' Health','').replace('Emotional Problems','Emot') for c in SC}
score['severity']=score[SC].mean(axis=1)
def load(layer): return pd.read_csv(B+layer+'.csv',index_col=0).T
def int_transform(X):
    """rank-based inverse normal per column (robust to outliers)"""
    R=X.rank(axis=0); n=X.notna().sum()
    return pd.DataFrame(stats.norm.ppf((R-0.5)/n),index=X.index,columns=X.columns)
def layer(name):
    X=load(name)
    if name=='specie':
        tp1=X.index.str.endswith('_tp1')
        keep=(X[tp1]>0).mean()>=0.10
        X=np.log10(X.loc[:,keep]+1e-5)
        return (X-X.mean())/X.std()
    # drop near-constant features
    mode_frac=X.apply(lambda s:s.value_counts(normalize=True).iloc[0])
    X=X.loc[:,mode_frac<0.5]
    return int_transform(X)
def residualise(X,M,cols=('age','male','bmi_i'),fit_idx=None):
    """regress out covariates (fit on fit_idx rows, apply to all rows)"""
    D=np.column_stack([np.ones(len(M))]+[M[c].values.astype(float) for c in cols])
    fi=np.arange(len(M)) if fit_idx is None else fit_idx
    beta,*_=np.linalg.lstsq(D[fi],X.values[fi],rcond=None)
    return pd.DataFrame(X.values-D@beta,index=X.index,columns=X.columns)
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); n=len(p)
    q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]
    out=np.empty(n); out[o]=np.minimum(q,1); return out
