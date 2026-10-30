import pandas as pd, numpy as np, warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
C=pd.read_pickle(OUT+'contrasts.pkl'); meta=pd.read_pickle(OUT+'meta.pkl'); fm=pd.read_pickle(OUT+'fm.pkl')
sub=pd.Series(fm['SUB PATHWAY'].values,index=fm.BIOCHEMICAL.astype(str))
rng=np.random.default_rng(1)
def fastols(Y,Xd,j):
    # Y n x p (no NaN), Xd n x k ; returns t, p for column j
    XtX=np.linalg.pinv(Xd.T@Xd); B=XtX@Xd.T@Y; Rz=Y-Xd@B; df=Xd.shape[0]-Xd.shape[1]
    s2=(Rz**2).sum(0)/df; se=np.sqrt(s2*XtX[j,j]); t=B[j]/np.where(se>0,se,np.nan)
    return B[j],t,2*stats.t.sf(np.abs(t),df)
def design(m,cols): return np.column_stack([np.ones(len(m))]+[m[c].values.astype(float) for c in cols])
cov=['case','male','AgeBin','BMIBin']
print('Permutation test (labels shuffled across persons, 500 perms): number of metabolites with p<0.05 for case vs control')
perm_rows=[]
for k in ['R1 (D1 post-pre)','R2 (D2 post-pre)','R2-R1','REC (D2pre-D1pre)','D1post->D2pre','CUM (D2post-D1pre)','LEVEL D1-PRE','LEVEL D2-PRE']:
    Y=C[k].dropna(); Y=Y.loc[:,Y.std()>1e-6]; m=meta.loc[Y.index]
    Xd=design(m,cov); _,_,p=fastols(Y.values,Xd,1); obs=(p<0.05).sum()
    null=[]
    for i in range(500):
        mm=m.copy(); mm['case']=rng.permutation(m.case.values)
        _,_,pp=fastols(Y.values,design(mm,cov),1); null.append((pp<0.05).sum())
    null=np.array(null); pv=(1+(null>=obs).sum())/501
    perm_rows.append((k,len(Y),Y.shape[1],obs,np.median(null),np.percentile(null,95),pv))
    print(f'{k:22s} n={len(Y)} obs={obs} null median={np.median(null):.0f} 95%={np.percentile(null,95):.0f} perm p={pv:.3f}')
pd.DataFrame(perm_rows,columns=['contrast','n','features','obs_p05','null_median','null_95','perm_p']).to_csv(OUT+'permutation_global.csv',index=False)

# consistency of trait differences D1PRE vs D2PRE
a=pd.read_pickle(OUT+'res_contrasts.pkl')
t1=a['LEVEL D1-PRE'].est/a['LEVEL D1-PRE'].se; t2=a['LEVEL D2-PRE'].est/a['LEVEL D2-PRE'].se
print('\ncorr of case-control t at D1PRE vs D2PRE:',round(np.corrcoef(t1.fillna(0),t2.fillna(0))[0,1],3))
tpo=a['LEVEL D1-POST'].est/a['LEVEL D1-POST'].se; tp2=a['LEVEL D2-POST'].est/a['LEVEL D2-POST'].se
print('corr D1PRE vs D2POST:',round(np.corrcoef(t1.fillna(0),tp2.fillna(0))[0,1],3))
# average of pre levels -> more power
Yavg=(C['LEVEL D1-PRE']+C['LEVEL D2-PRE'])/2; Yavg=Yavg.loc[:,Yavg.std()>1e-6]; m=meta.loc[Yavg.index]
b,t,p=fastols(Yavg.values,design(m,cov),1); R=pd.DataFrame({'est':b,'t':t,'p':p},index=Yavg.columns); R['q']=multipletests(R.p,method='fdr_bh')[1]; R['sub']=sub.reindex(R.index).values
print('\nMean of two resting samples, case vs ctrl: q<0.05',(R.q<0.05).sum()); R.sort_values('p').to_csv(OUT+'resting_mean_case_vs_ctrl.csv')
print(R.sort_values('p').head(40).round(4).to_string())
# subpathway enrichment of rest difference (mean t by subpathway)
R['absT']=R.t
g=R.groupby('sub').t.agg(['mean','count']); g=g[g['count']>=4]
g['p']=[stats.ttest_1samp(R.t[R['sub']==s],0).pvalue for s in g.index]; g['q']=multipletests(g.p,method='fdr_bh')[1]
print('\nSubpathways shifted at rest (mean t across members, q<0.1):'); print(g[g.q<0.1].sort_values('mean').round(4).to_string())
g.to_csv(OUT+'resting_subpathway_shift.csv')
