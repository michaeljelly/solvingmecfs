import pandas as pd, numpy as np, warnings
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
fm=pd.read_pickle(OUT+'fm.pkl'); gm=pd.read_pickle(OUT+'gm.pkl'); X=pd.read_pickle(OUT+'X_scaled.pkl'); Xo=pd.read_pickle(OUT+'X_orig.pkl')
print('in scaled not orig', set(X.columns)-set(Xo.columns))
det=Xo.notna().mean().reindex(X.columns).fillna(0)
sup=pd.Series(fm['SUPER PATHWAY'].values,index=fm.BIOCHEMICAL.astype(str)); sub=pd.Series(fm['SUB PATHWAY'].values,index=fm.BIOCHEMICAL.astype(str))
L=np.log2(X)
keep=det[det>=0.8].index
print('features total',X.shape[1],'detected>=80%',len(keep))
L=L[keep]; L=L.loc[:,~np.isinf(L).any()]
gm['case']=(gm.Phenotype=='CFS').astype(int); gm['male']=(gm.Sex=='M').astype(int)
gm['AgeBin']=gm.AgeBin.astype(float); gm['BMIBin']=gm.BMIBin.astype(float); gm['Bell']=gm.Bell.astype(float)
L.to_pickle(OUT+'L.pkl')
# wide per person
tp=['D1-PRE','D1-POST','D2-PRE','D2-POST']
W={t:L[gm.Timepoint.values==t].set_index(gm.ID[gm.Timepoint==t].values) for t in tp}
meta=gm.drop_duplicates('ID').set_index('ID')[['Phenotype','case','male','Sex','AgeBin','BMIBin','Bell']]
C={'R1 (D1 post-pre)':W['D1-POST']-W['D1-PRE'],
   'R2 (D2 post-pre)':W['D2-POST']-W['D2-PRE'],
   'R2-R1':(W['D2-POST']-W['D2-PRE'])-(W['D1-POST']-W['D1-PRE']),
   'REC (D2pre-D1pre)':W['D2-PRE']-W['D1-PRE'],
   'D1post->D2pre':W['D2-PRE']-W['D1-POST'],
   'CUM (D2post-D1pre)':W['D2-POST']-W['D1-PRE']}
for t in tp: C['LEVEL '+t]=W[t]
pd.to_pickle(C,OUT+'contrasts.pkl'); meta.to_pickle(OUT+'meta.pkl')
def ols_group(Y,m,term='case',extra=None,subset=None):
    rows=[]
    for f in Y.columns:
        y=Y[f].dropna(); mm=m.loc[y.index]
        if subset is not None: keepi=subset(mm); y=y[keepi]; mm=mm[keepi]
        Xd=mm[['case','male','AgeBin','BMIBin']].copy()
        if extra: 
            for e in extra: Xd[e[0]]=e[1](mm)
        Xd=sm.add_constant(Xd.astype(float))
        if Xd.shape[1]>1 and Xd.drop(columns='const').std().min()==0: Xd=Xd.loc[:,(Xd.std()>0)|(Xd.columns=='const')]
        r=sm.OLS(y.values,Xd.values).fit(cov_type='HC3')
        i=list(Xd.columns).index(term)
        rows.append((f,r.params[i],r.bse[i],r.pvalues[i],len(y)))
    R=pd.DataFrame(rows,columns=['feature','est','se','p','n']).set_index('feature')
    R.loc[R.se<1e-8,'p']=np.nan
    R['q']=np.nan; ok=R.p.notna(); R.loc[ok,'q']=multipletests(R.p[ok],method='fdr_bh')[1]
    return R
res={}
for k,Y in C.items():
    R=ols_group(Y,meta); R['sub']=sub.reindex(R.index).values; R['super']=sup.reindex(R.index).values
    # within-group mean change (for delta contrasts)
    if not k.startswith('LEVEL'):
        for g,gv in [('ctrl',0),('me',1)]:
            Yg=Y[meta.loc[Y.index,'case']==gv]
            t,p=stats.ttest_1samp(Yg,0,nan_policy='omit'); R['mean_'+g]=Yg.mean().values; R['p_'+g]=p; R['q_'+g]=multipletests(np.nan_to_num(p,nan=1),method='fdr_bh')[1]
    # Wilcoxon (unadjusted, paper-style)
    pw=[1.0 if Y[f].dropna().nunique()<3 else stats.mannwhitneyu(Y.loc[meta.loc[Y.index,'case']==1,f].dropna(),Y.loc[meta.loc[Y.index,'case']==0,f].dropna()).pvalue for f in R.index]
    R['p_wilcox']=pw; R['q_wilcox']=multipletests(pw,method='fdr_bh')[1]
    res[k]=R
    print(f'{k:22s} n={int(R.n.median())} group diff q<0.05: {(R.q<0.05).sum():4d}  q<0.15: {(R.q<0.15).sum():4d}  p<0.05: {(R.p<0.05).sum():4d} (exp {0.05*len(R):.0f})', end='')
    if 'q_ctrl' in R: print(f'  | within-ctrl q<.05 {(R.q_ctrl<0.05).sum()}  within-ME q<.05 {(R.q_me<0.05).sum()}')
    else: print()
pd.to_pickle(res,OUT+'res_contrasts.pkl')
with pd.ExcelWriter(OUT+'contrast_results.xlsx') as xw:
    for k,R in res.items(): R.sort_values('p').to_excel(xw,sheet_name=k.split(' ')[0][:30].replace('>','_'))
for k,R in res.items():
    print('\n==',k); print(R.sort_values('p').head(15)[['est','se','p','q','q_wilcox','sub']].round(4).to_string())
