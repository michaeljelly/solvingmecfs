import pandas as pd, numpy as np, warnings
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
L=pd.read_pickle(OUT+'L.pkl'); gm=pd.read_pickle(OUT+'gm.pkl'); fm=pd.read_pickle(OUT+'fm.pkl')
sub=pd.Series(fm['SUB PATHWAY'].values,index=fm.BIOCHEMICAL.astype(str))
d=pd.DataFrame({'ID':gm.ID,'case':(gm.Phenotype=='CFS').astype(int),'male':(gm.Sex=='M').astype(int),'age':gm.AgeBin.astype(float),'bmi':gm.BMIBin.astype(float),
 'day2':gm.Timepoint.str.startswith('D2').astype(int),'post':gm.Timepoint.str.endswith('POST').astype(int)})
terms=['case','case:post','case:day2','case:day2:post']
def fit(f):
    x=d.copy(); x['y']=L[f].values
    if not np.isfinite(x.y).all() or x.y.std()==0: return [f]+[np.nan]*8+['skip']
    x['y']=(x.y-x.y.mean())/x.y.std()
    fml='y ~ case*day2*post + male + age + bmi'
    for meth in ['lbfgs','powell','nm']:
        try:
            r=smf.mixedlm(fml,x,groups=x.ID).fit(reml=True,method=meth)
            if np.isfinite(r.pvalues[terms]).all(): return [f]+[v for t in terms for v in (r.params[t],r.pvalues[t])]+[meth]
        except Exception: pass
    r=smf.ols(fml,x).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(x.ID)[0]})
    return [f]+[v for t in terms for v in (r.params[t],r.pvalues[t])]+['ols_cluster']
rows=Parallel(n_jobs=2)(delayed(fit)(f) for f in L.columns)
R=pd.DataFrame(rows,columns=['feature']+[f'{a}_{t}' for t in terms for a in ('b','p')]+['fit']).set_index('feature')
for t in terms:
    ok=R['p_'+t].notna(); R.loc[ok,'q_'+t]=multipletests(R.loc[ok,'p_'+t],method='fdr_bh')[1]
    print(f'{t:16s} q<0.05: {(R["q_"+t]<0.05).sum():3d}  q<0.15: {(R["q_"+t]<0.15).sum():3d}  p<0.05: {(R["p_"+t]<0.05).sum():3d} of {ok.sum()}')
print(R.fit.value_counts()); R['sub']=sub.reindex(R.index).values
R.to_csv(OUT+'mixed_model_results.csv')
for t in terms[1:]:
    print('\n',t); print(R.sort_values('p_'+t).head(10)[['b_'+t,'p_'+t,'q_'+t,'sub']].round(4).to_string())
