from common import *
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_predict, KFold
T=pd.read_csv('patient_trajectories.csv',index_col=0)
M=meta.copy(); M['tp']=M.timepoints.str[2:].astype(int)
tp1=M[(M.tp==1)&(M.y==1)]; tp1=tp1.reset_index().set_index('subject')
T=T.join(tp1[['index','age','gender','bmi']])
# residualize change on baseline severity (regression to the mean)
b=np.polyfit(T.sev_first,T.delta_sev,1); T['delta_adj']=T.delta_sev-np.polyval(b,T.sev_first)
print('RTM slope',b[0].round(3), 'corr(first,delta)',stats.spearmanr(T.sev_first,T.delta_sev))
print('age vs improvement',stats.spearmanr(T.age,T.delta_adj)); print('duration',stats.mannwhitneyu(T[T.duration=='Short'].delta_adj,T[T.duration=='Long'].delta_adj))
print('gender',T.groupby('gender').delta_adj.mean())
# baseline symptoms predicting change
sc=score.loc[score.index.intersection(T['index'])]
for c in sc.columns:
    r,p=stats.spearmanr(sc.loc[T['index'],c].values,T.delta_adj.values); 
    if p<0.05: print('baseline symptom',c,round(r,2),round(p,4))
out=[]
for L in ['immune','specie','kegg','metabolome','quest']:
    X=load(L)
    if L=='specie': X=np.log10(X+1e-5); X=X.loc[:,(X>-4.9).mean()>0.2]
    if L=='kegg': X=np.log10(X+1)
    ids=[i for i in T['index'] if i in X.index]; TT=T[T['index'].isin(ids)]
    Xb=X.loc[TT['index']]; Xb=Xb.loc[:,Xb.std()>0]
    for f in Xb.columns:
        r,p=stats.spearmanr(Xb[f].values,TT.delta_adj.values); out.append((L,f,len(TT),r,p))
    # predictive check
    yv=TT.delta_adj.values; pr=cross_val_predict(RandomForestRegressor(300,min_samples_leaf=5,random_state=0,n_jobs=-1),Xb.fillna(0).values,yv,cv=KFold(5,shuffle=True,random_state=0))
    print(L,'n',len(TT),'CV prediction of improvement r=',round(stats.pearsonr(pr,yv)[0],3))
R=pd.DataFrame(out,columns=['layer','feature','n','rho','p']).dropna()
R['q']=np.nan
for L,g in R.groupby('layer'): R.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
R.sort_values('p').to_csv('baseline_predictors_of_change.csv',index=False)
print(R.groupby('layer').apply(lambda g:(g.q<0.1).sum()))
print(R[R.layer!='kegg'].sort_values('p').head(20).to_string())
