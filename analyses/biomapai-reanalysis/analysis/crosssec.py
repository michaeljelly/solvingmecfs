from common import *
from scipy import stats
from statsmodels.stats.multitest import multipletests
mm=pd.read_csv('/mnt/project-files/biomap/BioMapAI/codes/AI/input/feature_meta/metabolome_feature_meta.csv',index_col=0)
km=pd.read_csv('/mnt/project-files/biomap/BioMapAI/codes/AI/input/feature_meta/kegg_feature_meta.csv',index_col=0)
s=score.copy()
s['RAND36.Role Limitations due to Physical Health']/=4; s['RAND36.Role Limitations due to Emotional Problems']/=3
s['severity']=s.mean(axis=1)
M=meta.join(s,how='left')
M['tp']=M.timepoints.str[2:].astype(int)
def ann(L,f):
    if L=='metabolome' and f in mm.index: return f"{mm.loc[f,'SUPER_PATHWAY']} / {mm.loc[f,'SUB_PATHWAY']}"
    if L=='kegg' and f in km.index: return str(km.loc[f,'definition'])[:70]
    return ''
def prep(L):
    X=load(L)
    if L=='specie': X=np.log10(X+1e-5); X=X.loc[:,(X>-4.9).mean()>0.2]
    if L=='kegg': X=np.log10(X+1)
    return X.loc[:,X.std()>0]
res=[]
for L in ['immune','specie','kegg','metabolome','quest']:
    X=prep(L); idx=X.index.intersection(M[M.tp==1].index); X1=X.loc[idx]; m1=M.loc[idx]
    for a,b in [('Short','Control'),('Long','Control'),('Short','Long')]:
        A=X1[m1.illness_duration==a]; Bb=X1[m1.illness_duration==b]
        for f in X1.columns:
            u,pv=stats.mannwhitneyu(A[f],Bb[f])
            delta=2*u/(len(A)*len(Bb))-1  # Cliff's delta: + means higher in a
            res.append((L,f'{a}_vs_{b}',f,len(A),len(Bb),delta,pv,ann(L,f)))
    # severity among patients at tp1 (between-person)
    pm=m1[m1.y==1]
    for t in ['severity','RAND36.Physical Functioning','NLT.Orthostatic Intolerance','PSQI.sleeping_problem','RAND36.Energy/Fatigue']:
        ok=pm[t].notna()
        for f in X1.columns:
            r,pv=stats.spearmanr(X1.loc[pm.index[ok],f],pm.loc[ok,t])
            res.append((L,f'patients_tp1_corr_{t}',f,ok.sum(),0,r,pv,ann(L,f)))
R=pd.DataFrame(res,columns=['layer','comparison','feature','nA','nB','effect','p','annotation']).dropna(subset=['p'])
R['q']=np.nan
for k,g in R.groupby(['layer','comparison']): R.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
R.sort_values('p').to_csv('cross_sectional_results.csv',index=False)
print(R.groupby(['comparison','layer']).apply(lambda g:(g.q<0.05).sum()).unstack())
for c in ['Short_vs_Control','Long_vs_Control','Short_vs_Long']:
    print('\n#####',c)
    print(R[(R.comparison==c)].sort_values('p').groupby('layer').head(12)[['layer','feature','effect','p','q','annotation']].to_string())
for t in ['severity','NLT.Orthostatic Intolerance','PSQI.sleeping_problem']:
    print('\n##### patients corr',t)
    print(R[R.comparison==f'patients_tp1_corr_{t}'].sort_values('p').head(20)[['layer','feature','effect','p','q','annotation']].to_string())
