"""Part 3a: do patients keep their symptom partition over time? Later visits assigned to tp1 centroids (tp1 scaling)."""
from sub_common import *
from clust_tools import *
from sklearn.metrics import cohen_kappa_score
import statsmodels.formula.api as smf
A=pd.read_csv(OUT+'results/symptom_assignments_tp1.csv',index_col=0)
P1=meta.loc[A.index].join(score); S1=P1[SC]; mu,sd=S1.mean(),S1.std(); Z1=((S1-mu)/sd).values
allP=meta[meta.y==1].join(score,how='inner'); Zall=((allP[SC]-mu)/sd)
res=[]; rng=np.random.default_rng(0)
for k in (2,3):
    lab=A[f'sym_k{k}'].values; cent=np.vstack([Z1[lab==j].mean(0) for j in range(k)])
    d=((Zall.values[:,None,:]-cent[None])**2).sum(-1); allP[f'k{k}']=d.argmin(1)
    # sanity: tp1 reassignment agreement
    t1=allP[allP.tp==1]; print(k,'tp1 nearest-centroid agreement with kmeans',round((t1[f'k{k}'].values==A.loc[t1.index,f'sym_k{k}'].values).mean(),3))
    for tp in (2,3):
        a=allP[allP.tp==1].set_index('subject')[f'k{k}']; b=allP[allP.tp==tp].set_index('subject')[f'k{k}']; s=a.index.intersection(b.index)
        agree=(a[s]==b[s]).mean(); exp=sum(a[s].value_counts(normalize=True).reindex(range(k),fill_value=0)*b[s].value_counts(normalize=True).reindex(range(k),fill_value=0))
        kap=cohen_kappa_score(a[s],b[s]); ct=pd.crosstab(a[s],b[s])
        res.append(dict(k=k,tp=tp,n=len(s),agreement=agree,chance=exp,kappa=kap,table=ct.values.tolist()))
# benchmark: kappa expected if people simply keep their severity (continuous) with the observed test-retest correlation
sev1=allP[allP.tp==1].set_index('subject').severity; sev2=allP[allP.tp==2].set_index('subject').severity; s=sev1.index.intersection(sev2.index)
r=np.corrcoef(sev1[s],sev2[s])[0,1]; print('severity test-retest r tp1-tp2',round(r,3),'n',len(s))
for i in range(4):
    pc=[]
print(pd.DataFrame(res).round(3).to_string())
# symptom PCs test-retest (project later visits with tp1 PCA)
from sklearn.decomposition import PCA
pca=PCA(4).fit(Z1); PCall=pd.DataFrame(pca.transform(Zall.values),index=allP.index,columns=['PC1','PC2','PC3','PC4']); PCall['subject']=allP.subject; PCall['tp']=allP.tp
a=PCall[PCall.tp==1].set_index('subject'); b=PCall[PCall.tp==2].set_index('subject'); s=a.index.intersection(b.index)
print('symptom PC test-retest r (tp1 vs tp2):',{c:round(np.corrcoef(a.loc[s,c],b.loc[s,c])[0,1],2) for c in ['PC1','PC2','PC3','PC4']})
# trajectories: change in severity from tp1 to last visit by tp1 partition, adjusted for baseline severity (regression to the mean)
last=allP.sort_values('tp').groupby('subject').tail(1); last=last[last.tp>1].set_index('subject')
base=allP[allP.tp==1].set_index('subject'); s=last.index.intersection(base.index)
Tr=pd.DataFrame({'d_sev':last.loc[s,'severity']-base.loc[s,'severity'],'base_sev':base.loc[s,'severity'],'ntp':last.loc[s,'tp'],
   'k3':A.set_index(meta.loc[A.index,'subject']).loc[s,'sym_k3'].values,'age':base.loc[s,'age'],'male':base.loc[s,'male'],'dur':base.loc[s,'illness_duration']})
print('change by k3 (neg=improved):',Tr.groupby('k3').d_sev.agg(['mean','size']).round(3).to_dict())
m0=smf.ols('d_sev ~ C(k3)+ntp+age+male',Tr).fit(); m1=smf.ols('d_sev ~ C(k3)+base_sev+ntp+age+male',Tr).fit()
print('unadjusted k3 F-p',round(m0.f_test('C(k3)[T.1]=0, C(k3)[T.2]=0').pvalue,4),' adj for baseline severity F-p',round(m1.f_test('C(k3)[T.1]=0, C(k3)[T.2]=0').pvalue,4), 'base_sev beta',round(m1.params.base_sev,3))
pd.DataFrame(res).to_csv(OUT+'results/symptom_partition_temporal.csv',index=False); Tr.to_csv(OUT+'results/trajectory_table.csv')
print(m1.params.round(3).to_dict())
# sensitivity: is it regression to the mean on the items that define cluster 1 (emotional items)? recompute change without them
emo=['RAND36.Role Limitations due to Emotional Problems','RAND36.Emotional wellbeing']; oth=[c for c in SC if c not in emo]
Tr['d_sev_noemo']=last.loc[s,oth].mean(1)-base.loc[s,oth].mean(1); Tr['base_noemo']=base.loc[s,oth].mean(1)
Tr['d_emo']=last.loc[s,emo].mean(1)-base.loc[s,emo].mean(1)
m2=smf.ols('d_sev_noemo ~ C(k3)+base_noemo+ntp+age+male',Tr).fit()
print('non-emotional items: k3 F-p',round(m2.f_test('C(k3)[T.1]=0, C(k3)[T.2]=0').pvalue,4),m2.params.filter(like='k3').round(3).to_dict())
print('emotional items change by k3',Tr.groupby('k3').d_emo.mean().round(3).to_dict())
Tr.to_csv(OUT+'results/trajectory_table.csv')
