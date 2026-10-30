"""Follow-up on the only configurations that beat the null a little: species-only k=3 and SNF k=3."""
from sub_common import *
from clust_tools import *
from omics_space import *
from sklearn.cluster import SpectralClustering
from sklearn.metrics import cohen_kappa_score
OS=OmicsSpace(); per=OS.embed(OS.pt,per_layer=True); w=OS.models['specie'][4]
Xs=per['specie']*w; km=KMeans(3,n_init=50,random_state=0).fit(Xs); lab=km.labels_
P=meta.loc[OS.pt].join(score); P['gk']=lab
sp=OS.X['specie'].loc[OS.pt]
T=pd.DataFrame({f:[stats.kruskal(*[sp.loc[lab==j,f] for j in range(3)]).pvalue] for f in sp.columns}).T[0].nsmallest(8)
print('gut k3 sizes',np.bincount(lab)); print(sp[T.index].groupby(lab).mean().round(2).T.to_string())
ax=pd.read_csv(OUT+'results/axes_tp1_values.csv',index_col=0); J=P.join(ax[['axis_dysbiosis','butyrate_producers','R_gnavus']])
print(J.groupby('gk')[['axis_dysbiosis','butyrate_producers']].mean().round(2).to_dict(),'Rg detected',J.groupby('gk').R_gnavus.apply(lambda s:(s>-4.99).mean()).round(2).to_dict())
res=[]
for v in SC+['severity','age']:
    res.append((v,stats.kruskal(*[P.loc[P.gk==j,v] for j in range(3)]).pvalue,[round(P.loc[P.gk==j,v].median(),2) for j in range(3)]))
for v in ['illness_duration','gender','antivirals','probiotics','IBS']:
    res.append((v,stats.chi2_contingency(pd.crosstab(P[v],P.gk))[1],pd.crosstab(P[v],P.gk).values.tolist()))
R=pd.DataFrame(res,columns=['var','p','by']); R['q']=bh(R.p); print(R.round(3).to_string())
ctl=[i for i in OS.all_idx if meta.loc[i,'y']==0 and meta.loc[i,'tp']==1]
# controls with species data (not requiring all layers)
Xsp=OS.X['specie']; ctl=[i for i in Xsp.index if meta.loc[i,'y']==0 and meta.loc[i,'tp']==1]
beta,mu,sd,pca,wt=OS.models['specie']; M=meta.loc[ctl]; D=np.column_stack([np.ones(len(M)),M[COV].values.astype(float)])
Ec=pca.transform(((Xsp.loc[ctl].values-D@beta)-mu)/sd); lc=((Ec[:,None,:]-km.cluster_centers_[None])**2).sum(-1).argmin(1)
print('controls per gut cluster',np.bincount(lc,minlength=3),'patients',np.bincount(lab),'chi2 p',round(stats.chi2_contingency([np.bincount(lc,minlength=3),np.bincount(lab)])[1],4))
later=[i for i in Xsp.index if meta.loc[i,'y']==1 and meta.loc[i,'tp']==2]; M=meta.loc[later]; D=np.column_stack([np.ones(len(M)),M[COV].values.astype(float)])
El=pca.transform(((Xsp.loc[later].values-D@beta)-mu)/sd); ll=pd.Series(((El[:,None,:]-km.cluster_centers_[None])**2).sum(-1).argmin(1),index=M.subject)
b=pd.Series(lab,index=meta.loc[OS.pt,'subject']); s=b.index.intersection(ll.index)
print('gut k3 tp1->tp2 n',len(s),'agreement',round((b[s]==ll[s]).mean(),3),'kappa',round(cohen_kappa_score(b[s],ll[s]),3))
P[['gk']].to_csv(OUT+'results/gut_k3_assignments.csv')
# SNF k=3
W=np.load(OUT+'results/snf_W.npy'); sl=SpectralClustering(3,affinity='precomputed',random_state=0).fit_predict(W); P['snf']=sl
print('SNF3 sizes',np.bincount(sl),'ARI vs gut k3',round(ARI(sl,lab),3))
res=[]
for v in SC+['severity','age']: res.append((v,stats.kruskal(*[P.loc[P.snf==j,v] for j in range(3)]).pvalue))
for v in ['illness_duration','gender']: res.append((v,stats.chi2_contingency(pd.crosstab(P[v],P.snf))[1]))
R=pd.DataFrame(res,columns=['var','p']); R['q']=bh(R.p); print('SNF3 min q',R.q.min().round(3),R.nsmallest(3,'p').round(3).values.tolist())
