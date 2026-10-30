"""Part 2b + 3b: characterise the multi-omics k-means partition: symptoms/duration/sex/age/meds, control-likeness, temporal persistence, trajectory."""
import sys
from sub_common import *
from clust_tools import *
from omics_space import *
from sklearn.metrics import cohen_kappa_score
from sklearn.neighbors import NearestNeighbors
import statsmodels.formula.api as smf
k=int(sys.argv[1]) if len(sys.argv)>1 else 2
OS=OmicsSpace(); E=OS.embed(OS.pt)
km=KMeans(k,n_init=50,random_state=0).fit(E); lab=km.labels_; cent=km.cluster_centers_
P=meta.loc[OS.pt].join(score); P['omk']=lab
print(f'k={k} sizes',np.bincount(lab))
rows=[]
for v in SC+['severity','age','bmi']:
    g=[P.loc[P.omk==j,v].dropna() for j in range(k)]; p=stats.kruskal(*g).pvalue
    rows.append((v,'kruskal',p,' | '.join(f'{x.median():.2f}' for x in g)))
for v in ['illness_duration','gender','antivirals','probiotics','antifungals','antibiotics','IBS']:
    t=pd.crosstab(P[v],P.omk); p=stats.chi2_contingency(t)[1] if t.shape[0]>1 else np.nan
    rows.append((v,'chi2',p,' | '.join(f"{t[j].to_dict()}" for j in range(k))))
R=pd.DataFrame(rows,columns=['var','test','p','by_cluster']); R['q']=bh(R.p.fillna(1)); pd.set_option('display.width',250); pd.set_option('display.max_colwidth',120)
print(R.round(4).to_string()); R.to_csv(OUT+f'results/omics_k{k}_clinical.csv',index=False)
# which layers/PCs drive the split
dims=[f'{L}_PC{i+1}' for L in OS.layers for i in range(10)]
cd=pd.DataFrame(cent,columns=dims); sep=(cd.max()-cd.min()).sort_values(ascending=False)
print('dims separating clusters most:',sep.head(6).round(3).to_dict())
# top individual features per layer (descriptive; circular by construction)
top=[]
for L in OS.layers:
    X=OS.X[L].loc[OS.pt]
    for f in X.columns:
        g=[X.loc[lab==j,f] for j in range(k)]; top.append((L,f,stats.kruskal(*g).pvalue,*(x.mean() for x in g)))
Tp=pd.DataFrame(top,columns=['layer','feature','p']+[f'mean_c{j}' for j in range(k)]); Tp['q']=Tp.groupby('layer').p.transform(bh)
Tp.to_csv(OUT+f'results/omics_k{k}_features.csv',index=False)
print('features q<0.05 per layer (circular, descriptive):',Tp[Tp.q<0.05].groupby('layer').size().to_dict())
for L in OS.layers: print(L,'top:',Tp[Tp.layer==L].nsmallest(5,'p')[['feature']+[f'mean_c{j}' for j in range(k)]].round(2).values.tolist())
# controls
ctl=[i for i in OS.all_idx if meta.loc[i,'y']==0 and meta.loc[i,'tp']==1]; Ec=OS.embed(ctl)
lc=((Ec[:,None,:]-cent[None])**2).sum(-1).argmin(1)
print('controls n',len(ctl),'assigned to clusters',np.bincount(lc,minlength=k),'patients',np.bincount(lab),'chi2 p',round(stats.chi2_contingency([np.bincount(lc,minlength=k),np.bincount(lab)])[1],4))
ctr=Ec.mean(0); spread=np.sqrt(((E-E.mean(0))**2).sum(1).mean())
print('distance cluster centroid -> control centroid (rel. to patient RMS spread):',[round(np.linalg.norm(c-ctr)/spread,3) for c in cent],' patient centroid->control',round(np.linalg.norm(E.mean(0)-ctr)/spread,3))
Eall=np.vstack([E,Ec]); isc=np.r_[np.zeros(len(E)),np.ones(len(Ec))]
nn=NearestNeighbors(n_neighbors=16).fit(Eall); _,ind=nn.kneighbors(E); cfrac=isc[ind[:,1:]].mean(1)
print('mean fraction of controls among 15 nearest neighbours by cluster:',[round(cfrac[lab==j].mean(),3) for j in range(k)],'overall control share',round(isc.mean(),3),'kruskal p',round(stats.kruskal(*[cfrac[lab==j] for j in range(k)]).pvalue,4))
# temporal persistence (later patient visits)
later=[i for i in OS.all_idx if meta.loc[i,'y']==1 and meta.loc[i,'tp']>1]
El=OS.embed(later); ll=((El[:,None,:]-cent[None])**2).sum(-1).argmin(1)
base=pd.Series(lab,index=meta.loc[OS.pt,'subject'])
Lr=pd.DataFrame({'subject':meta.loc[later,'subject'].values,'tp':meta.loc[later,'tp'].values,'lab':ll})
for tp in (2,3):
    d=Lr[Lr.tp==tp].set_index('subject'); s=d.index.intersection(base.index)
    a,b=base[s],d.loc[s,'lab']; exp=sum(a.value_counts(normalize=True).reindex(range(k),fill_value=0)*b.value_counts(normalize=True).reindex(range(k),fill_value=0))
    print(f'tp1->tp{tp}: n={len(s)} agreement={(a==b).mean():.3f} chance={exp:.3f} kappa={cohen_kappa_score(a,b):.3f}')
# test-retest of the embedding dimensions themselves
El_df=pd.DataFrame(El,index=later); E_df=pd.DataFrame(E,index=OS.pt)
t2=Lr[Lr.tp==2].copy(); t2['idx']=[later[i] for i in np.where(Lr.tp==2)[0]]
pairs=[(meta.index[(meta.subject==s)&(meta.tp==1)][0],i) for s,i in zip(t2.subject,t2.idx) if s in base.index]
r=[np.corrcoef([E_df.loc[a,j] for a,_ in pairs],[El_df.loc[b,j] for _,b in pairs])[0,1] for j in range(40)]
print('test-retest r of embedding PCs (tp1 vs tp2) per layer (mean of 10 PCs):',{L:round(np.mean(r[i*10:(i+1)*10]),2) for i,L in enumerate(OS.layers)})
# trajectory by omics cluster
Tr=pd.read_csv(OUT+'results/trajectory_table.csv',index_col=0); Tr['omk']=base.reindex(Tr.index)
Tr=Tr.dropna(subset=['omk'])
m=smf.ols('d_sev ~ C(omk)+base_sev+ntp+age+male',Tr).fit()
print('severity change by omics cluster (mean, n):',Tr.groupby('omk').d_sev.agg(['mean','size']).round(3).to_dict(),' adj p',round(m.f_test(' , '.join(f'C(omk)[T.{j}.0]=0' for j in range(1,k))).pvalue,4))
out=P[['omk']].copy(); out.to_csv(OUT+f'results/omics_k{k}_assignments.csv')
