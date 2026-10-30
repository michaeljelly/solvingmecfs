import pandas as pd, numpy as np, warnings
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score, silhouette_score
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
C=pd.read_pickle(OUT+'contrasts.pkl'); meta=pd.read_pickle(OUT+'meta.pkl'); S=pd.read_pickle(OUT+'sets.pkl')
sd=C['LEVEL D1-PRE'].std(); sd=sd.where(sd>0.05)
def setscores(Y):
    out={}
    for s,mem in S.items():
        mem=[m for m in mem if pd.notna(sd.get(m))]
        if mem: out[s]=(Y[mem]/sd[mem]).mean(axis=1)
    return pd.DataFrame(out)
R1=setscores(C['R1 (D1 post-pre)']).dropna(); R2=setscores(C['R2 (D2 post-pre)']).dropna(); REC=setscores(C['D1post->D2pre']).dropna()
ids=R1.index.intersection(R2.index)
R1,R2=R1.loc[ids],R2.loc[ids]
# 1. per-set reproducibility of individual response (R1 vs R2 Spearman across people)
rep=pd.DataFrame({s:stats.spearmanr(R1[s],R2[s]) for s in R1.columns},index=['rho','p']).T.sort_values('rho',ascending=False)
rep.to_csv(OUT+'response_reproducibility_R1_vs_R2.csv')
print('Median day1-vs-day2 response rho across sets: %.2f; sets with rho>0.5: %d of %d'%(rep.rho.median(),(rep.rho>0.5).sum(),len(rep)))
print(rep.head(12).round(3).to_string()); print(rep.loc[[i for i in rep.index if i.startswith('LEAD')]].round(3).to_string())
# 2. metabolite-level: self vs other profile correlation
A=C['R1 (D1 post-pre)'].loc[ids]; B=C['R2 (D2 post-pre)'].loc[ids]
keep=A.columns[(A.std()>0.05)&(B.std()>0.05)]
A=((A[keep]-A[keep].mean())/A[keep].std()).fillna(0); B=((B[keep]-B[keep].mean())/B[keep].std()).fillna(0)
M=np.corrcoef(A.values,B.values)[:len(ids),len(ids):]
self_=np.diag(M); other=M[~np.eye(len(ids),dtype=bool)]
print('\nWhole-profile day1 vs day2 response, same person r=%.3f (IQR %.2f-%.2f) vs different person r=%.3f'%(self_.mean(),*np.percentile(self_,[25,75]),other.mean()))
rank=[(M[i]>M[i,i]).sum() for i in range(len(ids))]; print('fraction of people whose own day2 profile is their best match: %.2f'%np.mean(np.array(rank)==0))
# 3. clustering patients by exercise response (R1 set scores), stability
pat=[i for i in ids if meta.loc[i,'case']==1]
def clus(D,k,seed=0):
    Z=(D-D.mean())/D.std(); P=PCA(n_components=8,random_state=0).fit_transform(Z.fillna(0))
    km=KMeans(k,n_init=20,random_state=seed).fit(P); return km.labels_,P
rng=np.random.default_rng(0)
for k in [2,3,4]:
    lab1,P1=clus(R1.loc[pat],k); lab2,P2=clus(R2.loc[pat],k)
    sil=silhouette_score(P1,lab1)
    # null silhouette: columnwise permuted data
    nulls=[]
    for b in range(50):
        Dp=R1.loc[pat].apply(lambda c: rng.permutation(c.values)); l,Pp=clus(Dp,k); nulls.append(silhouette_score(Pp,l))
    # bootstrap stability
    aris=[]
    for b in range(100):
        bi=rng.choice(len(pat),len(pat),replace=True); ub=np.unique(bi)
        lb,_=clus(R1.loc[pat].iloc[bi],k,seed=b)
        # map back: compare labels on unique sampled people
        lbu=pd.Series(lb,index=bi).groupby(level=0).first()
        aris.append(adjusted_rand_score(lab1[lbu.index],lbu.values))
    print(f'k={k}: silhouette {sil:.3f} (null {np.mean(nulls):.3f}, 95% {np.percentile(nulls,95):.3f}); bootstrap ARI median {np.median(aris):.2f}; day1-clusters vs day2-clusters ARI {adjusted_rand_score(lab1,lab2):.2f}; sizes {np.bincount(lab1)}')
# 4. simple 1-D response phenotypes: lactate/glycolysis rise, relation to Bell and sex
d=meta.loc[pat].copy(); d['lac1']=R1.loc[pat,'LEAD lactate']; d['lac2']=R2.loc[pat,'LEAD lactate']
print('\nlactate response rho day1 vs day2 patients: %.2f'%stats.spearmanr(d.lac1,d.lac2)[0])
print('lactate rise vs Bell: day1 rho %.2f p %.3f; day2 rho %.2f p %.3f'%(*stats.spearmanr(d.Bell,d.lac1),*stats.spearmanr(d.Bell,d.lac2)))
c=meta.loc[[i for i in ids if meta.loc[i,'case']==0]]; 
print('controls lactate rise day1 median %.2f; patients %.2f; patients below control 10th pct: %d/%d'%(R1.loc[c.index,'LEAD lactate'].median(),d.lac1.median(),(d.lac1<R1.loc[c.index,'LEAD lactate'].quantile(.1)).sum(),len(d)))
# 5. what drives k=2 split
lab1,P1=clus(R1.loc[pat],2)
diff=(R1.loc[pat][lab1==1].mean()-R1.loc[pat][lab1==0].mean()).sort_values()
print('\nk=2 split, largest set differences:'); print(diff.head(6).round(2).to_string()); print(diff.tail(6).round(2).to_string())
g=pd.Series(lab1,index=pat); dm=meta.loc[pat]
print('cluster vs sex:',pd.crosstab(g,dm.Sex).values.tolist(),' Bell by cluster:',dm.groupby(g).Bell.median().to_dict(),' MWU p=%.3f'%stats.mannwhitneyu(dm.Bell[g==0],dm.Bell[g==1]).pvalue)
ctl=[i for i in ids if meta.loc[i,'case']==0]; lc,Pc=clus(R1.loc[ctl],2)
print('controls k=2 silhouette %.3f sizes %s'%(silhouette_score(Pc,lc),np.bincount(lc)))
