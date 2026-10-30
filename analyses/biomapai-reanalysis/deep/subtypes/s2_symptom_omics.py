"""Part 1b: which omics features differ between symptom-based partitions (k=2,3) and track symptom PCs; OLS adjusted for age, sex, BMI; BH-FDR."""
from sub_common import *
from clust_tools import *
from sklearn.decomposition import PCA
P=meta[(meta.tp==1)&(meta.y==1)].join(score,how='inner')
S=P[SC]; Z=((S-S.mean())/S.std()).values
pca=PCA(4).fit(Z); PCs=pca.transform(Z)
if np.corrcoef(PCs[:,0],P.severity)[0,1]<0: PCs[:,0]*=-1
load_=pd.DataFrame(pca.components_.T,index=[SHORT[c] for c in SC],columns=['PC1','PC2','PC3','PC4'])
print(load_.round(2).to_string()); load_.to_csv(OUT+'results/symptom_pca_loadings.csv')
lab={k:fit(Z,k,'kmeans') for k in (2,3)}
# order clusters by severity
for k in lab:
    order=pd.Series(P.severity.values).groupby(lab[k]).mean().sort_values().index
    lab[k]=pd.Series(lab[k]).map({o:i for i,o in enumerate(order)}).values
P['sym_k2']=lab[2]; P['sym_k3']=lab[3]
for i in range(4): P[f'symPC{i+1}']=PCs[:,i]
P[['sym_k2','sym_k3','symPC1','symPC2','symPC3','symPC4','severity']].to_csv(OUT+'results/symptom_assignments_tp1.csv')
print(pd.DataFrame(Z,columns=[SHORT[c] for c in SC]).groupby(lab[3]).mean().round(2).T.to_string())
print('k3 sizes',np.bincount(lab[3]),'severity by k3',P.groupby('sym_k3').severity.mean().round(3).to_dict())
def ftest(Y,cov,G):
    """Y n x p; cov n x c; G n x g design of the tested term. returns F p-values and partial R2"""
    D0=np.column_stack([np.ones(len(Y)),cov]); D1=np.column_stack([D0,G])
    r0=Y-D0@np.linalg.lstsq(D0,Y,rcond=None)[0]; r1=Y-D1@np.linalg.lstsq(D1,Y,rcond=None)[0]
    ss0=(r0**2).sum(0); ss1=(r1**2).sum(0); df1=G.shape[1]; df2=len(Y)-D1.shape[1]
    F=((ss0-ss1)/df1)/(ss1/df2); return stats.f.sf(F,df1,df2),1-ss1/ss0
rows=[]
for L in ['immune','specie','metabolome','quest']:
    X=layer(L); idx=P.index.intersection(X.index); Y=X.loc[idx].values; M=P.loc[idx]
    cov=M[['age','male','bmi_i','bmi_miss']].values.astype(float)
    for term in ['sym_k2','sym_k3','symPC1','symPC2','symPC3']:
        G=pd.get_dummies(M[term],drop_first=True).values.astype(float) if term.startswith('sym_k') else M[[term]].values
        p,r2=ftest(Y,cov,G)
        # direction: mean diff highest-vs-lowest cluster (z units) or correlation sign
        if term.startswith('sym_k'): eff=Y[M[term].values==M[term].max()].mean(0)-Y[M[term].values==0].mean(0)
        else: eff=np.array([stats.spearmanr(Y[:,j],M[term])[0] for j in range(Y.shape[1])])
        for j,f in enumerate(X.columns): rows.append((L,term,f,len(idx),p[j],r2[j],eff[j]))
R=pd.DataFrame(rows,columns=['layer','term','feature','n','p','partial_R2','effect'])
R['q_layer']=R.groupby(['layer','term']).p.transform(bh); R['q_all']=R.groupby('term').p.transform(bh)
R.to_csv(OUT+'results/symptom_partition_omics_tests.csv',index=False)
print(R.groupby(['term','layer']).apply(lambda d:pd.Series(dict(n_feat=len(d),q05=(d.q_layer<0.05).sum(),q10=(d.q_layer<0.1).sum(),p01=(d.p<0.01).sum(),exp01=0.01*len(d)))).to_string())
pd.set_option('display.width',220)
print(R[R.q_layer<0.1].sort_values('p').head(40).round(4).to_string())
