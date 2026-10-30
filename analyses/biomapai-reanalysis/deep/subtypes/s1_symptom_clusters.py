"""Part 1a: cluster tp1 patients on 12 standardised symptom scores; choose k by stability vs a Gaussian (no-cluster) null."""
from sub_common import *
from clust_tools import *
from sklearn.decomposition import PCA
P=meta[(meta.tp==1)&(meta.y==1)].join(score,how='inner')
S=P[SC]; mu,sd=S.mean(),S.std(); Z=((S-mu)/sd).values
shape=Z-Z.mean(1,keepdims=True)  # profile shape: remove each person's overall level
pca=PCA().fit(Z); print('symptom PC var explained',pca.explained_variance_ratio_[:4].round(3))
pc1=pca.transform(Z)[:,0]; print('PC1 vs severity rho',round(stats.spearmanr(pc1,P.severity)[0],3))
rows=[]; rng=np.random.default_rng(1)
for space,X in [('level',Z),('shape',shape)]:
    nulls=[gaussian_null(X,rng) for _ in range(20)]
    for method in ['kmeans','ward','gmm']:
        for k in range(2,7):
            s1,s2,lab=stability(X,k,method,B=40)
            sil=silhouette_score(X,lab)
            ns=[];nst=[]
            for N in nulls[:10]:
                l=fit(N,k,method); ns.append(silhouette_score(N,l))
                if method=='kmeans': nst.append(stability(N,k,method,B=10)[0])
            rows.append(dict(space=space,method=method,k=k,ARI_vs_full=s1,ARI_pairwise=s2,silhouette=sil,
                 null_sil_mean=np.mean(ns),null_sil_95=np.quantile(ns,0.95),null_stab=np.mean(nst) if nst else np.nan,
                 min_cluster=np.bincount(lab).min()))
    # GMM BIC k=1..6 (k=1 = single continuous Gaussian)
    for cov in ['full','diag']:
        b=[GaussianMixture(k,covariance_type=cov,n_init=3,random_state=0).fit(X).bic(X) for k in range(1,7)]
        print(space,cov,'GMM BIC k=1..6',np.round(b,0),'best k',int(np.argmin(b))+1)
R=pd.DataFrame(rows); R.to_csv(OUT+'results/symptom_cluster_selection.csv',index=False)
pd.set_option('display.width',200); print(R.round(3).to_string())
