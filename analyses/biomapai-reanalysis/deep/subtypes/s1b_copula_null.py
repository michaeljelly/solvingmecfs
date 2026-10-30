"""Stricter null: Gaussian copula keeping each score's empirical (discrete/skewed) marginal and the rank correlations, but no clusters."""
from sub_common import *
from clust_tools import *
P=meta[(meta.tp==1)&(meta.y==1)].join(score,how='inner')
def copula_null(S,rng):
    R=stats.norm.ppf((S.rank()-0.5)/len(S)); C=np.corrcoef(R,rowvar=False)
    G=rng.multivariate_normal(np.zeros(S.shape[1]),C,size=len(S)); U=stats.norm.cdf(G)
    return np.column_stack([np.quantile(S.iloc[:,j],U[:,j],method='inverted_cdf') for j in range(S.shape[1])])
rows=[]; rng=np.random.default_rng(7)
for label,cols in [('all12',SC),('no_rolelim10',[c for c in SC if 'Role' not in c])]:
    S=P[cols]; Z=((S-S.mean())/S.std()).values
    nulls=[copula_null(S,rng) for _ in range(20)]; nulls=[(N-N.mean(0))/N.std(0,ddof=1) for N in nulls]
    for k in range(2,6):
        s1,s2,lab=stability(Z,k,'kmeans',B=40); sil=silhouette_score(Z,lab)
        ns=[silhouette_score(N,fit(N,k,'kmeans')) for N in nulls]; nst=[stability(N,k,'kmeans',B=10)[0] for N in nulls[:10]]
        rows.append(dict(scores=label,k=k,stab=s1,null_stab_mean=np.mean(nst),null_stab_95=np.quantile(nst,.95),
            sil=sil,null_sil_mean=np.mean(ns),null_sil_95=np.quantile(ns,.95),p_sil=(1+sum(np.array(ns)>=sil))/(1+len(ns))))
R=pd.DataFrame(rows); R.to_csv(OUT+'results/symptom_cluster_copula_null.csv',index=False); print(R.round(3).to_string())
