import numpy as np, pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score as ARI, silhouette_score
def fit(X,k,method,seed=0):
    if method=='kmeans': return KMeans(k,n_init=10,random_state=seed).fit_predict(X)
    if method=='ward': return AgglomerativeClustering(k,linkage='ward').fit_predict(X)
    if method=='gmm': return GaussianMixture(k,covariance_type='diag',n_init=3,random_state=seed).fit(X).predict(X)
def stability(X,k,method,B=50,frac=0.8,seed=0):
    """mean ARI between clustering of an 80% subsample and the full-data clustering restricted to that subsample,
    plus mean pairwise ARI between subsamples on their overlap"""
    rng=np.random.default_rng(seed); n=len(X); full=fit(X,k,method,seed)
    labs=[]; a=[]
    for b in range(B):
        idx=np.sort(rng.choice(n,int(frac*n),replace=False))
        l=fit(X[idx],k,method,seed+b+1); a.append(ARI(full[idx],l))
        v=np.full(n,-1); v[idx]=l; labs.append(v)
    pw=[]
    for i in range(0,B,2):
        j=i+1; m=(labs[i]>=0)&(labs[j]>=0); pw.append(ARI(labs[i][m],labs[j][m]))
    return np.mean(a),np.mean(pw),full
def consensus(X,k,method='kmeans',B=60,frac=0.8,seed=0):
    rng=np.random.default_rng(seed); n=len(X); M=np.zeros((n,n)); I=np.zeros((n,n))
    for b in range(B):
        idx=np.sort(rng.choice(n,int(frac*n),replace=False)); l=fit(X[idx],k,method,seed+b+1)
        same=(l[:,None]==l[None,:]).astype(float); M[np.ix_(idx,idx)]+=same; I[np.ix_(idx,idx)]+=1
    C=np.divide(M,I,out=np.zeros_like(M),where=I>0)
    iu=np.triu_indices(n,1); pac=np.mean((C[iu]>0.1)&(C[iu]<0.9))
    return C,pac
def gaussian_null(X,rng):
    """same mean/covariance, no clusters"""
    return rng.multivariate_normal(X.mean(0),np.cov(X,rowvar=False),size=len(X))
