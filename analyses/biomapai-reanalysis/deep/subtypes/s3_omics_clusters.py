"""Part 2a: patients-only multi-omics clustering (residualised on age/sex/BMI; per-layer PCA, equal block weights).
Stability via subsampling + consensus PAC, compared with a Gaussian null of the same covariance; also SNF."""
from sub_common import *
from clust_tools import *
from omics_space import *
from sklearn.cluster import SpectralClustering
OS=OmicsSpace(); E=OS.embed(OS.pt); per=OS.embed(OS.pt,per_layer=True)
print('tp1 patients with all 4 layers',len(OS.pt),'embedding',E.shape)
for L in OS.layers: print(L,'10 PCs explain',OS.models[L][3].explained_variance_ratio_.sum().round(2),'PC1',OS.models[L][3].explained_variance_ratio_[0].round(3))
rng=np.random.default_rng(3); rows=[]
def run(space,X,methods=('kmeans','ward'),nnull=10):
    nulls=[gaussian_null(X,rng) for _ in range(nnull)]
    for method in methods:
        for k in range(2,7):
            s1,s2,lab=stability(X,k,method,B=30); sil=silhouette_score(X,lab)
            C,pac=consensus(X,k,method,B=40) if method=='kmeans' else (None,np.nan)
            ns=[];nst=[];npac=[]
            for N in nulls:
                l=fit(N,k,method); ns.append(silhouette_score(N,l))
            for N in nulls[:5]:
                nst.append(stability(N,k,method,B=10)[0])
                if method=='kmeans': npac.append(consensus(N,k,method,B=20)[1])
            rows.append(dict(space=space,method=method,k=k,stab=s1,null_stab=np.mean(nst),sil=sil,null_sil=np.mean(ns),null_sil95=np.quantile(ns,.95),
                 pac=pac,null_pac=np.mean(npac) if npac else np.nan,sizes=str(sorted(np.bincount(lab)))))
run('concat4',E)
for L in OS.layers: run('only_'+L,per[L]*OS.models[L][4],methods=('kmeans',),nnull=8)
# SNF
def affinity(X,K=20,mu=0.5):
    D=((X[:,None,:]-X[None,:,:])**2).sum(-1); Dn=np.sqrt(D)
    T=np.sort(Dn,1)[:,1:K+1].mean(1); eps=(T[:,None]+T[None,:]+Dn)/3
    return np.exp(-D/(2*(mu*eps)**2))
def snf(Xs,K=20,t=20):
    Ws=[affinity(X,K) for X in Xs]
    def norm(W): W=W.copy(); np.fill_diagonal(W,0); W=W/(2*W.sum(1,keepdims=True)); np.fill_diagonal(W,0.5); return W
    def knn(W):
        S=np.zeros_like(W); idx=np.argsort(-W,1)[:,:K]
        for i in range(len(W)): S[i,idx[i]]=W[i,idx[i]]
        return S/S.sum(1,keepdims=True)
    P=[norm(W) for W in Ws]; S=[knn(W) for W in Ws]
    for _ in range(t):
        P=[S[v]@(sum(P[u] for u in range(len(P)) if u!=v)/(len(P)-1))@S[v].T for v in range(len(P))]; P=[norm(p) for p in P]
    W=sum(P)/len(P); return (W+W.T)/2
Xs=[per[L] for L in OS.layers]
W=snf(Xs); n=len(W)
def spec(W,k,seed=0): return SpectralClustering(k,affinity='precomputed',random_state=seed).fit_predict(W)
for k in range(2,6):
    full=spec(W,k); a=[]
    for b in range(25):
        idx=np.sort(rng.choice(n,int(.8*n),replace=False)); Wb=snf([x[idx] for x in Xs]); a.append(ARI(full[idx],spec(Wb,k,b)))
    # null: joint gaussian of concatenated per-layer PCs
    na=[]
    for r in range(4):
        N=gaussian_null(E,rng); Ns=[N[:,i*10:(i+1)*10] for i in range(4)]; Wn=snf(Ns); fn=spec(Wn,k)
        for b in range(6):
            idx=np.sort(rng.choice(n,int(.8*n),replace=False)); na.append(ARI(fn[idx],spec(snf([x[idx] for x in Ns]),k,b)))
    rows.append(dict(space='SNF4',method='spectral',k=k,stab=np.mean(a),null_stab=np.mean(na),sizes=str(sorted(np.bincount(full)))))
R=pd.DataFrame(rows); R.to_csv(OUT+'results/omics_cluster_selection.csv',index=False)
pd.set_option('display.width',220); print(R.round(3).to_string())
np.save(OUT+'results/snf_W.npy',W)
