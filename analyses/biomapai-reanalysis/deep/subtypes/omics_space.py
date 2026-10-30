"""Builds the patients-only multi-omics embedding (fit on tp1 patients) and can project any other samples into it."""
from sub_common import *
from sklearn.decomposition import PCA
LAYERS=['immune','specie','metabolome','quest']
NPC=10
COV=['age','male','bmi_i','bmi_miss']
class OmicsSpace:
    def __init__(self,layers=LAYERS,npc=NPC):
        self.layers=layers; self.npc=npc; self.X={L:layer(L) for L in layers}
        common=None
        for L in layers: common=self.X[L].index if common is None else common.intersection(self.X[L].index)
        self.all_idx=common
        self.pt=[i for i in common if meta.loc[i,'tp']==1 and meta.loc[i,'y']==1]
        self.models={}
        for L in layers:
            Xp=self.X[L].loc[self.pt]; Mp=meta.loc[self.pt]
            D=np.column_stack([np.ones(len(Mp)),Mp[COV].values.astype(float)])
            beta=np.linalg.lstsq(D,Xp.values,rcond=None)[0]
            R=Xp.values-D@beta; mu=R.mean(0); sd=R.std(0)+1e-9
            pca=PCA(npc,random_state=0).fit((R-mu)/sd)
            T=pca.transform((R-mu)/sd); w=np.sqrt((T**2).sum(1).mean())  # block weight: equal total variance per layer
            self.models[L]=(beta,mu,sd,pca,w)
    def embed(self,idx,per_layer=False):
        M=meta.loc[idx]; D=np.column_stack([np.ones(len(M)),M[COV].values.astype(float)]); out={}
        for L in self.layers:
            beta,mu,sd,pca,w=self.models[L]
            out[L]=pca.transform(((self.X[L].loc[idx].values-D@beta)-mu)/sd)/w
        return out if per_layer else np.hstack([out[L] for L in self.layers])
