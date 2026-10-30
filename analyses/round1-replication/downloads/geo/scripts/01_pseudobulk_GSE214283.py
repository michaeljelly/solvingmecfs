# Pseudobulk + crude marker-based cell typing for GSE214283 (10x PBMC, 30 ME/CFS vs 28 controls, D1 pre-CPET and D2 24h post)
import h5py, numpy as np, pandas as pd, scipy.sparse as sp, glob, os, re, sys
D='/tmp/work_geo/sc'; OUT='/tmp/work_geo/pb'; os.makedirs(OUT,exist_ok=True)
files=sorted(glob.glob(D+'/*.h5'))
MARK={'T':['CD3D','CD3E','CD3G'],'NK':['NKG7','GNLY','KLRD1','KLRF1','PRF1'],'B':['MS4A1','CD79A','CD79B','BANK1'],
      'Mono':['LYZ','S100A8','S100A9','FCN1','CST3'],'pDC':['LILRA4','IL3RA','CLEC4C'],'Platelet':['PPBP','PF4'],'Eryth':['HBB','HBA1','HBA2']}
def load(f):
    h=h5py.File(f,'r')['matrix']
    X=sp.csc_matrix((h['data'][:],h['indices'][:],h['indptr'][:]),shape=h['shape'][:])  # genes x cells
    feat=h['features']; names=[x.decode() for x in feat['name'][:]]; ftype=[x.decode() for x in feat['feature_type'][:]]
    keep=np.array([t=='Gene Expression' for t in ftype])
    X=X[keep,:].tocsr(); names=np.array(names)[keep]
    return X,names
res={}; props=[]
genes_ref=None
for i,f in enumerate(files):
    gsm,title=re.match(r'.*/(GSM\d+)_(.*)_filtered',f).groups()
    X,names=load(f)
    if genes_ref is None: genes_ref=names
    assert (names==genes_ref).all()
    lib=np.asarray(X.sum(0)).ravel(); ngen=np.asarray((X>0).sum(0)).ravel()
    mt=np.array([n.startswith('MT-') for n in names]); pmt=np.asarray(X[mt].sum(0)).ravel()/np.maximum(lib,1)
    ok=(ngen>=200)&(pmt<0.2)&(lib>=500)
    Xn=X[:,ok].multiply(1e4/lib[ok]).tocsr(); Xn.data=np.log1p(Xn.data)
    idx={n:j for j,n in enumerate(names)}
    sc={}
    for k,g in MARK.items():
        rows=[idx[x] for x in g if x in idx]; sc[k]=np.asarray(Xn[rows].mean(0)).ravel()
    S=pd.DataFrame(sc); lab=S.idxmax(1).values; lab[S.max(1).values<0.3]='Unassigned'
    cd8=np.asarray(Xn[[idx['CD8A'],idx['CD8B']]].mean(0)).ravel()
    cd4like=np.asarray(Xn[[idx[x] for x in ['CD4','IL7R','CCR7','LTB'] if x in idx]].mean(0)).ravel()
    lab=np.where((lab=='T')&(cd8>0.3),'CD8T',np.where(lab=='T','CD4T_other',lab))
    mono16=np.asarray(Xn[[idx['FCGR3A'],idx['MS4A7']]].mean(0)).ravel(); cd14=np.asarray(Xn[[idx['CD14']]].mean(0)).ravel()
    lab=np.where((lab=='Mono')&(mono16>cd14),'ncMono',np.where(lab=='Mono','cMono',lab))
    Xok=X[:,ok]
    for ct in ['ALL','CD4T_other','CD8T','NK','B','cMono','ncMono']:
        m=np.ones(ok.sum(),bool) if ct=='ALL' else (lab==ct)
        res.setdefault(ct,{})[gsm]=np.asarray(Xok[:,m].sum(1)).ravel()
    vc=pd.Series(lab).value_counts(); vc['n_cells']=ok.sum(); vc['gsm']=gsm; vc['title']=title; props.append(vc)
    print(i,gsm,title,ok.sum(),dict(pd.Series(lab).value_counts(normalize=True).round(3)),flush=True)
for ct,d in res.items():
    df=pd.DataFrame(d,index=genes_ref); df=df.groupby(level=0).sum()
    df.to_parquet(f'{OUT}/pb_{ct}.parquet')
pd.DataFrame(props).fillna(0).to_csv(f'{OUT}/celltype_counts.csv',index=False)
