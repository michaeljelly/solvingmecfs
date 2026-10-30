# Person-level pseudobulk per fine cell type; DE baseline (D1) and 24h change (D2-D1), case vs control, sex+age bin; BH per cell type.
import pandas as pd, numpy as np, scipy.sparse as sp, glob, sys
sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/singlecell_methyl/scripts')
from scipy import stats
from statsmodels.stats.multitest import multipletests
O='/tmp/work_scm/sc_out'; R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
src=open('/mnt/project-files/mecfs-research/deep2/singlecell_methyl/scripts/02_sc_collect.py').read(); exec(src[src.index('MAP={'):src.index("f=C['fine']")])
genes=pd.read_csv(f'{O}/genes.csv').iloc[:,0].values
meta=pd.read_csv(f'{R}/sc_sample_meta.csv').set_index('geo_accession'); meta['case']=(meta['disease state']=='case').astype(int); meta['age']=meta['age bin'].astype(float)
PB={}
for gsm in meta.index:
    M=sp.load_npz(f'{O}/{gsm}_pb.npz').toarray(); idx=pd.read_csv(f'{O}/{gsm}_pbidx.csv')
    for i,r in idx.iterrows():
        t=r.ct if r.ct.startswith(('Mono_','NK_')) or r.ct=='Doublet' else MAP.get(r.ct,'Other')
        PB.setdefault(t,{}).setdefault(gsm,[np.zeros(len(genes)),0]); PB[t][gsm][0]+=M[i]; PB[t][gsm][1]+=r.n
def ols(Y,X):
    XtX=np.linalg.pinv(X.T@X); Bh=Y@X@XtX.T; res=Y-Bh@X.T; df=X.shape[0]-X.shape[1]; s2=(res**2).sum(1)/df
    se=np.sqrt(s2*XtX[1,1]); t=Bh[:,1]/se; return Bh[:,1],t,2*stats.t.sf(np.abs(t),df)
summary=[]; allres=[]
for t,d in PB.items():
    if t in('Doublet','Other','Erythroid','HSPC','cDC1','NK_CD56bright'): continue
    n=pd.Series({g:v[1] for g,v in d.items()})
    ok=n[n>=20].index
    if len(ok)<60: continue
    Cn=pd.DataFrame({g:d[g][0] for g in ok},index=genes); Cn=Cn.groupby(level=0).sum()
    lib=Cn.sum(0); cpm=Cn/lib*1e6; keepg=(cpm>=5).mean(1)>=0.5; L=np.log2(cpm[keepg]+1)
    np.save if False else None
    L.astype('float32').to_parquet(f'/tmp/work_scm/pb_{t.replace("/","_")}.parquet')
    m=meta.loc[ok]; m.index.name='geo_accession'
    # baseline
    b=m[m.day=='D1']; X=np.column_stack([np.ones(len(b)),b.case,(b.Sex=='M').astype(float),b.age,np.log(n[b.index])])
    lf,tt,p=ols(L[b.index].values,X); q=multipletests(p,method='fdr_bh')[1]
    r1=pd.DataFrame({'celltype':t,'contrast':'baseline','gene':L.index,'log2FC':lf,'t':tt,'p':p,'q':q})
    # change
    w1=m[m.day=='D1'].reset_index().set_index('person'); w2=m[m.day=='D2'].reset_index().set_index('person'); pp=w1.index.intersection(w2.index)
    Dl=L[w2.loc[pp,'geo_accession']].values-L[w1.loc[pp,'geo_accession']].values
    X=np.column_stack([np.ones(len(pp)),w1.loc[pp,'case'],(w1.loc[pp,'Sex']=='M').astype(float),w1.loc[pp,'age']])
    lf2,tt2,p2=ols(Dl,X); q2=multipletests(p2,method='fdr_bh')[1]
    r2=pd.DataFrame({'celltype':t,'contrast':'change_D2minusD1','gene':L.index,'log2FC':lf2,'t':tt2,'p':p2,'q':q2})
    allres+= [r1,r2]
    lam=lambda p:np.median(stats.chi2.isf(p,1))/0.4549
    summary.append(dict(celltype=t,n_samples=len(ok),n_base=len(b),n_pairs=len(pp),median_cells=int(n[ok].median()),genes=int(keepg.sum()),q05_baseline=int((q<0.05).sum()),q10_baseline=int((q<0.1).sum()),lambda_base=lam(p),
                        q05_change=int((q2<0.05).sum()),q10_change=int((q2<0.1).sum()),lambda_change=lam(p2)))
    print(summary[-1],flush=True)
A=pd.concat(allres); A.to_parquet('/tmp/work_scm/sc_de_all.parquet')
A[A.q<0.1].sort_values('p').to_csv(f'{R}/sc_pseudobulk_DE_q10.csv',index=False)
A.sort_values('p').groupby(['celltype','contrast']).head(30).to_csv(f'{R}/sc_pseudobulk_DE_top30_per_celltype.csv',index=False)
pd.DataFrame(summary).to_csv(f'{R}/sc_pseudobulk_DE_summary.csv',index=False)
print(A[A.q<0.1].sort_values('p').head(60).to_string())
