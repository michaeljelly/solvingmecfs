# Per-sample QC, CellTypist (Immune_All_Low, majority voting) + marker gating, per-cell scores, pseudobulk by fine type.
import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['NUMBA_NUM_THREADS']='1'
import h5py, numpy as np, pandas as pd, scipy.sparse as sp, glob, re, sys, warnings
warnings.filterwarnings('ignore')
import anndata as ad, scanpy as sc, celltypist
D='/tmp/work_geo/sc'; OUT='/tmp/work_scm/sc_out'
files=sorted(glob.glob(D+'/*.h5'))
model=celltypist.models.Model.load('Immune_All_Low.pkl')
SETS={
 'ISG':['IFI44L','IFI6','ISG15','MX1','RSAD2','OAS1','IFIT1','IFIT3','XAF1','HERC5','IFI44','EPSTI1','OAS3','LY6E'],
 'hypoxia':['VEGFA','SLC2A1','PGK1','LDHA','BNIP3','ENO1','HK2','PDK1','ADM','EGLN3','ALDOA','P4HA1','ANKRD37','MXI1','NDRG1','BHLHE40','DDIT4'],
 'glycolysis':['HK1','HK2','GPI','PFKL','PFKP','ALDOA','TPI1','GAPDH','PGK1','PGAM1','ENO1','PKM','LDHA'],
 'oxphos':['NDUFA4','NDUFB8','SDHB','UQCRC1','UQCRC2','COX4I1','COX5A','ATP5F1A','ATP5F1B','ATP5MC1','CYCS','NDUFS1'],
 'FAO':['CPT1A','CPT2','ACADM','ACADVL','HADHA','HADHB','ECHS1','ACAA2','ETFA','ETFB'],
 'platelet':['PPBP','PF4','GP9','ITGA2B','TUBB1','GNG11','CLU','NRGN','SPARC','CAVIN2'],
 'IEG_stress':['FOS','FOSB','JUN','JUNB','DUSP1','IER2','EGR1','ATF3','KLF6','NR4A1'],
 'heatshock':['HSPA1A','HSPA1B','HSPH1','DNAJB1','HSP90AA1','HSPD1','HSPE1'],
 'ABC':['TBX21','ITGAX','FCRL5','ZEB2','FCRL3','CD86','SOX5','FGR','CR2_neg'],
 'cyto_effector':['KLRG1','CX3CR1','GZMB','GZMH','PRF1','FGFBP2','ASSP_placeholder','ZEB2','TBX21','S1PR5','B3GAT1'],
 'exhaustion':['PDCD1','TOX','HAVCR2','LAG3','TIGIT','CTLA4','ENTPD1','CXCL13','LAYN'],
 'Th17':['RORC','CCR6','KLRB1','IL23R','RORA','CCR4','IL17RE','CTSH'],
 'adaptiveNK':['KLRC2','KLRC3','CD3E_neg'],
 'inflam':['IL1B','TNF','CXCL8','CCL3','CCL4','NFKBIA','IL6','PTGS2','CCL3L1','CCL4L2'],
 'MHCII':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','CD74'],
 'AP1_TNF':['NFKB1','NFKB2','REL','RELB','TNFAIP3','NFKBIZ'],
 'adrenergic':['ADRB2','ADRB1','ADRA2A','ADRB3'],
 'kynurenine':['IDO1','KMO','KYNU','HAAO','QPRT','AFMID','TDO2'],
 'ribo':None,'mito':None,
}
def load(f):
    h=h5py.File(f,'r')['matrix']
    X=sp.csc_matrix((h['data'][:],h['indices'][:],h['indptr'][:]),shape=h['shape'][:])
    names=np.array([x.decode() for x in h['features']['name'][:]]); ft=np.array([x.decode() for x in h['features']['feature_type'][:]])
    ids=np.array([x.decode() for x in h['features']['id'][:]])
    bc=np.array([x.decode() for x in h['barcodes'][:]])
    k=ft=='Gene Expression'
    return X[k,:].T.tocsr().astype(np.float32), names[k], ids[k], bc
genes_ref=None
for i,f in enumerate(files):
    gsm,title=re.match(r'.*/(GSM\d+)_(.*)_filtered',f).groups()
    if os.path.exists(f'{OUT}/{gsm}_cells.parquet'): continue
    X,names,ids,bc=load(f)
    if genes_ref is None:
        genes_ref=names; pd.Series(names).to_csv(f'{OUT}/genes.csv',index=False)
    lib=np.asarray(X.sum(1)).ravel(); ngen=np.diff(X.indptr)
    mt=np.char.startswith(names.astype(str),'MT-'); rb=np.char.startswith(names.astype(str),'RPS')|np.char.startswith(names.astype(str),'RPL')
    pmt=np.asarray(X[:,mt].sum(1)).ravel()/np.maximum(lib,1); prb=np.asarray(X[:,rb].sum(1)).ravel()/np.maximum(lib,1)
    ok=(ngen>=200)&(pmt<0.2)&(lib>=500)&(ngen<6000)
    A=ad.AnnData(X[ok].copy()); A.var_names=pd.Index(names).astype(str); A.var_names_make_unique(); A.obs_names=bc[ok]
    A.obs['lib']=lib[ok]; A.obs['ngenes']=ngen[ok]; A.obs['pct_mt']=pmt[ok]; A.obs['pct_ribo']=prb[ok]
    try:
        sc.pp.scrublet(A, verbose=False); dbl=A.obs['predicted_doublet'].values.astype(bool); A.obs['dbl_score']=A.obs['doublet_score'].values
    except Exception as e:
        print('scrublet fail',e); dbl=np.zeros(A.n_obs,bool); A.obs['dbl_score']=np.nan
    A.layers['counts']=A.X.copy()
    sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
    pr=celltypist.annotate(A,model=model,majority_voting=True)
    lab=pr.predicted_labels
    A.obs['ct_pred']=lab['predicted_labels'].values; A.obs['ct_mv']=lab['majority_voting'].values
    A.obs['ct_conf']=pr.probability_matrix.max(1).values
    Xn=A.X.tocsc(); vi={g:j for j,g in enumerate(A.var_names)}
    def gexp(g): return np.asarray(Xn[:,vi[g]].todense()).ravel() if g in vi else np.zeros(A.n_obs)
    for k,gs in SETS.items():
        if gs is None: continue
        gg=[g for g in gs if g in vi]
        A.obs['s_'+k]=np.asarray(Xn[:,[vi[g] for g in gg]].mean(1)).ravel() if gg else np.nan
    for g in ['CD14','FCGR3A','KLRC2','FCER1G','NCAM1','KLRG1','CX3CR1','B3GAT1','PDCD1','TOX','HAVCR2','TBX21','ITGAX','FCRL5','CR2','FOXP3','RORC','CCR6','KLRB1','TRDC','TRGC1','SLC4A10','ZBTB16','GZMK','GZMB','IL7R','CD27','MKI67','PPBP','PF4','HBB','IFNG','TNF','IL17A','CD38','JCHAIN','LILRA4','HLA-DRA','ADRB2','SELL','CCR7','ZNF683','KLRC1']:
        A.obs['g_'+g]=gexp(g)
    # curated fine labels
    m=A.obs['ct_mv'].astype(str).values.copy()
    cd14=A.obs.g_CD14.values; cd16=A.obs.g_FCGR3A.values
    fine=m.copy()
    mono=np.isin(m,['Classical monocytes','Non-classical monocytes','Monocytes','Intermediate macrophages','Macrophages','Mono-mac'])
    fine[mono & (cd16<=0)]='Mono_classical'
    fine[mono & (cd16>0) & (cd14>0)]='Mono_intermediate'
    fine[mono & (cd16>0) & (cd14<=0)]='Mono_nonclassical'
    nk=np.isin(m,['CD16+ NK cells','NK cells'])
    fine[nk]='NK_CD56dim'
    fine[nk & (A.obs.g_KLRC2.values>0) & (A.obs.g_FCER1G.values<=0)]='NK_adaptive'
    fine[np.isin(m,['CD16- NK cells','Transitional NK'])]='NK_CD56bright'
    fine[dbl]='Doublet'
    A.obs['fine']=fine
    A.obs['doublet']=dbl
    obs=A.obs.copy()
    for c in obs.columns:
        if obs[c].dtype==np.float64: obs[c]=obs[c].astype(np.float32)
    obs['gsm']=gsm; obs['title']=title
    obs.to_parquet(f'{OUT}/{gsm}_cells.parquet')
    # pseudobulk counts (raw) by fine and by ct_mv
    C=A.layers['counts'].tocsr()
    out={}
    for lvl in ['fine']:
        for ct in pd.unique(obs[lvl]):
            idx=np.where(obs[lvl].values==ct)[0]
            out[ct]=(np.asarray(C[idx].sum(0)).ravel(), len(idx))
    cts=list(out.keys())
    M=sp.csr_matrix(np.vstack([out[c][0] for c in cts]).astype(np.float32))
    sp.save_npz(f'{OUT}/{gsm}_pb.npz',M); pd.DataFrame({'ct':cts,'n':[out[c][1] for c in cts]}).to_csv(f'{OUT}/{gsm}_pbidx.csv',index=False)
    print(i,gsm,title,A.n_obs,'dbl',dbl.sum(),obs.fine.value_counts().head(8).to_dict(),flush=True)
print('DONE')
