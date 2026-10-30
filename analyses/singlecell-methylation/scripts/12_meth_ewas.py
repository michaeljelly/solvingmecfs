# EWAS per dataset (M-values, OLS), with and without cell-composition covariates; cross-dataset replication + Stouffer meta; simple DMRs; targeted genes.
import pandas as pd, numpy as np, sys
from scipy import stats
from statsmodels.stats.multitest import multipletests
W='/tmp/work_scm/meth'; R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
ann450=pd.read_parquet(f'{W}/ann450.parquet').set_index('IlmnID'); annE=pd.read_parquet(f'{W}/annEPIC.parquet').set_index('IlmnID')
snp=set(ann450.index[ann450.Probe_SNPs.notna()|ann450.Probe_SNPs_10.notna()])
CFG={'GSE156792':dict(ann=ann450,base=['sex_inf','chip'],comp=['frac_CD4nv','frac_CD8nv','frac_CD8mem','frac_NK']),
     'GSE93266':dict(ann=ann450,base=[],comp=['rein_neutrophil','rein_monocyte','rein_b_cell','rein_nk_cell','rein_cd8_t_cell']),
     'GSE111183':dict(ann=annE,base=[],comp=['rein_neutrophil','rein_monocyte','rein_b_cell','rein_nk_cell','rein_cd8_t_cell']),
     'GSE327255':dict(ann=annE,base=['sex_inf','age'],comp=['frac_Neu','frac_Mono'])}
def design(D,cols):
    X=[np.ones(len(D)),(D.group=='ME').astype(float).values]
    for c in cols:
        v=D[c]
        if v.dtype==object:
            X+= [ (v==lvl).astype(float).values for lvl in sorted(v.unique())[1:] ]
        else: X.append(v.astype(float).values)
    return np.column_stack(X)
def ols_t(Y,X):  # Y: probes x samples
    XtX=np.linalg.pinv(X.T@X); Bh=Y@X@XtX.T  # probes x k
    res=Y-Bh@X.T; df=X.shape[0]-np.linalg.matrix_rank(X); s2=(res**2).sum(1)/df
    se=np.sqrt(s2*XtX[1,1]); t=Bh[:,1]/se; p=2*stats.t.sf(np.abs(t),df)
    return Bh[:,1],t,p,df
out={}
for g,cfg in CFG.items():
    B=pd.read_parquet(f'{W}/{g}_beta.parquet'); D=pd.read_csv(f'{R}/meth_{g}_sample_scores.csv',index_col=0)
    B=B[D.index]; ann=cfg['ann']
    if 'chip' in D: D['chip']=D['chip'].astype(str)
    keep=ann.index[(~ann.CHR.isin(['X','Y']))&ann.index.str.startswith('cg')].intersection(B.index)
    keep=keep[~keep.isin(snp)]
    B=B.loc[keep]; B=B[B.notna().mean(1)>0.9]
    B=B.apply(lambda r:r.fillna(r.mean()),axis=1)
    Bv=B.values.astype(np.float64).clip(0.001,0.999); M=np.log2(Bv/(1-Bv))
    dbeta=Bv[:,(D.group=='ME').values].mean(1)-Bv[:,(D.group=='HC').values].mean(1)
    base=[c for c in cfg['base'] if D[c].nunique()>1]
    res=pd.DataFrame(index=B.index); res['dbeta']=dbeta
    rs=np.random.default_rng(0).choice(M.shape[0],min(50000,M.shape[0]),replace=False); Ms=M[rs]-M[rs].mean(1,keepdims=True)
    u,sv,vt=np.linalg.svd(Ms,full_matrices=False)
    for k in range(4): D[f'PC{k+1}']=vt[k]
    from scipy.stats import mannwhitneyu
    print(g,'PC var',np.round(sv[:4]**2/(sv**2).sum(),3),'PC-group MWU p',[round(mannwhitneyu(D.loc[D.group=='ME',f'PC{k+1}'],D.loc[D.group=='HC',f'PC{k+1}']).pvalue,4) for k in range(4)],
          'mean beta ME %.4f HC %.4f'%(Bv[:,(D.group=='ME').values].mean(),Bv[:,(D.group=='HC').values].mean()))
    for lab,cols in [('A',base),('B',base+cfg['comp']),('C',base+cfg['comp']+['PC1','PC2','PC3','PC4'])]:
        X=design(D,cols); b,t,p,df=ols_t(M,X)
        res[f't_{lab}']=t; res[f'p_{lab}']=p; res[f'q_{lab}']=multipletests(p,method='fdr_bh')[1]
        lam=np.median(stats.chi2.isf(p,1))/0.4549
        print(g,lab,'n',len(D),'df',df,'probes',len(p),'FDR<0.05:',(res[f'q_{lab}']<0.05).sum(),'lambda %.2f'%lam,flush=True)
    res['n']=len(D); res.to_parquet(f'{W}/ewas_{g}.parquet'); out[g]=res
    del B,Bv,M
# replication among blood datasets (model B)
blood=['GSE156792','GSE93266','GSE111183']
common=out[blood[0]].index
for g in blood[1:]: common=common.intersection(out[g].index)
rep=[]
for dg in blood+['GSE327255']:
    top=out[dg].sort_values('p_B').head(1000)
    for vg in blood+['GSE327255']:
        if vg==dg: continue
        s=top.index.intersection(out[vg].index)
        conc=(np.sign(top.loc[s,'t_B'])==np.sign(out[vg].loc[s,'t_B'])).mean()
        nom=((np.sign(top.loc[s,'t_B'])==np.sign(out[vg].loc[s,'t_B']))&(out[vg].loc[s,'p_B']<0.05)).sum()
        rep.append(dict(discovery=dg,validation=vg,n_top=len(s),sign_concordance=conc,binom_p=stats.binomtest(int(round(conc*len(s))),len(s)).pvalue,n_same_sign_p05=nom,expected_p05=len(s)*0.025))
rep=pd.DataFrame(rep); rep.to_csv(f'{R}/meth_ewas_replication_top1000.csv',index=False); print(rep.round(3).to_string())
# Stouffer meta across blood datasets (model B), weights sqrt(n)
Z=pd.DataFrame({g:np.sign(out[g].loc[common,'t_B'])*stats.norm.isf(out[g].loc[common,'p_B']/2) for g in blood})
w=np.array([np.sqrt(out[g]['n'].iloc[0]) for g in blood]); mz=(Z.values*w).sum(1)/np.sqrt((w**2).sum())
meta=pd.DataFrame({'z_meta':mz},index=common); meta['p_meta']=2*stats.norm.sf(np.abs(mz)); meta['q_meta']=multipletests(meta.p_meta,method='fdr_bh')[1]
meta['n_same_sign']=(np.sign(Z).eq(np.sign(mz),axis=0)).sum(1)
for g in blood: meta[f'z_{g}']=Z[g]; meta[f'dbeta_{g}']=out[g].loc[common,'dbeta']
sal=out['GSE327255']; meta['z_saliva']=(np.sign(sal['t_B'])*stats.norm.isf(sal['p_B']/2)).reindex(common)
a=ann450.reindex(common); meta['CHR']=a.CHR; meta['POS']=a.MAPINFO; meta['gene']=a.UCSC_RefGene_Name; meta['region']=a.UCSC_RefGene_Group; meta['island']=a.Relation_to_UCSC_CpG_Island
meta=meta.sort_values('p_meta'); meta.head(5000).to_csv(f'{R}/meth_ewas_meta_top5000.csv')
print('meta probes',len(meta),'q<0.05',(meta.q_meta<0.05).sum(),'q<0.05 & all 3 same sign',((meta.q_meta<0.05)&(meta.n_same_sign==3)).sum())
print(meta.head(25)[['z_meta','p_meta','q_meta','n_same_sign','z_GSE156792','z_GSE93266','z_GSE111183','z_saliva','gene','CHR']].round(4).to_string())
meta.to_parquet(f'{W}/ewas_meta.parquet')
# per-dataset top hits
for g in out:
    o=out[g].sort_values('p_B').head(2000); a=(ann450 if g in('GSE156792','GSE93266') else annE).reindex(o.index)
    o['gene']=a.UCSC_RefGene_Name; o['CHR']=a.CHR; o['POS']=a.MAPINFO; o.to_csv(f'{R}/meth_ewas_{g}_top2000.csv')
