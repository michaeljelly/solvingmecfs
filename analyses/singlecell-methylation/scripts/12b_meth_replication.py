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
out={g:pd.read_parquet(f'{W}/ewas_{g}.parquet') for g in CFG}
import sys
MOD=sys.argv[1] if len(sys.argv)>1 else 'C'
for g in out:
    out[g]['t_B']=out[g][f't_{MOD}']; out[g]['p_B']=out[g][f'p_{MOD}']

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
rep=pd.DataFrame(rep); rep.to_csv(f'{R}/meth_ewas_replication_top1000_model'+MOD+'.csv',index=False); print(rep.round(3).to_string())
# Stouffer meta across blood datasets (model B), weights sqrt(n)
Z=pd.DataFrame({g:np.sign(out[g].loc[common,'t_B'])*stats.norm.isf(out[g].loc[common,'p_B']/2) for g in blood})
w=np.array([np.sqrt(out[g]['n'].iloc[0]) for g in blood]); mz=(Z.values*w).sum(1)/np.sqrt((w**2).sum())
meta=pd.DataFrame({'z_meta':mz},index=common); meta['p_meta']=2*stats.norm.sf(np.abs(mz)); meta['q_meta']=multipletests(meta.p_meta,method='fdr_bh')[1]
meta['n_same_sign']=(np.sign(Z).eq(np.sign(mz),axis=0)).sum(1)
for g in blood: meta[f'z_{g}']=Z[g]; meta[f'dbeta_{g}']=out[g].loc[common,'dbeta']
sal=out['GSE327255']; meta['z_saliva']=(np.sign(sal['t_B'])*stats.norm.isf(sal['p_B']/2)).reindex(common)
a=ann450.reindex(common); meta['CHR']=a.CHR; meta['POS']=a.MAPINFO; meta['gene']=a.UCSC_RefGene_Name; meta['region']=a.UCSC_RefGene_Group; meta['island']=a.Relation_to_UCSC_CpG_Island
meta=meta.sort_values('p_meta'); meta.head(5000).to_csv(f'{R}/meth_ewas_meta_top5000_model'+MOD+'.csv')
print('meta probes',len(meta),'q<0.05',(meta.q_meta<0.05).sum(),'q<0.05 & all 3 same sign',((meta.q_meta<0.05)&(meta.n_same_sign==3)).sum())
print(meta.head(25)[['z_meta','p_meta','q_meta','n_same_sign','z_GSE156792','z_GSE93266','z_GSE111183','z_saliva','gene','CHR']].round(4).to_string())
meta.to_parquet(f'{W}/ewas_meta_'+MOD+'.parquet')
