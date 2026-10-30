import numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy import stats
GENESETS={
 'Th17_IL17':['IL17A','IL17F','RORC','IL23R','CCR6','IL22','IL21','KLRB1'],
 'T_exhaustion':['PDCD1','HAVCR2','LAG3','TIGIT','CTLA4','TOX','ENTPD1','CD160','CD244'],
 'T_naive_CM':['CCR7','SELL','LEF1','TCF7','IL7R','CD27','CD28'],
 'T_effector_mem':['GZMB','GZMH','PRF1','KLRG1','CX3CR1','GZMK','FGFBP2'],
 'Kynurenine':['IDO1','IDO2','TDO2','KMO','KYNU','HAAO','QPRT','AFMID','KYAT1','AADAT'],
 'AHR_targets':['CYP1A1','CYP1B1','AHRR','TIPARP','NQO1'],
 'FAO':['CPT1A','CPT2','ACADVL','ACADM','HADHA','HADHB','SLC25A20','ETFDH','ACAA2','ECHS1'],
 'Ketone':['HMGCS2','OXCT1','BDH1','ACAT1','HMGCL','SLC16A1'],
 'PPARA_targets':['PDK4','ANGPTL4','CPT1A','ACADVL','PLIN2','ACOX1','FABP4','CD36','PPARA','TXNIP'],
 'OXPHOS_core':['NDUFA4','NDUFB8','NDUFS1','SDHA','SDHB','UQCRC1','UQCRC2','COX4I1','COX5A','ATP5F1A','ATP5F1B','CYCS'],
 'IFN_signature':['IFI27','IFI44L','IFIT1','ISG15','RSAD2','SIGLEC1','MX1','OAS1','OAS3','IFI6','IFI44','IFIT3'],
 'B_cell':['CD19','MS4A1','CD79A','CD79B','BANK1','PAX5','CD22','FCRL1','TCL1A','BLK'],
 'NK_cell':['NCR1','KLRF1','KLRD1','SH2D1B','NCAM1','KIR2DL4','KIR3DL1','GNLY','NKG7','PRF1'],
 'Adrenergic_receptors':['ADRB1','ADRB2','ADRB3','ADRA1A','ADRA1B','ADRA1D','ADRA2A','ADRA2B','ADRA2C'],
 'Catecholamine_metab':['TH','DDC','DBH','PNMT','COMT','MAOA','MAOB','SLC6A2','SLC6A3'],
 'cAMP_NR4A_response':['NR4A1','NR4A2','NR4A3','FOSB','CREM','PDE4B','RGS2','DUSP1'],
 'DecodeME_named':['OLFM4','RABGAP1L','FBXL4','CA10','BTN2A2','CSE1L','ARFGEF2','SUDS3','CCPG1'],
}
TARGET_GENES=sorted(set(g for v in GENESETS.values() for g in v)|{'AHR','ADRB2','IL17RA','IL23A','IL6','TNF','IL1B'})
def bh(p):
    p=np.asarray(p,float); out=np.full(p.shape,np.nan); ok=~np.isnan(p)
    if ok.sum(): out[ok]=multipletests(p[ok],method='fdr_bh')[1]
    return out
def logcpm(counts,prior=1.0):
    lib=counts.sum(0); return np.log2((counts+prior).div(lib+2*prior,axis=1)*1e6)
def filter_expr(logm,counts=None,min_cpm=1.0,min_frac=0.25):
    if counts is not None:
        cpm=counts.div(counts.sum(0),axis=1)*1e6; keep=(cpm>=min_cpm).mean(1)>=min_frac
    else:
        keep=logm.mean(1)>=np.nanpercentile(logm.mean(1),25)
    return logm.loc[keep]
def zscore_rows(m): return m.sub(m.mean(1),axis=0).div(m.std(1).replace(0,np.nan),axis=0)
def set_scores(logm,sets=GENESETS,min_genes=2):
    z=zscore_rows(logm); rows={}; used={}
    for k,g in sets.items():
        gg=[x for x in g if x in z.index]
        if len(gg)>=min_genes: rows[k]=z.loc[gg].mean(0); used[k]=len(gg)
    return pd.DataFrame(rows), used
def ols_test(y,meta,group_col='case',covars=()):
    """y: Series indexed by sample. Returns beta (case - control in y units), t, p, standardized d"""
    d=meta.copy(); d['y']=y.loc[d.index].values
    d=d.dropna(subset=['y'])
    f='y ~ '+' + '.join([group_col]+list(covars))
    try: r=smf.ols(f,d).fit()
    except Exception as e: return dict(beta=np.nan,t=np.nan,p=np.nan,d=np.nan,n_case=np.nan,n_ctrl=np.nan)
    b=r.params.get(group_col,np.nan); sd=d['y'].std()
    return dict(beta=b,t=r.tvalues.get(group_col,np.nan),p=r.pvalues.get(group_col,np.nan),d=b/sd if sd>0 else np.nan,
                n_case=int((d[group_col]==1).sum()),n_ctrl=int((d[group_col]==0).sum()))
def test_matrix(mat,meta,covars=(),group_col='case'):
    out=[]
    for g in mat.index:
        r=ols_test(mat.loc[g],meta,group_col,covars); r['feature']=g; out.append(r)
    df=pd.DataFrame(out).set_index('feature'); df['q']=bh(df['p']); return df
def paired_change(mat,meta,subject_col,time_col,t0,t1):
    """returns delta matrix (features x subjects) and per-subject meta"""
    a=meta[meta[time_col]==t0]; b=meta[meta[time_col]==t1]
    common=sorted(set(a[subject_col])&set(b[subject_col]))
    ia=a.reset_index().set_index(subject_col).loc[common,a.index.name or 'index']
    ib=b.reset_index().set_index(subject_col).loc[common,b.index.name or 'index']
    delta=pd.DataFrame(mat[ib.values].values-mat[ia.values].values,index=mat.index,columns=common)
    sm_=a.reset_index().set_index(subject_col).loc[common]
    return delta,sm_

def run_case_control(name,logm,meta,covars=(),outdir='.',extra_label=''):
    """logm: genes x samples (already filtered); meta indexed by sample with 'case' 0/1."""
    meta=meta.loc[[c for c in logm.columns if c in meta.index]]; logm=logm[meta.index]
    tg=[g for g in TARGET_GENES if g in logm.index]
    gres=test_matrix(logm.loc[tg],meta,covars); gres['dataset']=name; gres['contrast']='case_vs_control'+extra_label
    sc,used=set_scores(logm); sres=test_matrix(sc.T,meta,covars); sres['n_genes']=pd.Series(used); sres['dataset']=name; sres['contrast']='case_vs_control'+extra_label
    missing=[g for g in TARGET_GENES if g not in logm.index]
    return gres,sres,sc,missing
def run_paired(name,logm,meta,subject_col,time_col,t0,t1,covars=(),label=None):
    """change t0->t1 per subject; tests (a) case vs control difference in change, (b) within-group mean change"""
    meta=meta.loc[[c for c in logm.columns if c in meta.index]]; logm=logm[meta.index]
    meta=meta.copy(); meta.index.name='sample'
    tg=[g for g in TARGET_GENES if g in logm.index]
    sc,used=set_scores(logm)
    feats=pd.concat([logm.loc[tg],sc.T])
    delta,sm_=paired_change(feats,meta,subject_col,time_col,t0,t1)
    lab=label or f'{t0}->{t1}'
    res=test_matrix(delta,sm_,covars); res['contrast']=f'change_{lab}_case_minus_ctrl'
    rows=[]
    for grp,gv in [('case',1),('control',0)]:
        cols=sm_.index[sm_['case']==gv]
        for f in delta.index:
            x=delta.loc[f,cols].dropna()
            if len(x)>=3:
                t,p=stats.ttest_1samp(x,0); rows.append(dict(feature=f,group=grp,mean_change=x.mean(),d=x.mean()/x.std() if x.std()>0 else np.nan,p=p,n=len(x)))
    within=pd.DataFrame(rows)
    if len(within): within['q']=within.groupby('group')['p'].transform(lambda s: pd.Series(bh(s.values),index=s.index))
    res['dataset']=name; within['dataset']=name; within['contrast']=f'within_{lab}'
    res['is_set']=res.index.isin(sc.columns)
    return res,within

def competitive(logm,meta,covars=(),nperm=300,seed=0,sets=GENESETS):
    """Compare each set's case-control d to random gene sets of the same size (expression-matched by decile)."""
    rng=np.random.default_rng(seed)
    meta=meta.loc[[c for c in logm.columns if c in meta.index]]; lm=logm[meta.index]
    z=zscore_rows(lm).dropna()
    # per-gene d via fast OLS residualisation
    X=np.column_stack([np.ones(len(meta)),meta['case'].values]+[meta[c].values.astype(float) for c in covars])
    B=np.linalg.lstsq(X,z.T.values,rcond=None)[0]; gd=pd.Series(B[1],index=z.index)  # z-units so ~d
    dec=pd.qcut(lm.loc[z.index].mean(1),10,labels=False)
    out={}
    for k,g in sets.items():
        gg=[x for x in g if x in z.index]
        if len(gg)<2: continue
        obs=gd[gg].mean(); null=[]
        pools=[dec[dec==dec[x]].index.values for x in gg]
        for _ in range(nperm): null.append(gd[[rng.choice(p) for p in pools]].mean())
        null=np.array(null); out[k]=dict(mean_gene_d=obs,null_mean=null.mean(),null_sd=null.std(),p_comp=(np.sum(np.abs(null-null.mean())>=abs(obs-null.mean()))+1)/(nperm+1))
    df=pd.DataFrame(out).T; df['z_comp']=(df.mean_gene_d-df.null_mean)/df.null_sd; df['q_comp']=bh(df.p_comp); return df
