# Epigenetic clocks, DNAm surrogates (GrimAge components, smoking, BMI, alcohol), cell composition (Salas 12-cell NNLS),
# per dataset; group comparison adjusted for sex (+ chip batch, + age where available).
import pandas as pd, numpy as np, re, os, sys
from scipy.optimize import nnls
from scipy.stats import rankdata
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
W='/tmp/work_scm/meth'; R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
BL='/usr/local/lib/python3.11/dist-packages/biolearn/data'
ann450=pd.read_parquet(f'{W}/ann450.parquet').set_index('IlmnID'); annE=pd.read_parquet(f'{W}/annEPIC.parquet').set_index('IlmnID')
ref_avg=pd.read_csv(f'{BL}/biolearn_averages_450k.csv',index_col=0).iloc[:,0]
def kv(M):
    out=pd.DataFrame(index=M.index)
    for c in M.columns:
        if c.startswith('characteristics'):
            for i,v in M[c].items():
                if ':' in str(v):
                    k,val=str(v).split(':',1); out.loc[i,k.strip()]=val.strip()
    out['gsm']=M['geo_accession']; out['title']=M['title']; return out
def load(g):
    B=pd.read_parquet(f'{W}/{g}_beta.parquet'); M=kv(pd.read_csv(f'{W}/{g}_meta_raw.csv'))
    if g=='GSE327255':
        M.index=B.columns  # same order as series (54 EM then 21 controls); checked by title numbers below
        M['group']=np.where(M['subject group'].str.upper().str.startswith('EM'),'ME','HC')
        M['age']=M['age at time of test'].astype(float); M['sex']=M['Sex'].map({'F':'F','M':'M'}); M['bmi']=M['bmi'].astype(float)
    else:
        M.index=M['gsm']; B=B[M.index]
        col=[c for c in M.columns if c in('sample group','disease state')][0]
        M['group']=np.where(M[col].str.contains('CFS|ME/CFS|Myalgic',regex=True),'ME','HC')
        if g=='GSE156792': M['chip']=M['title'].str.extract(r'(\d{12})_')[0]
    return B,M
def infer_sex(B,ann):
    x=ann.index[(ann.CHR=='X')&(ann.Relation_to_UCSC_CpG_Island=='Island')].intersection(B.index)
    y=ann.index[ann.CHR=='Y'].intersection(B.index)
    xm=B.loc[x].mean(); ym=B.loc[y].mean() if len(y) else pd.Series(np.nan,index=B.columns)
    return xm,ym
def fill(B,cpgs):
    X=B.reindex(cpgs)
    miss=X.isna().all(1).sum()
    X=X.apply(lambda r:r.fillna(r.mean()),axis=1)
    X=X.apply(lambda r:r.fillna(ref_avg.get(r.name,0.5)),axis=1)
    return X,miss
def lin(B,file,transform=None):
    c=pd.read_csv(f'{BL}/{file}',index_col=0).iloc[:,0]; ic=c.get('intercept',0.0); c=c.drop('intercept',errors='ignore')
    X,miss=fill(B,c.index); s=(X.mul(c,axis=0)).sum()+ic
    return (transform(s) if transform else s), miss, len(c)
def anti(x,a=20): x=np.asarray(x,float); return np.where(x<0,(1+a)*np.exp(x)-1,(1+a)*x+a)
def dunedin(B):
    gold=pd.read_csv(f'{BL}/DunedinPACE_Gold_Means.csv',index_col=0)['mean']
    c=pd.read_csv(f'{BL}/DunedinPACE.csv',index_col=0).iloc[:,0]; ic=c.get('intercept'); c=c.drop('intercept')
    X=B.reindex(gold.index); miss=X.isna().all(1).sum()
    X=X.apply(lambda r:r.fillna(gold[r.name]) if r.isna().all() else r.fillna(r.mean()),axis=1)
    st=np.sort(gold.values); out=X.copy()
    for col in X.columns:
        r=rankdata(X[col].values); f=np.floor(r).astype(int); hi=(r-f)>0.4
        v=np.where(hi,0.5*(st[np.clip(f-1,0,len(st)-1)]+st[np.clip(f,0,len(st)-1)]),st[f-1]); out[col]=v
    return (out.reindex(c.index).mul(c,axis=0)).sum()+ic, miss, len(c)
def grim_components(B):
    G=pd.read_csv(f'{BL}/GrimAgeV1.csv'); res={}
    for name,g in G.groupby('Y.pred'):
        if name in('COX','transform'): continue
        g=g.set_index('var')['beta']; cp=g.drop([k for k in ['Intercept','Age','Female'] if k in g.index])
        X,miss=fill(B,cp.index); res[name]=(X.mul(cp,axis=0)).sum()+g.get('Intercept',0)  # CpG part + intercept (age/sex terms handled as covariates)
    return pd.DataFrame(res)
def deconv(B):
    ref=pd.read_csv(f'{BL}/twelve_cell_deconv.csv',index_col=0)
    cp=ref.index.intersection(B.index); X=B.loc[cp]; X=X.apply(lambda r:r.fillna(r.mean()),axis=1)
    out={}
    for s in X.columns:
        w,_=nnls(ref.loc[cp].values,X[s].values); out[s]=w/ w.sum() if w.sum()>0 else w
    return pd.DataFrame(out,index=ref.columns).T, len(cp), len(ref)
rows=[]; allscores={}
for g in ['GSE156792','GSE93266','GSE111183','GSE327255']:
    B,M=load(g); ann=annE if g in('GSE111183','GSE327255') else ann450
    xm,ym=infer_sex(B,ann); M['xisland_mean']=xm; M['y_mean']=ym
    M['sex_inf']=np.where(xm>0.35,'F','M')
    if g=='GSE327255': M['sex_inf']=M['sex']
    S=pd.DataFrame(index=B.columns); cov={}
    for nm,f,tr in [('Horvath',"Horvath1.csv",lambda s:anti(s+0.696)),('HorvathSkinBlood','Horvath2.csv',lambda s:anti(s-0.447119319)),
                    ('Hannum','Hannum.csv',None),('PhenoAge','PhenoAge.csv',None),('DNAmTL','DNAmTL.csv',None),('Zhang10_mortality','Zhang_10.csv',None),
                    ('Smoking_McCartney','Smoking.csv',None),('Alcohol_McCartney','Alcohol.csv',None),('BMI_McCartney','BMI_McCartney.csv',None),('EpiTOC2_like_Lin','Lin.csv',None)]:
        v,miss,n=lin(B,f,tr); S[nm]=np.asarray(v); cov[nm]=f'{n-miss}/{n}'
    v,miss,n=dunedin(B); S['DunedinPACE']=np.asarray(v); cov['DunedinPACE']=f'{n}(gold missing {miss})'
    Gc=grim_components(B); S=S.join(Gc)
    Cf,nc,nr=deconv(B); S=S.join(Cf.add_prefix('frac_')); cov['deconv']=f'{nc}/{nr}'
    rr=pd.read_csv(f'{BL}/450K_reinius_12_reference.csv',index_col=0); cp=rr.index.intersection(B.index); X=B.loc[cp].apply(lambda r:r.fillna(r.mean()),axis=1)
    R6=pd.DataFrame({s_:nnls(rr.loc[cp].values,X[s_].values)[0] for s_ in X.columns},index=rr.columns).T; R6=R6.div(R6.sum(1),axis=0)
    S=S.join(R6.add_prefix('rein_')); cov['reinius']=f'{len(cp)}/{len(rr)}'
    S['cg05575921_AHRR']=B.reindex(['cg05575921']).T.iloc[:,0] if 'cg05575921' in B.index else np.nan
    D=M.join(S)
    D.to_csv(f'{R}/meth_{g}_sample_scores.csv')
    allscores[g]=D
    print(g,B.shape,dict(D.group.value_counts()),'sex_inf',dict(D.sex_inf.value_counts()),'coverage',cov,flush=True)
    covars=['C(sex_inf)'] if D.sex_inf.nunique()>1 else []
    if g=='GSE156792': covars.append('C(chip)')
    if g=='GSE327255': covars.append('age')
    for col in S.columns:
        d=D.dropna(subset=[col]).copy(); d['y']=(d[col]-d[col].mean())/d[col].std()
        if d['y'].std()==0 or d['y'].isna().all(): continue
        f='y ~ C(group, Treatment("HC"))'+(' + '+' + '.join(covars) if covars else '')
        try: m=smf.ols(f,data=d).fit()
        except Exception as e: print(col,e); continue
        k=[x for x in m.params.index if 'group' in x][0]
        rows.append(dict(dataset=g,score=col,n_ME=(d.group=='ME').sum(),n_HC=(d.group=='HC').sum(),mean_ME=d.loc[d.group=='ME',col].mean(),mean_HC=d.loc[d.group=='HC',col].mean(),
                         d_adj=m.params[k],p=m.pvalues[k],covars='+'.join(covars),coverage=cov.get(col,'')))
    if g=='GSE327255':  # age acceleration + illness duration
        for col in ['Horvath','HorvathSkinBlood','Hannum','PhenoAge','DNAmTL','DunedinPACE']:
            d=D.copy(); d['acc']=smf.ols(f'{col} ~ age',data=d[d.group=='HC']).fit().predict(d); d['acc']=d[col]-d['acc']
            print(' ',col,'r(age)=%.2f'%np.corrcoef(d[col],d.age)[0,1],'acc ME %.2f HC %.2f'%(d.loc[d.group=='ME','acc'].mean(),d.loc[d.group=='HC','acc'].mean()))
res=pd.DataFrame(rows); res['q_within_dataset']=res.groupby('dataset').p.transform(lambda p:multipletests(p,method='fdr_bh')[1])
res.to_csv(f'{R}/meth_clocks_composition_group_tests.csv',index=False)
print(res.sort_values('p').head(40).to_string())
