# Summarise BCR (PRJEB83659) and TCR (PRJNA940635) repertoires per sample; group tests; VDJdb EBV/CMV matches.
import pandas as pd, numpy as np, glob, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'; rng=np.random.default_rng(1)
def vgrp(v):
    v=v.replace('D','') if v.endswith('D') else v
    return 'IGHV3-30/3-33' if v in('IGHV3-30','IGHV3-30-3','IGHV3-30-5','IGHV3-33') else v
def clon(x,n):
    if len(x)<n: return dict(top10=np.nan,gini_simpson=np.nan,uniq=np.nan)
    s=pd.Series(rng.choice(x,n,replace=False)).value_counts()
    p=s/n; return dict(top10=s.head(10).sum()/n, simpson=(p**2).sum(), uniq=len(s)/n)
def tests(D,feats,groups,label):
    out=[]
    for f in feats:
        d=D.dropna(subset=[f]).copy(); d['z']=(d[f]-d[f].mean())/d[f].std()
        for gname,sub in groups.items():
            dd=d[d.grp.isin(sub[0]+sub[1])].copy(); dd['x']=dd.grp.isin(sub[0]).astype(int)
            if dd.x.nunique()<2 or dd.z.std()==0: continue
            m=smf.ols('z~x'+('+C(covar)' if 'covar' in dd and dd.covar.nunique()>1 else ''),dd).fit()
            out.append(dict(data=label,feature=f,contrast=gname,n1=int(dd.x.sum()),n0=int((dd.x==0).sum()),mean1=dd.loc[dd.x==1,f].mean(),mean0=dd.loc[dd.x==0,f].mean(),d=m.params['x'],p=m.pvalues['x']))
    O=pd.DataFrame(out); O['q']=O.groupby('contrast').p.transform(lambda p:multipletests(p,method='fdr_bh')[1]); return O
# ---------- BCR ----------
A=pd.read_csv('/tmp/work_scm/bcr/sample_attrs.csv').set_index('alias')
B=pd.concat([pd.read_parquet(f) for f in glob.glob('/tmp/work_scm/bcr/parsed/*.parquet')])
RU=pd.read_csv('/tmp/work_scm/bcr/runs_ftp.tsv',sep='\t').sort_values('run_accession').groupby('library_name').head(1)
B=B[B.run.isin(RU.run_accession)]
B['vg']=B.v.map(vgrp); B['clone']=B.vg+'_'+B.cdr3.fillna('')
rows=[]
for s,d in B.groupby('sample'):
    r=dict(sample=s,n=len(d),frac_IgG=(d.iso=='IgG').mean())
    for iso in ['IgM','IgG']:
        di=d[d.iso==iso]; dc=di[di.cdr3.notna()]
        r[f'{iso}_cdr3len']=dc.cdr3.str.len().mean(); r[f'{iso}_shm_median']=di.shm.median()
        for v in ['IGHV3-30/3-33','IGHV4-34','IGHV1-69','IGHV3-23','IGHV4-39','IGHV3-7','IGHV1-2','IGHV5-51']:
            r[f'{iso}_{v}']=(di.vg==v).mean()
        c=clon(dc.clone.values,1500 if iso=='IgM' else 300); r.update({f'{iso}_{k}':v for k,v in c.items()})
        # clones seen in both replicate PCR libraries (A and B)
        la=set(dc[dc.lib.str.endswith('A')].clone); lb=set(dc[dc.lib.str.endswith('B')].clone)
        r[f'{iso}_shared_AB_frac']=len(la&lb)/max(1,min(len(la),len(lb))) if la and lb else np.nan
    rows.append(r)
SB=pd.DataFrame(rows).set_index('sample').join(A[['study_group']]); SB['grp']=SB.study_group
SB.to_csv(f'{R}/bcr_sample_metrics.csv')
print(SB.groupby('grp').n.describe())
feats=[c for c in SB.columns if c not in('n','study_group','grp')]
G={'ME_all_vs_HC':(['ME/CFSmm','ME/CFSsa'],['Healthy_control']),'MEsevere_vs_HC':(['ME/CFSsa'],['Healthy_control']),'MEsevere_vs_MEmild':(['ME/CFSsa'],['ME/CFSmm']),'ME_all_vs_MS':(['ME/CFSmm','ME/CFSsa'],['MS_control'])}
OB=tests(SB,feats,G,'BCR'); OB.to_csv(f'{R}/bcr_group_tests.csv',index=False)
pd.set_option('display.width',220); print(OB.sort_values('p').head(25).round(4).to_string())
print(OB[OB.feature.str.contains('3-30/3-33|4-34')].round(4).to_string())
# ---------- TCR ----------
TA=pd.read_csv('/tmp/work_scm/tcr/sample_attrs.csv').set_index('sample_accession')
T=pd.concat([pd.read_parquet(f) for f in glob.glob('/tmp/work_scm/tcr/parsed/*.parquet')]); T=T[T['prod']==True]
vd=pd.read_csv('/tmp/work_scm/tcr/vdjdb/vdjdb-2026-06-03/vdjdb.slim.txt',sep='\t',low_memory=False)
vd=vd[(vd.species=='HomoSapiens')&(vd.gene=='TRB')]
spec={}
for sp in ['EBV','CMV','InfluenzaA','SARS-CoV-2','HomoSapiens']:
    spec[sp]=set(vd.loc[vd['antigen.species']==sp,'cdr3'])
rows=[]
for s,d in T.groupby('sample'):
    x=d.cdr3.values; r=dict(sample=s,n=len(x))
    c=clon(x,800); r.update(c)
    sub=rng.choice(x,800,replace=False) if len(x)>=800 else x
    su=pd.Series(sub)
    for sp,st in spec.items():
        m=su.isin(st); r[f'{sp}_reads_frac']=m.mean(); r[f'{sp}_uniq_frac']=su[m].nunique()/su.nunique()
        top=su.value_counts().head(20).index; r[f'{sp}_in_top20']=np.isin(top,list(st)).sum()
    r['EBV_minus_flu_uniq']=r['EBV_uniq_frac']-r['InfluenzaA_uniq_frac']
    rows.append(r)
ST=pd.DataFrame(rows).set_index('sample').join(TA[['isolate','tissue','age']]); ST['grp']=ST.isolate; ST['cell']=np.where(ST.tissue.str.contains('CD8'),'CD8','CD4'); ST['covar']=ST.age
ST.to_csv(f'{R}/tcr_sample_metrics.csv'); print(ST.groupby(['cell','grp']).n.describe())
G2={'ME_all_vs_HC':(['ME/CFSmm','ME/CFSsa'],['Healthy control']),'MEsevere_vs_HC':(['ME/CFSsa'],['Healthy control']),'MEsevere_vs_MEmild':(['ME/CFSsa'],['ME/CFSmm']),'ME_all_vs_MS':(['ME/CFSmm','ME/CFSsa'],['MS control'])}
feats=['top10','simpson','uniq']+[c for c in ST.columns if c.endswith(('_reads_frac','_uniq_frac','_in_top20'))]+['EBV_minus_flu_uniq']
OT=pd.concat([tests(ST[ST.cell==c],feats,G2,'TCR_'+c) for c in ['CD4','CD8']]); OT.to_csv(f'{R}/tcr_group_tests.csv',index=False)
print(OT.sort_values('p').head(25).round(4).to_string())
