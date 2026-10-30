# MHC region: DecodeME signal by subregion, DQ2.5 (DQA1*05:01-DQB1*02:01) tag, and overlap with EBV/CMV/IM MHC associations.
import pandas as pd, numpy as np, glob, os, sys
from scipy import stats
sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/genetics_mr/scripts'); import mrlib
G='/tmp/work_gmr/'; R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
d1=pd.read_parquet(G+'mhc_dme_gwas_1.parquet'); di=pd.read_parquet(G+'mhc_dme_gwas_1_infectious_onset.parquet'); mv=pd.read_parquet(G+'mhc_mvp.parquet')
ni=pd.read_parquet(G+'mhc_dme_noninf_raw.parquet'); ni=ni.merge(d1[['rs','pos38','ea','oa']],left_on=['GENPOS','ALLELE1','ALLELE0'],right_on=['pos38','ea','oa'])
ni=ni.rename(columns={'BETA':'beta','SE':'se','A1FREQ':'eaf'})[['rs','ea','oa','eaf','beta','se']]
OUTS={'GWAS1':d1[['rs','ea','oa','eaf','beta','se']],'infection':di[['rs','ea','oa','eaf','beta','se']],'non_infection':ni,'MVP':mv[['rs','ea','oa','eaf','beta','se']]}
# 1. DQ2.5 tags and classical-allele tag SNPs
tags={'rs2040410':'DQ2.5 (DR3-DQA1*05:01-DQB1*02:01) proxy r2=0.98 with rs2187668','rs2854275':'DQ2.5 proxy r2=0.97','rs3135388':'DRB1*15:01 tag (MS/EBNA-1)','rs9271366':'DRB1*15:01 tag',
      'rs2395029':'HLA-B*57:01 tag','rs6457110':'HLA-A*01/IM tag region','rs9275596':'DQB1 region','rs9358913':'6p22.2 DecodeME lead (HIST1/BTN)'}
rows=[]
for rs,lab in tags.items():
    r=dict(rs=rs,tag=lab)
    for k,o in OUTS.items():
        x=o[o.rs==rs]
        if len(x): r[k+'_ea']=x.ea.iloc[0]; r[k+'_OR']=np.exp(x.beta.iloc[0]); r[k+'_p']=2*stats.norm.sf(abs(x.beta.iloc[0]/x.se.iloc[0]))
    rows.append(r)
T=pd.DataFrame(rows)
# heterogeneity infection vs non-infection onset at DQ2.5 tag (shared controls -> test is conservative-ish approximation)
for rs in ['rs2040410','rs2854275']:
    a=di[di.rs==rs].iloc[0]; b=OUTS['non_infection'][OUTS['non_infection'].rs==rs].iloc[0]
    zb=(a.beta-b.beta)/np.sqrt(a.se**2+b.se**2); print(rs,'inf vs noninf beta',a.beta,b.beta,'het p',2*stats.norm.sf(abs(zb)))
T.to_csv(R+'hla_tag_lookup.csv',index=False); print(T.round(4).to_string())
# 2. per EBV/CMV/IM trait: MHC lead SNPs (clumped r2<0.1, p<5e-8 or 1e-5 if none), DecodeME lookup + MHC-restricted IVW
res=[]; corr=[]
for f in sorted(glob.glob(G+'exp/*.mhc.tsv')):
    eid=os.path.basename(f).replace('.mhc.tsv','')
    e=pd.read_csv(f,sep='\t',dtype=str)
    for c in ['pos','eaf','beta','se','p']: e[c]=pd.to_numeric(e[c],errors='coerce')
    e=e[e.rs.astype(str).str.startswith('rs')].dropna(subset=['beta','se','p']); e['ea']=e.ea.str.upper(); e['oa']=e.oa.str.upper()
    if len(e)==0: continue
    thr=5e-8 if (e.p<5e-8).sum()>0 else 1e-5
    lead=e.sort_values('p').iloc[0]
    for r2c in [0.1,0.01]:
        cl=mrlib.clump(e[['rs','p']],p1=thr,r2=r2c,kb=3000)
        inst=e[e.rs.isin(cl)]
        for k,o in OUTS.items():
            x=mrlib.harmonise(inst[['rs','ea','oa','eaf','beta','se']],o)
            if len(x)==0: continue
            agree=(np.sign(x.beta_e)==np.sign(x.beta_o)).sum()
            b,se,Q,Qp=mrlib.ivw(x.beta_e.values,x.beta_o.values,x.se_o.values) if len(x)>1 else (x.beta_o.iloc[0]/x.beta_e.iloc[0],abs(x.se_o.iloc[0]/x.beta_e.iloc[0]),np.nan,np.nan)
            # all-SNP z correlation across MHC (descriptive; SNPs are in strong LD so no valid p)
            y=mrlib.harmonise(e[['rs','ea','oa','eaf','beta','se']],o); z1_=(y.beta_e/y.se_e).values; z2_=(y.beta_o/y.se_o).values; okz=np.isfinite(z1_)&np.isfinite(z2_); r=np.corrcoef(z1_[okz],z2_[okz])[0,1]
            lk=o[o.rs==lead.rs]
            res.append(dict(trait=eid,group=('EBV lytic' if eid.split('_')[0] in ('ZEBRA','EAD','VCAp18') else 'EBV latent/other EBV' if eid.split('_')[0] in ('EBNA1','EBVsero','EBVload','IM','EBVreadpos') else 'CMV' if eid.startswith('CMV') else 'control pathogen'),clump_r2=r2c,outcome=k,thr=thr,trait_lead=lead.rs,trait_lead_pos=lead.pos,trait_lead_p=lead.p,outcome_p_at_trait_lead=(2*stats.norm.sf(abs(lk.beta.iloc[0]/lk.se.iloc[0])) if len(lk) else np.nan),
                            n_inst=len(x),n_same_sign=int(agree),ivw_b=b,ivw_se=se,ivw_p=2*stats.norm.sf(abs(b/se)),Q_p=Qp,z_corr_allMHC=r,n_mhc_snps=len(y)))
H=pd.DataFrame(res); H.to_csv(R+'hla_ebv_cmv_overlap.csv',index=False); print(H.round(4).to_string())
