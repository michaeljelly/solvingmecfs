# Build compact DecodeME outcome tables (rsID-keyed) + derive MVP-only estimates from the DME-1+MVP METAL meta.
import pyarrow.parquet as pq, pyarrow as pa, pyarrow.csv as pc, pandas as pd, numpy as np
W='/tmp/work_decodeme/pq/'; O='/tmp/work_gmr/'
m=pc.read_csv(W+'id_rs_map.tsv.gz',parse_options=pc.ParseOptions(delimiter='\t'))
for g in ['gwas_1','gwas_1_infectious_onset']:
    t=pq.read_table(W+f'{g}.parquet',columns=['ID','CHROM','GENPOS','ALLELE0','ALLELE1','A1FREQ','BETA','SE','LOG10P','N_CASES','N_CONTROLS'])
    t=t.set_column(0,'ID',t.column('ID').cast(pa.string()))
    j=t.join(m,'ID',join_type='inner')
    d=j.to_pandas()
    d=d.rename(columns={'CHROM':'chr38','GENPOS':'pos38','ALLELE1':'ea','ALLELE0':'oa','A1FREQ':'eaf','BETA':'beta','SE':'se','LOG10P':'mlogp'})
    d=d[['rs','chr','pos19','pos38','ea','oa','eaf','beta','se','mlogp','N_CASES','N_CONTROLS']]
    d=d.drop_duplicates('rs',keep=False)
    d.to_parquet(O+f'dme_{g}.parquet',index=False); print(g,len(d),flush=True)
    del d,j,t
# MVP-only by subtraction. The METAL file used the sample-size (Z) scheme; BETA/SE there are back-derived from Z.
# FRQ and Z refer to the allele in column A2 (checked against DecodeME eaf/sign at all leads), so effect allele = A2.
# Z_meta = (sqrt(Nd) Zd + sqrt(Nm) Zm)/sqrt(Nd+Nm); N is effective N (DME 58,792; MVP 15,427 = 3,891 cases/439,202 controls).
dm=pd.read_parquet(O+'dme_gwas_1.parquet',columns=['rs','ea','oa','beta','se'])
mt=pd.read_csv('/tmp/work_decodeme/meta/GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz',sep='\t',usecols=['SNP','CHR','BP','A1','A2','FRQ','N','Z','DIRECTION'],dtype={'CHR':str})
mt=mt[mt.DIRECTION.str.len()==2]; mt=mt[mt.DIRECTION.str[1]!='?']
mt['ea']=mt.A2.str.upper(); mt['oa']=mt.A1.str.upper()
x=mt.merge(dm,left_on='SNP',right_on='rs',how='left',suffixes=('','_d'))
has=x.beta.notna()&(x.DIRECTION.str[0]!='?')
same=(x.ea==x.ea_d)&(x.oa==x.oa_d); flip=(x.ea==x.oa_d)&(x.oa==x.ea_d)
zd=np.where(same,x.beta/x.se,np.where(flip,-x.beta/x.se,np.nan))
Nd=58792.0
Nm=np.where(has,x.N-Nd,x.N)
zm=np.where(has,(x.Z*np.sqrt(x.N)-np.sqrt(Nd)*zd)/np.sqrt(np.clip(Nm,1,None)),x.Z)
ok=np.isfinite(zm)&(Nm>10000)
p=x.FRQ.clip(0.001,0.999)
se=1/np.sqrt(Nm/4*2*p*(1-p))
r=pd.DataFrame({'rs':x.SNP,'chr':x.CHR,'pos38':x.BP,'ea':x.ea,'oa':x.oa,'eaf':x.FRQ,'z':zm,'beta':zm*se,'se':se,'Neff':Nm})[ok]
r.to_parquet(O+'mvp_only_derived.parquet',index=False); print('mvp',len(r))
chk={'rs4810894':1.2408,'rs7165327':0.8702,'rs2503773':-1.7935,'rs12579722':3.14,'rs1982935':2.835,'rs12071663':-0.404,'rs9358913':-0.154,'rs34626694':-0.645}
for s_,z in chk.items():
    q=r[r.rs==s_]; print(s_,z, q.z.round(3).tolist(), q.ea.tolist())
