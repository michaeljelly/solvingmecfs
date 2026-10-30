# Reverse MR: genetic liability to ME/CFS (DecodeME GWAS-1 genome-wide significant SNPs) -> metabolites / infection / blood traits.
# Instruments are few (<=8) so power is low; results are reported as such.
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/genetics_mr/scripts'); import mrlib, tabixlib
G='/tmp/work_gmr/'; R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
D=pd.read_parquet(G+'dme_gwas_1.parquet'); D['p']=10**-D.mlogp
S=D[D.p<1e-6]
mrlib.BFILE='/tmp/work_decodeme/magma/g1000_eur'
keep=mrlib.clump(S[['rs','p']],p1=5e-8,r2=0.001,kb=10000)
# indel leads absent from 1000G: add best SNV proxy per locus not already represented (lead SNVs at 6q16 and 12q24 are p>5e-8 in GWAS-1 -> not added)
E=D[D.rs.isin(keep)][['rs','chr','pos38','ea','oa','eaf','beta','se','p']]
print('ME/CFS instruments',len(E)); print(E.to_string())
E.to_csv(R+'reverse_mr_instruments.csv',index=False)
targets={'IPA (Chen)':'GCST90199686','Indolelactate (Chen)':'GCST90199650','Kynurenine (Chen)':'GCST90199636','Tryptophan (Chen)':'GCST90200441','3-hydroxybutyrate (Chen)':'GCST90200367',
 'Carnitine (Chen)':'GCST90199621','Glycocholate (Chen)':'GCST90200387','EBV DNA load (Kamitaki)':'GCST90809803','Reticulocyte count (Vuckovic)':'GCST90002405',
 'Immature reticulocyte fraction (Vuckovic)':'GCST90002387','Haemoglobin (Chen MH)':'GCST90002310','EBV EBNA-1 IgG':'GCST90006899','EBV ZEBRA IgG':'GCST90006901','25(OH)D (Revez)':'GCST90000618'}
res=[]
def run(lab,o):
    x=mrlib.harmonise(E[['rs','ea','oa','eaf','beta','se']],o)
    if len(x)==0: return
    r=mrlib.mr_all(x,nboot=300); r['outcome']=lab; res.append(r); print(lab,{k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()},flush=True)
for lab,acc in targets.items():
    try:
        u=tabixlib.acc_url(acc); rows=[]
        for _,e in E.iterrows():
            q=tabixlib.std(tabixlib.region(u,e.chr,e.pos38-1,e.pos38+1)); rows.append(q[q.rs==e.rs])
        o=pd.concat(rows); run(lab,o[['rs','ea','oa','eaf','beta','se']])
    except Exception as ex:
        # no tabix index: fall back to the locally streamed HapMap3 / p<5e-6 extracts
        loc={'GCST90002405':'RetCount_Vuckovic2020','GCST90002387':'IRF_Vuckovic2020','GCST90002310':'Hb_Chen2020','GCST90000618':'VitD_Revez2020'}.get(acc)
        o=None
        if loc:
            import os
            fs=[f for f in [G+f'exp/{loc}.hm3.tsv.gz',G+f'exp/{loc}.tsv'] if os.path.exists(f)]
            if fs:
                o=pd.concat([pd.read_csv(f,sep='\t',dtype=str) for f in fs]).drop_duplicates('rs')
                for c in ['eaf','beta','se']: o[c]=pd.to_numeric(o[c],errors='coerce')
                o['ea']=o.ea.str.upper(); o['oa']=o.oa.str.upper(); o=o[o.rs.isin(set(E.rs))]
        if o is not None and len(o): run(lab+' [local extract, %d/5 SNPs]'%len(o),o[['rs','ea','oa','eaf','beta','se']])
        else: print(lab,'fail',ex)
# local full files
for lab,path,cols in [('Multisite chronic pain','/tmp/work_decodeme/other/chronic_pain-bgen.stats.gz',dict(SNP='rs',ALLELE1='ea',ALLELE0='oa',A1FREQ='eaf',BETA='beta',SE='se')),
                      ('Major depression','/tmp/work_decodeme/other/GCST90726344.h.tsv.gz',dict(rsid='rs',effect_allele='ea',other_allele='oa',effect_allele_frequency='eaf',beta='beta',standard_error='se')),
                      ('Infectious mononucleosis FinnGen','/tmp/work_gmr/dl/fg_AB1_EBV.tsv.gz',dict(rsids='rs',alt='ea',ref='oa',af_alt='eaf',beta='beta',sebeta='se'))]:
    sep='\t'
    parts=[c.rename(columns=cols)[list(cols.values())].pipe(lambda z: z[z.rs.isin(set(E.rs))]) for c in pd.read_csv(path,sep=r'\s+' if 'chronic' in path else sep,usecols=list(cols.keys()),chunksize=2_000_000)]
    run(lab,pd.concat(parts))
pd.DataFrame(res).to_csv(R+'reverse_mr_results.csv',index=False)
