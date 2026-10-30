# Load DecodeME REGENIE sumstats, keep QC-passed variants (per README), save parquet
import pandas as pd, numpy as np, sys
W='/tmp/work_decodeme/'
qc=set(pd.read_csv(W+'ss/gwas_qced.var.gz',header=0)['ID'])
print('QC variants',len(qc))
for g in ['gwas_1','gwas_2','gwas_1_female','gwas_1_male','gwas_1_infectious_onset','gwas_1_non_infectious_onset']:
    parts=[]
    for ch in pd.read_csv(W+f'ss/{g}.regenie.gz',sep=' ',chunksize=2_000_000,
        usecols=['CHROM','GENPOS','ID','ALLELE0','ALLELE1','A1FREQ','A1FREQ_CASES','N','N_CASES','N_CONTROLS','BETA','SE','LOG10P','EXTRA'],
        dtype={'CHROM':str}):
        parts.append(ch[ch.ID.isin(qc)])
    d=pd.concat(parts); d['CHROM']=d.CHROM.astype(int)
    print(g,'rows QC',len(d),'Ncases',d.N_CASES.max(),'Nctrl',d.N_CONTROLS.max(),'firth fail',d.EXTRA.notna().sum(),
          'GW sig',(d.LOG10P>-np.log10(5e-8)).sum(),flush=True)
    d.drop(columns='EXTRA').to_parquet(W+f'pq/{g}.parquet',index=False)
