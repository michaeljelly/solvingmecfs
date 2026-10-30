# LD (1000G EUR, MAGMA reference) between each locus' GWAS-1 top SNP and the top eQTLGen SNP for candidate genes
import pandas as pd, numpy as np
from bed_reader import open_bed
W='/tmp/work_decodeme/'
R=pd.read_csv(W+'results/coloc_eqtlgen_blood.csv'); E=pd.read_csv(W+'results/eqtlgen_region_extract.tsv.gz',sep='\t')
bim=pd.read_csv(W+'magma/g1000_eur.bim',sep='\t',header=None,usecols=[1],names=['rs']).rs
pos={r:i for i,r in enumerate(bim.values)}
bed=open_bed(W+'magma/g1000_eur.bed',count_A1=True,iid_count=503,sid_count=len(bim))
out=[]
for _,r in R[R.gene.isin(['CSE1L','ARFGEF2','ZNFX1','STAU1','BTN2A2','BTN3A2','TRIM38','RABGAP1L','MRPS14','SUDS3','TAOK3','PEBP1','OLFM4','FBXL4','ZNF184','HIST1H2BF','CACYBP'])].iterrows():
    e=E[(E.locus==r.locus)&(E.gene==r.gene)].sort_values('P').iloc[0]
    a,b=r.gwas_top_snp,e.rs
    if a in pos and b in pos:
        X=bed.read(index=np.s_[:,[pos[a],pos[b]]],dtype='float32'); X=X[~np.isnan(X).any(1)]
        r2=np.corrcoef(X.T)[0,1]**2
    else: r2=np.nan
    out.append(dict(locus=r.locus,gene=r.gene,gwas_top=a,eqtl_top=b,eqtl_top_p=e.P,r2_1000G_EUR=round(r2,3),eqtl_z_at_gwas_top=r.eqtl_z_at_gwas_top,PP_H4=r.PP_H4,PP_H3=r.PP_H3,risk_allele_raises_expression=r.risk_allele_raises_expression))
O=pd.DataFrame(out).sort_values(['locus','r2_1000G_EUR'],ascending=[True,False]); O.to_csv(W+'results/coloc_ld_check.csv',index=False); print(O.to_string())
