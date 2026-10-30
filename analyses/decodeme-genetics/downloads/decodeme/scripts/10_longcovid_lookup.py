# Long COVID (HGI, Lammi 2025, GCST90454541, GRCh38) vs DecodeME GWAS-1: lead lookups both ways and sign concordance of independent suggestive DecodeME SNPs
import pandas as pd, numpy as np
from scipy.stats import binomtest, norm
W='/tmp/work_decodeme/'
g=pd.read_parquet(W+'pq/gwas_1.parquet')
# independent suggestive SNPs (p<1e-5), greedy distance clumping 500kb
s=g[g.LOG10P>5].sort_values('LOG10P',ascending=False); keep=[]
for i,r in s.iterrows():
    if not any((k[0]==r.CHROM) and abs(k[1]-r.GENPOS)<500_000 for k in keep): keep.append((r.CHROM,r.GENPOS,r.ID))
ids=[k[2] for k in keep]; sub=g[g.ID.isin(ids)].copy()
key=set(zip(sub.CHROM,sub.GENPOS))|{(6,41515652)}
parts=[]
for c in pd.read_csv(W+'lc/GCST90454541.h.tsv.gz',sep='\t',chunksize=2_000_000,usecols=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','p_value','rsid'],dtype={'chromosome':str}):
    c=c[c.chromosome.str.isdigit()]; c['chromosome']=c.chromosome.astype(int)
    parts.append(c[[ (a,b) in key for a,b in zip(c.chromosome,c.base_pair_location)]])
lc=pd.concat(parts)
x=sub.merge(lc,left_on=['CHROM','GENPOS'],right_on=['chromosome','base_pair_location'])
same=(x.effect_allele==x.ALLELE1)&(x.other_allele==x.ALLELE0); flip=(x.effect_allele==x.ALLELE0)&(x.other_allele==x.ALLELE1)
x=x[same|flip].copy(); x['lc_beta_aligned']=np.where(same[same|flip],x.beta,-x.beta)
x['concordant']=np.sign(x.lc_beta_aligned)==np.sign(x.BETA)
x=x.sort_values('LOG10P',ascending=False)
x[['ID','rsid','LOG10P','BETA','lc_beta_aligned','standard_error','p_value','concordant']].to_csv(W+'results/longcovid_lookup_decodeme_suggestive.csv',index=False)
bt=binomtest(int(x.concordant.sum()),len(x),0.5,alternative='greater')
print('independent DecodeME p<1e-5 SNPs found in LC:',len(x),'concordant',int(x.concordant.sum()),'binomial p',bt.pvalue)
gw=x[x.LOG10P>7.3]; print(gw[['ID','rsid','BETA','lc_beta_aligned','p_value']].to_string())
# inverse-variance weighted combined z of LC effects in direction of DecodeME
zz=np.sign(x.BETA)*x.lc_beta_aligned/x.standard_error; print('mean aligned LC z',zz.mean().round(3),'sum z/sqrt(n) p',norm.sf(zz.sum()/np.sqrt(len(zz))))
# FOXP4 lead in DecodeME
f=g[(g.CHROM==6)&(g.GENPOS==41515652)]; print('FOXP4 rs9367106 in DecodeME GWAS-1:'); print(f[['ID','A1FREQ','BETA','SE','LOG10P']])
f2=g[(g.CHROM==6)&g.GENPOS.between(41_400_000,41_650_000)].sort_values('LOG10P',ascending=False).head(3); print('best DecodeME SNP in FOXP4 region',f2[['ID','BETA','LOG10P']])
for gg in ['gwas_1_infectious_onset']:
    h=pd.read_parquet(W+f'pq/{gg}.parquet'); print(gg,h[(h.CHROM==6)&(h.GENPOS==41515652)][['ID','BETA','SE','LOG10P']])
