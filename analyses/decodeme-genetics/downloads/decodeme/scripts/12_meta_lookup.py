# DecodeME GWAS-1 + Million Veteran Program meta-analysis (Maccallini 2026, Zenodo 20204356): GW-significant loci and DecodeME lead lookups
import pandas as pd, numpy as np
W='/tmp/work_decodeme/'
d=pd.read_csv(W+'meta/GWAS_METAL_DME_1_MVP_GRCh38.tsv.gz',sep='\t',usecols=['SNP','CHR','BP','A1','A2','FRQ','N','Z','P','DIRECTION','BETA','SE'],dtype={'CHR':str})
print('rows',len(d),'N range',d.N.min(),d.N.max())
s=d[d.P<5e-8].sort_values('P'); keep=[]
for _,r in s.iterrows():
    if not any(k[0]==r.CHR and abs(k[1]-r.BP)<500_000 for k in keep): keep.append((r.CHR,r.BP,r.SNP,r.P,r.DIRECTION,r.A1,r.BETA))
K=pd.DataFrame(keep,columns=['CHR','BP','SNP','P','DIRECTION','A1','BETA']); print(K.to_string())
L=pd.read_csv(W+'results/lead_loci_all_gwas.csv')
out=[]
for _,r in L.drop_duplicates('locus').iterrows():
    x=d[(d.CHR==str(r.chr))&(abs(d.BP-r.pos)<250_000)].sort_values('P').head(1)
    y=d[(d.CHR==str(r.chr))&(d.BP==r.pos)]
    out.append(dict(locus=r.locus,decodeme_lead=r.lead,nearest=r.nearest_pc_genes.split('(')[0],meta_p_at_lead=y.P.min() if len(y) else np.nan,meta_dir_at_lead=y.DIRECTION.iloc[0] if len(y) else None,best_meta_p_250kb=x.P.iloc[0],best_meta_snp=x.SNP.iloc[0],best_dir=x.DIRECTION.iloc[0]))
O=pd.DataFrame(out); print(O.to_string()); O.to_csv(W+'results/meta_mvp_lookup.csv',index=False); K.to_csv(W+'results/meta_mvp_gw_loci.csv',index=False)
