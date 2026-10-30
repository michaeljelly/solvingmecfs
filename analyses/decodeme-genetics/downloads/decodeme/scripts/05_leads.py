# Define genome-wide significant loci (distance-merged, +/-500kb) for each DecodeME GWAS; nearest genes; cross-GWAS lookup
import pandas as pd, numpy as np, os
W='/tmp/work_decodeme/'; OUT=W+'results/'; os.makedirs(OUT,exist_ok=True)
loc=pd.read_csv(W+'magma/NCBI38.gene.loc',sep='\t',header=None,names=['id','chr','s','e','str','sym'])
loc=loc[loc.chr.astype(str).str.isdigit()]; loc['chr']=loc.chr.astype(int)
gi=pd.read_csv(W+'ref/Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','type_of_gene'])
pc=set(gi.GeneID[gi.type_of_gene=='protein-coding']); loc=loc[loc.id.isin(pc)]
def near(c,p,k=3):
    g=loc[loc.chr==c].copy(); g['d']=np.where((g.s<=p)&(g.e>=p),0,np.minimum(abs(g.s-p),abs(g.e-p)))
    g=g.sort_values('d').head(k); return ';'.join(f"{s}({d//1000}kb)" for s,d in zip(g.sym,g.d))
def genes_in(c,a,b):
    g=loc[(loc.chr==c)&(loc.e>=a)&(loc.s<=b)]; return ';'.join(g.sym)
G=['gwas_1','gwas_2','gwas_1_female','gwas_1_male','gwas_1_infectious_onset','gwas_1_non_infectious_onset']
D={g:pd.read_parquet(W+f'pq/{g}.parquet').set_index('ID') for g in G}
rows=[]
for g in G:
    d=D[g]; s=d[d.LOG10P>-np.log10(5e-8)].sort_values('LOG10P',ascending=False).copy()
    used=pd.Series(False,index=s.index)
    for i,r in s.iterrows():
        if used[i]: continue
        m=(s.CHROM==r.CHROM)&(abs(s.GENPOS-r.GENPOS)<=500_000)&(~used)
        blk=s[m]; used[m]=True
        a,b=blk.GENPOS.min(),blk.GENPOS.max()
        rows.append(dict(gwas=g,lead=i,chr=r.CHROM,pos=r.GENPOS,effect_allele=r.ALLELE1,eaf=r.A1FREQ,beta=r.BETA,se=r.SE,
            OR=np.exp(r.BETA),OR_lo=np.exp(r.BETA-1.96*r.SE),OR_hi=np.exp(r.BETA+1.96*r.SE),log10p=r.LOG10P,n_gw_snps=len(blk),
            region=f"{r.CHROM}:{a}-{b}",nearest_pc_genes=near(r.CHROM,r.GENPOS),pc_genes_in_sig_region=genes_in(r.CHROM,a,b)))
L=pd.DataFrame(rows)
# merge loci across GWAS (+/-500kb) into locus IDs
L=L.sort_values(['chr','pos']); lid=[];cur=None
for _,r in L.iterrows():
    if cur is None or r.chr!=cur[0] or r.pos-cur[1]>1_000_000: cur=[r.chr,r.pos,len(set(lid))+1]
    lid.append(cur[2])
L['locus']=lid
for g in G:
    L['log10p_'+g]=[D[g].LOG10P.get(i,np.nan) for i in L.lead]
L.to_csv(OUT+'lead_loci_all_gwas.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(L[['locus','gwas','lead','log10p','OR','nearest_pc_genes']].to_string())
# paper Table 3 SNPs
paper=['1:173846152:T:C','6:26239176:A:G','6:97984426:C:CA','13:53194927:GT:G','15:54866724:A:G','17:52183006:C:T','20:48914387:T:TA','6:26329103:A:T','20:47529913:C:T']
T=[]
for sid in paper:
    row={'snp':sid}
    for g in G:
        r=D[g].loc[sid] if sid in D[g].index else None
        row[g]=None if r is None else f"OR={np.exp(r.BETA):.3f} -log10p={r.LOG10P:.2f}"
    T.append(row)
# SUDS3 indel id differs; find the lead near 12:118202773
d=D['gwas_1']; x=d[(d.CHROM==12)&(abs(d.GENPOS-118202773)<50)]
print(x[['GENPOS','ALLELE0','ALLELE1','BETA','LOG10P']])
for sid in x.index:
    row={'snp':sid}
    for g in G:
        r=D[g].loc[sid] if sid in D[g].index else None
        row[g]=None if r is None else f"OR={np.exp(r.BETA):.3f} -log10p={r.LOG10P:.2f}"
    T.append(row)
T=pd.DataFrame(T); T.to_csv(OUT+'paper_table3_lookup.csv',index=False); print(T.to_string())
for g in G:
    d=D[g]; chi=(d.BETA/d.SE)**2
    print(g,'lambdaGC',round(np.median(chi)/0.4549,3),'N snps',len(d))
