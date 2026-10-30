# Approximate coloc.abf (Giambartolomei 2014; coloc v5 default priors) between DecodeME GWAS-1 loci and eQTLGen whole-blood cis-eQTLs (Vosa 2021, n~31,684)
import pandas as pd, numpy as np, gzip, sys
W='/tmp/work_decodeme/'
L=pd.read_csv(W+'results/lead_loci_all_gwas.csv')
L=L.sort_values('log10p',ascending=False).drop_duplicates('locus')
m=pd.read_csv(W+'pq/id_rs_map.tsv.gz',sep='\t')
g1=pd.read_parquet(W+'pq/gwas_1.parquet',columns=['ID','CHROM','GENPOS','A1FREQ','BETA','SE','LOG10P','ALLELE0','ALLELE1'])
leads=L.merge(m,left_on='lead',right_on='ID',how='left')
# region in hg19 from lead rs position; fallback: nearest mapped SNP
regs=[]
for _,r in L.iterrows():
    sub=g1[(g1.CHROM==r.chr)&(abs(g1.GENPOS-r.pos)<=500_000)].merge(m,on='ID')
    lead19=sub.iloc[(sub.GENPOS-r.pos).abs().argsort()[:1]].pos19.values[0]
    regs.append((r.locus,r.gwas,r.lead,int(r.chr),lead19,sub))
# stream eQTLGen full
want={}
for loc_,g,lead,c,p19,sub in regs: want.setdefault(c,[]).append((p19-500_000,p19+500_000,loc_))
import subprocess, io
cond=' || '.join(f'($3=={c} && $4>={a} && $4<={b})' for c,v in want.items() for a,b,l in v)
cmd=f"zcat {W}eqtl/cis-eQTLs_full_20180905.txt.gz | awk -F'\\t' 'NR>1 && ({cond})'"
raw=pd.read_csv(io.StringIO(subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout),sep='\t',header=None,
    names=['Pvalue','SNP','SNPChr','SNPPos','Zscore','AssessedAllele','OtherAllele','Gene','GeneSymbol','GeneChr','GenePos','NrCohorts','NrSamples','FDR'])
rows=[]
for c,v in want.items():
    for a,b,l in v:
        x=raw[(raw.SNPChr==c)&raw.SNPPos.between(a,b)]
        rows+=list(zip([l]*len(x),x.SNP,x.AssessedAllele,x.OtherAllele,x.Zscore,x.GeneSymbol,x.NrSamples,x.Pvalue))
E=pd.DataFrame(rows,columns=['locus','rs','eA','eO','Z','gene','N','P'])
E.to_csv(W+'results/eqtlgen_region_extract.tsv.gz',sep='\t',index=False)
def labf(z,V,W2): r=W2/(V+W2); return 0.5*(np.log(1-r)+r*z*z)
def lsum(x): mx=np.max(x); return mx+np.log(np.sum(np.exp(x-mx)))
p1,p2,p12=1e-4,1e-4,5e-6
out=[]
for loc_,g,lead,c,p19,sub in regs:
    sub=sub.dropna(subset=['rs'])
    e=E[E.locus==loc_]
    for gene,eg in e.groupby('gene'):
        if eg.P.min()>1e-5: continue
        mm=sub.merge(eg,on='rs')
        same=(mm.eA==mm.ALLELE1)&(mm.eO==mm.ALLELE0); flip=(mm.eA==mm.ALLELE0)&(mm.eO==mm.ALLELE1)
        mm=mm[same|flip].copy(); mm['Ze']=np.where(same[same|flip],mm.Z,-mm.Z)
        if len(mm)<50: continue
        maf=mm.A1FREQ.clip(0.01,0.99)
        # GWAS: case-control, prior sd 0.2 on log-OR
        l1=labf(mm.BETA/mm.SE,mm.SE**2,0.2**2)
        # eQTL: standardised beta from Z, prior sd 0.15
        se_e=1/np.sqrt(2*maf*(1-maf)*(mm.N+mm.Ze**2)); b_e=mm.Ze*se_e
        l2=labf(mm.Ze,se_e**2,0.15**2)
        L1,L2,L12=lsum(l1),lsum(l2),lsum(l1+l2)
        d=min(L12-L1-L2,-1e-12)
        lH=np.array([0,np.log(p1)+L1,np.log(p2)+L2,np.log(p1)+np.log(p2)+L1+L2+np.log(-np.expm1(d)),np.log(p12)+L12])
        pp=np.exp(lH-lsum(lH))
        top=mm.iloc[np.argmax(l1+l2)]
        # direction: sign of eQTL Z at GWAS lead-ish SNP (max GWAS labf)
        tg=mm.iloc[np.argmax(l1)]
        out.append(dict(locus=loc_,gwas_lead=lead,gene=gene,n_snps=len(mm),min_eqtl_p=eg.P.min(),PP_H3=pp[3],PP_H4=pp[4],
            top_shared_snp=top.rs,gwas_top_snp=tg.rs,eqtl_z_at_gwas_top=tg.Ze,gwas_beta_at_top=tg.BETA,
            risk_allele_raises_expression=bool(np.sign(tg.Ze)==np.sign(tg.BETA))))
R=pd.DataFrame(out).sort_values(['locus','PP_H4'],ascending=[True,False])
R.to_csv(W+'results/coloc_eqtlgen_blood.csv',index=False)
print(R[R.PP_H4>0.5].to_string()); print(len(R),'gene-locus tests')
