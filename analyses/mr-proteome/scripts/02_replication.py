# Part 1: independent replication of DecodeME in MVP-only (derived), FinnGen R12 G6_POSTVIRFAT, Pan-UKB 20002_1482 (EUR)
import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/genetics_mr/scripts'); import mrlib
W='/tmp/work_gmr/'; R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
import pyarrow.parquet as pq, pyarrow.dataset as ds, gc
D=pd.read_parquet(W+'dme_gwas_1.parquet'); D['p']=10**-D.mlogp
snv=D.ea.str.len().eq(1)&D.oa.str.len().eq(1)
leads_pos=[(1,173846152),(6,26239176),(6,97984426),(12,118202773),(13,53194927),(15,54866724),(17,52183006),(20,48914387)]
near=np.zeros(len(D),bool)
for c,pos in leads_pos: near|=(D.chr.astype(str)==str(c)).values&(abs(D.pos38-pos)<250000).values
D=D[((D.p<1e-3)|near)&snv].copy(); snv=pd.Series(True,index=D.index); gc.collect()
print('D subset',len(D),flush=True)
RS=set(D.rs)
M=pd.read_parquet(W+'mvp_only_derived.parquet',columns=['rs','ea','oa','eaf','beta','se']); M=M[M.rs.isin(RS)]
parts=[]
for c in pd.read_csv(W+'dl/fg_G6_POSTVIRFAT.tsv.gz',sep='\t',chunksize=2_000_000):
    c.columns=['chr','pos38','oa','ea','rs','pval','beta','se','eaf']; c['rs']=c.rs.astype(str).str.split(',').str[0]; parts.append(c[c.rs.isin(RS)])
fg=pd.concat(parts)[['rs','ea','oa','eaf','beta','se']]
key=D[['rs','chr','pos19']].copy(); key['chr']=key.chr.astype(str); parts=[]
for c in pd.read_csv(W+'dl/panukb_20002_1482.tsv.gz',sep='\t',dtype={'chr':str},chunksize=2_000_000):
    c.columns=['chr','pos19','oa','ea','af_case','eaf','beta','se','mlogp']; parts.append(c.merge(key,on=['chr','pos19']))
pu=pd.concat(parts); pu=pu[pu.eaf.between(0.005,0.995)][['rs','ea','oa','eaf','beta','se']]
REP={'MVP_only (3,891 cases; derived)':M,'FinnGen_R12_G6_POSTVIRFAT (283 cases)':fg,'PanUKB_20002_1482_EUR (UKB self-report CFS)':pu}
print({k:len(v) for k,v in REP.items()},flush=True)
Dx=D[['rs','ea','oa','eaf','beta','se','p']]
# --- A. 8 loci: best SNV per locus (within 250kb of lead) in DecodeME GWAS-1
leads=[('1q25.1 RABGAP1L',1,173846152),('6p22.2 BTN2A2',6,26239176),('6q16.1 FBXL4',6,97984426),('12q24.23 SUDS3',12,118202773),('13q14.3 OLFM4',13,53194927),('15q21.3 CCPG1',15,54866724),('17q22 CA10',17,52183006),('20q13.13 ARFGEF2/CSE1L',20,48914387)]
rows=[]
for name,c,pos in leads:
    s=D[(D.chr.astype(str)==str(c))&(abs(D.pos38-pos)<250000)&snv].sort_values('p').head(1).iloc[0]
    row=dict(locus=name,rs=s.rs,ea=s.ea,oa=s.oa,eaf=round(s.eaf,3),dme_OR=np.exp(s.beta),dme_p=s.p)
    for k,v in REP.items():
        h=mrlib.harmonise(Dx[Dx.rs==s.rs].rename(columns={'p':'pe'}),v[v.rs==s.rs])
        if len(h): 
            b,se=h.beta_o.iloc[0],h.se_o.iloc[0]; row[k+' OR']=np.exp(b); row[k+' p1side']=stats.norm.sf(b/se*np.sign(s.beta)); row[k+' z_aligned']=b/se*np.sign(s.beta)
    rows.append(row)
L=pd.DataFrame(rows); L.to_csv(R+'rep_8loci.csv',index=False); print(L.round(4).to_string(),flush=True)
# combined sign test across 8 loci per cohort
for k in REP:
    z=L[k+' z_aligned'].dropna(); print(k,'n',len(z),'same sign',(z>0).sum(),'binom p',stats.binomtest(int((z>0).sum()),len(z),0.5,alternative='greater').pvalue,'Stouffer z',z.sum()/np.sqrt(len(z)),'p1',stats.norm.sf(z.sum()/np.sqrt(len(z))))
# --- B. polygenic sign concordance at p thresholds (clumped r2<0.01, 1000 kb, 1000G EUR)
cl=mrlib.clump(D[snv][['rs','p']],p1=1e-3,r2=0.01,kb=1000)
C=D[D.rs.isin(cl)].copy(); print('clumped',len(C),flush=True)
out=[]
for thr in [5e-8,1e-6,1e-5,1e-4,1e-3]:
    s=C[C.p<thr]
    for k,v in REP.items():
        h=mrlib.harmonise(s[['rs','ea','oa','eaf','beta','se']],v)
        if len(h)<3: continue
        agree=np.sign(h.beta_e)==np.sign(h.beta_o); n=len(h)
        # IVW slope of replication beta on discovery beta (not an MR; a concordance slope)
        w=1/h.se_o**2; slope=np.sum(w*h.beta_e*h.beta_o)/np.sum(w*h.beta_e**2); sse=1/np.sqrt(np.sum(w*h.beta_e**2))
        out.append(dict(threshold=thr,cohort=k,n_snps=n,n_same_sign=int(agree.sum()),frac=agree.mean(),binom_p=stats.binomtest(int(agree.sum()),n,0.5,alternative='greater').pvalue,slope=slope,slope_se=sse,slope_p=2*stats.norm.sf(abs(slope/sse))))
S=pd.DataFrame(out); S.to_csv(R+'rep_sign_concordance.csv',index=False); print(S.round(4).to_string())
