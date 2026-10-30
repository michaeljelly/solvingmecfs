# Genetic correlation (numpy LDSC, same implementation as downloads/decodeme/scripts/09_cross_trait.py) of DecodeME GWAS-1
# with independent ME/CFS datasets and infection / haematology traits. N = effective N for binary traits.
import pandas as pd, numpy as np, glob, os, gzip
from scipy.stats import norm
W='/tmp/work_decodeme/'; G='/tmp/work_gmr/'; R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
ld=pd.concat([pd.read_csv(f,sep='\t',usecols=['CHR','SNP','BP','L2']) for f in sorted(glob.glob(W+'ldsc/eur_w_ld_chr/*.l2.ldscore.gz'))])
M=sum(float(open(f).read().split()[0]) for f in glob.glob(W+'ldsc/eur_w_ld_chr/*.l2.M_5_50'))
hm=set(ld.SNP)
comp={'A':'T','T':'A','C':'G','G':'C'}
def neff(ca,co): return 4/(1/ca+1/co)
def dme(g,ne):
    d=pd.read_parquet(G+f'dme_{g}.parquet',columns=['rs','ea','oa','beta','se']); d=d[d.rs.isin(hm)]
    return pd.DataFrame({'SNP':d.rs,'A1':d.ea,'A2':d.oa,'Z':d.beta/d.se,'N':ne})
def align(a,b):
    x=a.merge(b,on='SNP',suffixes=('1','2'))
    x=x[x.A11.str.len().eq(1)&x.A21.str.len().eq(1)&x.A12.str.len().eq(1)&x.A22.str.len().eq(1)]
    x=x[~(x.A11.map(comp)==x.A21)]
    same=(x.A11==x.A12)&(x.A21==x.A22); flip=(x.A11==x.A22)&(x.A21==x.A12)
    x=x[same|flip].copy(); x.loc[flip[same|flip].values,'Z2']*=-1
    return x.merge(ld,on='SNP')
def wreg(y,X,w):
    Xw=X*np.sqrt(w)[:,None]; yw=y*np.sqrt(w); return np.linalg.lstsq(Xw,yw,rcond=None)[0]
def ldsc_h2(z,N,l):
    chi=z**2; ok=chi<max(80,0.001*N.max())
    h=np.clip((chi.mean()-1)*M/(N.mean()*l.mean()),0,1)
    for _ in range(2):
        w=1/(np.maximum(l,1)*(1+N*h*l/M)**2)
        b=wreg(chi[ok],np.c_[N[ok]*l[ok]/M,np.ones(ok.sum())],w[ok]); h=np.clip(b[0],0,1)
    return b,w,ok
def rg(x):
    x=x.sort_values(['CHR','BP']).reset_index(drop=True)
    z1,z2,l=x.Z1.values,x.Z2.values,x.L2.values; N1,N2=x.N1.values.astype(float),x.N2.values.astype(float)
    b1,w1,ok1=ldsc_h2(z1,N1,l); b2,w2,ok2=ldsc_h2(z2,N2,l)
    h1,h2=b1[0],b2[0]; ok=ok1&ok2
    rho0=np.clip(np.mean(z1*z2)*M/(np.sqrt(N1*N2).mean()*l.mean()),-1,1)*np.sqrt(max(h1,1e-4)*max(h2,1e-4))
    w=1/(np.maximum(l,1)*((1+N1*h1*l/M)*(1+N2*h2*l/M)+(np.sqrt(N1*N2)*rho0*l/M)**2))
    blocks=np.array_split(np.arange(len(x)),200)
    def stats(y,X,wt):
        o=[]
        for bl in blocks:
            k=bl[ok[bl]]; Xw=X[k]*wt[k][:,None]; o.append((Xw.T@X[k], Xw.T@y[k]))
        return o
    S=[stats(z1**2,np.c_[N1*l/M,np.ones(len(l))],w1),stats(z2**2,np.c_[N2*l/M,np.ones(len(l))],w2),stats(z1*z2,np.c_[np.sqrt(N1*N2)*l/M,np.ones(len(l))],w)]
    def est(drop=None):
        r=[]
        for st in S:
            A=sum(t[0] for i,t in enumerate(st) if i!=drop); B=sum(t[1] for i,t in enumerate(st) if i!=drop); r.append(np.linalg.solve(A,B))
        return np.array([r[0][0],r[1][0],r[2][0],r[2][1],r[2][0]/np.sqrt(abs(r[0][0])*abs(r[1][0])),r[1][1]])
    full=est(); ps=np.array([est(i) for i in range(200)]); pv=200*full-199*ps
    se=np.sqrt(np.var(pv,axis=0,ddof=1)/200)
    return dict(n_snps=len(x),h2_dme=full[0],h2_dme_se=se[0],h2_other=full[1],h2_other_se=se[1],h2_other_z=full[1]/se[1],int_other=full[5],gcov_int=full[3],gcov_int_se=se[3],rg=full[4],rg_se=se[4])
D=dme('gwas_1',neff(15579,259909)); DI=dme('gwas_1_infectious_onset',neff(9738,259909))
key=pd.read_parquet(G+'dme_gwas_1.parquet',columns=['rs','chr','pos19']); key=key[key.rs.isin(hm)]; key['chr']=key.chr.astype(str)
T={}
m=pd.read_parquet(G+'mvp_only_derived.parquet'); m=m[m.rs.isin(hm)]
T['MVP-only ME/CFS (derived, 3,891 cases)']=pd.DataFrame({'SNP':m.rs,'A1':m.ea,'A2':m.oa,'Z':m.z,'N':m.Neff})
parts=[]
for c in pd.read_csv(G+'dl/fg_G6_POSTVIRFAT.tsv.gz',sep='\t',chunksize=2_000_000):
    c.columns=['chr','pos38','oa','ea','rs','pval','beta','se','eaf']; c['rs']=c.rs.astype(str).str.split(',').str[0]; parts.append(c[c.rs.isin(hm)])
c=pd.concat(parts); T['FinnGen G6_POSTVIRFAT (283 cases)']=pd.DataFrame({'SNP':c.rs,'A1':c.ea,'A2':c.oa,'Z':c.beta/c.se,'N':neff(283,463029)})
parts=[]
for c in pd.read_csv(G+'dl/fg_AB1_EBV.tsv.gz',sep='\t',chunksize=2_000_000):
    c.columns=['chr','pos38','oa','ea','rs','pval','beta','se','eaf']; c['rs']=c.rs.astype(str).str.split(',').str[0]; parts.append(c[c.rs.isin(hm)])
c=pd.concat(parts); T['FinnGen infectious mononucleosis (3,838 cases)']=pd.DataFrame({'SNP':c.rs,'A1':c.ea,'A2':c.oa,'Z':c.beta/c.se,'N':neff(3838,486087)})
parts=[]
for c in pd.read_csv(G+'dl/panukb_20002_1482.tsv.gz',sep='\t',dtype={'chr':str},chunksize=2_000_000):
    c.columns=['chr','pos19','oa','ea','af_case','eaf','beta','se','mlogp']; parts.append(c.merge(key,on=['chr','pos19']))
c=pd.concat(parts); T['Pan-UKB self-reported CFS 20002_1482 (EUR; shares UKB controls)']=pd.DataFrame({'SNP':c.rs,'A1':c.ea,'A2':c.oa,'Z':c.beta/c.se,'N':neff(1882,418591)})
Nmap={'EBVload_Kamitaki2026':654938,'HHV6Bload_Kamitaki2026':645249,'SleepDur_Dashti2019':446118,'met_IPA_Chen2023':8153,'met_ILA_Chen2023':8233,'met_KYN_Chen2023':8244,'met_TRP_Chen2023':8235,'RetCount_Vuckovic2020':408112,'IRF_Vuckovic2020':408112,'Hb_Chen2020':563946,'IM_FinnGen':None,'EBNA1_BL2020':7972,'VitD_Revez2020':None}
for eid,N in Nmap.items():
    f=G+f'exp/{eid}.hm3.tsv.gz'
    if N is None or not os.path.exists(f): continue
    c=pd.read_csv(f,sep='\t',dtype={'chr':str}); c['beta']=pd.to_numeric(c.beta,errors='coerce'); c['se']=pd.to_numeric(c.se,errors='coerce')
    c=c.dropna(subset=['beta','se']); c=c[c.se>0]
    T[eid]=pd.DataFrame({'SNP':c.rs,'A1':c.ea.str.upper(),'A2':c.oa.str.upper(),'Z':c.beta/c.se,'N':N})
res=[]
for t,d in T.items():
    for lab,dd in [('GWAS-1',D),('infection-onset',DI)]:
        if lab=='infection-onset' and not any(k in t for k in ['EBV','mononucleosis','MVP','POSTVIRFAT']): continue
        x=align(dd,d); r=rg(x); r['trait']=t; r['decodeme']=lab; res.append(r); print({k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items()},flush=True)
Rr=pd.DataFrame(res); Rr['rg_p']=2*norm.sf(abs(Rr.rg/Rr.rg_se)); Rr.to_csv(R+'ldsc_rg_replication_infection_haem.csv',index=False); print(Rr[['trait','decodeme','n_snps','h2_other_z','rg','rg_se','rg_p','gcov_int']].round(4).to_string())
