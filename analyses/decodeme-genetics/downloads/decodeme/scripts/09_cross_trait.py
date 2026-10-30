# Cross-trait: (a) lookup of DecodeME leads in Long COVID (HGI, Lammi 2025) and FOXP4 in DecodeME;
# (b) minimal numpy LDSC (Bulik-Sullivan 2015) for h2 and genetic correlation, HapMap3 LD scores (eur_w_ld_chr), 200-block jackknife.
import pandas as pd, numpy as np, glob
W='/tmp/work_decodeme/'
m=pd.read_csv(W+'pq/id_rs_map.tsv.gz',sep='\t')
ld=pd.concat([pd.read_csv(f,sep='\t',usecols=['CHR','SNP','BP','L2']) for f in sorted(glob.glob(W+'ldsc/eur_w_ld_chr/*.l2.ldscore.gz'))])
M=sum(float(open(f).read().split()[0]) for f in glob.glob(W+'ldsc/eur_w_ld_chr/*.l2.M_5_50'))
hm=set(ld.SNP)
def dme(g):
    d=pd.read_parquet(W+f'pq/{g}.parquet',columns=['ID','ALLELE0','ALLELE1','BETA','SE','N','A1FREQ','N_CASES','N_CONTROLS'])
    d=d.merge(m[['ID','rs']],on='ID'); d=d[d.rs.isin(hm)]
    return pd.DataFrame({'SNP':d.rs,'A1':d.ALLELE1,'A2':d.ALLELE0,'Z':d.BETA/d.SE,'N':d.N,'MAF':d.A1FREQ,'ncase':d.N_CASES,'nctrl':d.N_CONTROLS})
def ext(path,cols,sep='\t',Nfix=None):
    it=pd.read_csv(path,sep=sep,usecols=list(cols.values()),chunksize=2_000_000,low_memory=False)
    parts=[]
    for c in it:
        c=c.rename(columns={v:k for k,v in cols.items()}); c=c[c.SNP.isin(hm)]; parts.append(c)
    d=pd.concat(parts); d['Z']=d.BETA/d.SE
    if Nfix: d['N']=Nfix
    return d[['SNP','A1','A2','Z','N']]
comp={'A':'T','T':'A','C':'G','G':'C'}
def align(a,b):
    x=a.merge(b,on='SNP',suffixes=('1','2'))
    x=x[x.A11.str.len().eq(1)&x.A21.str.len().eq(1)]
    amb=((x.A11.map(comp)==x.A21)); x=x[~amb]
    same=(x.A11==x.A12)&(x.A21==x.A22); flip=(x.A11==x.A22)&(x.A21==x.A12)
    x=x[same|flip].copy(); x.loc[flip[same|flip].values,'Z2']*=-1
    return x.merge(ld,on='SNP')
def wreg(y,X,w):
    Xw=X*np.sqrt(w)[:,None]; yw=y*np.sqrt(w); return np.linalg.lstsq(Xw,yw,rcond=None)[0]
def ldsc_h2(z,N,l):
    chi=z**2; ok=chi<max(80,0.001*N.max()); 
    h=np.clip((chi.mean()-1)*M/(N.mean()*l.mean()),0,1)
    for _ in range(2):
        w=1/(np.maximum(l,1)*(1+N*h*l/M)**2)
        b=wreg(chi[ok],np.c_[N[ok]*l[ok]/M,np.ones(ok.sum())],w[ok]); h=np.clip(b[0],0,1)
    return b,w,ok
def jk(fun,n,nb=200):
    idx=np.array_split(np.arange(n),nb); full=fun(np.ones(n,bool)); ps=[]
    for i in idx:
        k=np.ones(n,bool); k[i]=False; ps.append(fun(k))
    ps=np.array(ps); pv=nb*full-(nb-1)*ps; return full,np.sqrt(np.var(pv,axis=0,ddof=1)/nb)
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
        return np.array([r[0][0],r[1][0],r[2][0],r[2][1],r[2][0]/np.sqrt(abs(r[0][0])*abs(r[1][0]))])
    full=est(); ps=np.array([est(i) for i in range(200)]); pv=200*full-199*ps
    est,se=full,np.sqrt(np.var(pv,axis=0,ddof=1)/200)
    return dict(n_snps=len(x),h2_1=est[0],h2_1_se=se[0],h2_2=est[1],h2_2_se=se[1],gcov=est[2],gcov_int=est[3],rg=est[4],rg_se=se[4])
from scipy.stats import norm
res=[]
D=dme('gwas_1')
# liability-scale h2 for DecodeME (K=0.0065 approx UK prevalence per paper? report observed scale and liability with K=0.01)
traits={
 'LongCOVID_HGI_GCST90454541':(W+'lc/GCST90454541.h.tsv.gz',dict(SNP='rsid',A1='effect_allele',A2='other_allele',BETA='beta',SE='standard_error'),None,1096374),
 'MultisiteChronicPain_Johnston2019':(W+'other/chronic_pain-bgen.stats.gz',dict(SNP='SNP',A1='ALLELE1',A2='ALLELE0',BETA='BETA',SE='SE'),None,387649),
 'MDD_GCST90726344':(W+'other/GCST90726344.h.tsv.gz',dict(SNP='rsid',A1='effect_allele',A2='other_allele',BETA='beta',SE='standard_error'),None,514455),
}
for a,b in [('gwas_1_infectious_onset','gwas_1_non_infectious_onset'),('gwas_1_female','gwas_1_male')]:
    x=align(dme(a)[['SNP','A1','A2','Z','N']],dme(b)[['SNP','A1','A2','Z','N']]); r=rg(x); r['trait']=f'DecodeME {a} vs {b}'; res.append(r); print(r,flush=True)
for t,(p,cols,_,N) in traits.items():
    try:
        T=ext(p,cols,Nfix=N); x=align(D[['SNP','A1','A2','Z','N']],T); r=rg(x); r['trait']=t; res.append(r); print(r,flush=True)
    except Exception as e: print(t,'failed',e)
R=pd.DataFrame(res); R['rg_p']=2*norm.sf(abs(R.rg/R.rg_se)); R.to_csv(W+'results/ldsc_rg.csv',index=False); print(R.to_string())
