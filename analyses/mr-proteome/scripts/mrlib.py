# Minimal TwoSampleMR-style toolkit in python (harmonise, IVW, MR-Egger, weighted median, Steiger, clumping via plink 1.9 + 1000G EUR)
import numpy as np, pandas as pd, subprocess, os, tempfile
from scipy import stats
PLINK='/tmp/work_gmr/ref/plink'; BFILE='/tmp/work_decodeme/magma/g1000_eur'
COMP={'A':'T','T':'A','C':'G','G':'C'}

def clump(df, p1=5e-8, r2=0.001, kb=10000, p2=1.0, snp='rs', pcol='p'):
    """df with rs and p; returns index SNPs present in 1000G EUR (others dropped, as in TwoSampleMR/OpenGWAS)."""
    d=df[[snp,pcol]].dropna(); d=d[d[pcol]<=p1]
    if len(d)==0: return []
    if len(d)==1: return d[snp].tolist()
    tmp=tempfile.mkdtemp(dir='/tmp/work_gmr')
    d.rename(columns={snp:'SNP',pcol:'P'}).to_csv(tmp+'/in.txt',sep='\t',index=False)
    subprocess.run([PLINK,'--bfile',BFILE,'--clump',tmp+'/in.txt','--clump-p1',str(p1),'--clump-p2',str(p2),'--clump-r2',str(r2),'--clump-kb',str(kb),'--out',tmp+'/o','--threads','1','--memory','1500'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    f=tmp+'/o.clumped'
    out=[] if not os.path.exists(f) else pd.read_csv(f,sep=r'\s+').SNP.tolist()
    subprocess.run(['rm','-rf',tmp]); return out

def ld_r(snps):
    """pairwise r (signed) matrix among snps from 1000G EUR; returns DataFrame"""
    tmp=tempfile.mkdtemp(dir='/tmp/work_gmr'); open(tmp+'/s.txt','w').write('\n'.join(snps))
    subprocess.run([PLINK,'--bfile',BFILE,'--extract',tmp+'/s.txt','--r','square','spaces','--write-snplist','--out',tmp+'/o','--threads','1','--memory','1500'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    try:
        s=open(tmp+'/o.snplist').read().split(); m=np.loadtxt(tmp+'/o.ld').reshape(len(s),len(s)); r=pd.DataFrame(m,index=s,columns=s)
    except Exception: r=None
    subprocess.run(['rm','-rf',tmp]); return r

def harmonise(e, o):
    """e,o: DataFrames with rs, ea, oa, eaf, beta, se. Returns merged with outcome aligned to exposure effect allele.
    Palindromic SNPs kept only if MAF<0.42 and eaf consistent (TwoSampleMR action 2)."""
    x=e.merge(o,on='rs',suffixes=('_e','_o'))
    for c in ['ea_e','oa_e','ea_o','oa_o']: x[c]=x[c].astype(str).str.upper()
    same=(x.ea_e==x.ea_o)&(x.oa_e==x.oa_o); swap=(x.ea_e==x.oa_o)&(x.oa_e==x.ea_o)
    fl=lambda s: s.map(lambda a: ''.join(COMP.get(c,'N') for c in a))
    same_f=(x.ea_e==fl(x.ea_o))&(x.oa_e==fl(x.oa_o)); swap_f=(x.ea_e==fl(x.oa_o))&(x.oa_e==fl(x.ea_o))
    pal=(x.ea_e.map(lambda a: COMP.get(a,'?'))==x.oa_e)
    sign=np.where(same|same_f,1.0,np.where(swap|swap_f,-1.0,np.nan))
    x['beta_o']=x.beta_o*sign
    if 'eaf_o' in x: x['eaf_o']=np.where(sign==-1,1-x.eaf_o,x.eaf_o)
    x=x[~np.isnan(sign)].copy()
    palx=pal[~np.isnan(sign)].values
    if 'eaf_e' in x and 'eaf_o' in x:
        maf=np.minimum(x.eaf_e,1-x.eaf_e)
        ambig=palx&((maf>0.42)|x.eaf_e.isna())
        # for palindromic with clear MAF, check that allele frequencies agree; if opposite sides of 0.5, flip outcome
        wrong=palx&~ambig&(((x.eaf_e-0.5)*(x.eaf_o-0.5))<0)
        x.loc[wrong,'beta_o']*=-1
        x=x[~ambig]
    else:
        x=x[~palx]
    return x.reset_index(drop=True)

def ivw(bx,by,sy):
    w=1/sy**2; b=np.sum(w*bx*by)/np.sum(w*bx**2); se_f=1/np.sqrt(np.sum(w*bx**2))
    k=len(bx); res=by-b*bx; Q=np.sum(w*res**2); phi=Q/(k-1) if k>1 else 1
    se=se_f*max(1,np.sqrt(phi)) if k>1 else se_f
    return b,se,Q,(stats.chi2.sf(Q,k-1) if k>1 else np.nan)

def egger(bx,by,sy):
    s=np.sign(bx); s[s==0]=1; bx=bx*s; by=by*s; w=1/sy**2
    X=np.c_[np.ones_like(bx),bx]; W=np.diag(w); XtW=X.T*w
    beta=np.linalg.solve(XtW@X,XtW@by); res=by-X@beta; k=len(bx)
    sig2=max(1,np.sum(w*res**2)/(k-2)); cov=np.linalg.inv(XtW@X)*sig2; se=np.sqrt(np.diag(cov))
    t=beta/se; p=2*stats.t.sf(np.abs(t),k-2)
    return beta[1],se[1],p[1],beta[0],se[0],p[0]

def wmedian_est(b,w):
    o=np.argsort(b); b=b[o]; w=w[o]/w.sum(); cs=np.cumsum(w)-0.5*w
    i=np.searchsorted(cs,0.5); 
    if i==0: return b[0]
    if i>=len(b): return b[-1]
    return b[i-1]+(b[i]-b[i-1])*(0.5-cs[i-1])/(cs[i]-cs[i-1])

def wmedian(bx,by,sx,sy,nboot=1000,seed=1):
    br=by/bx; sr=np.sqrt(sy**2/bx**2+by**2*sx**2/bx**4); w=1/sr**2
    est=wmedian_est(br,w); rng=np.random.default_rng(seed); bs=[]
    for _ in range(nboot):
        bxb=rng.normal(bx,sx); byb=rng.normal(by,sy); bs.append(wmedian_est(byb/bxb,w))
    return est,np.std(bs)

def r2_cont(beta,se,n):
    z2=(beta/se)**2; return z2/(z2+n-2)
def r2_binary(lor,eaf,ncase,ncont,prev):
    # liability-scale variance explained per SNP (as TwoSampleMR get_r_from_lor, logit model, Lee 2012 ascertainment ignored)
    vg=lor**2*2*eaf*(1-eaf); return vg/(vg+np.pi**2/3)

def mr_all(x, n_e=None, steiger_o=None, nboot=500):
    """x: harmonised frame. Returns dict of results."""
    bx,sx,by,sy=x.beta_e.values,x.se_e.values,x.beta_o.values,x.se_o.values
    k=len(x); r={'nsnp':k}
    if k==0: return r
    r['F_mean']=np.mean((bx/sx)**2)
    if k==1:
        b=by[0]/bx[0]; se=abs(sy[0]/bx[0]); r.update(method_main='Wald ratio',b=b,se=se,p=2*stats.norm.sf(abs(b/se)))
    else:
        b,se,Q,Qp=ivw(bx,by,sy); r.update(method_main='IVW',b=b,se=se,p=2*stats.norm.sf(abs(b/se)),Q=Q,Q_p=Qp)
    if k>=3:
        eb,ese,ep,ei,eise,eip=egger(bx.copy(),by.copy(),sy); r.update(egger_b=eb,egger_se=ese,egger_p=ep,egger_int=ei,egger_int_p=eip)
        wb,wse=wmedian(bx,by,sx,sy,nboot=nboot); r.update(wm_b=wb,wm_se=wse,wm_p=2*stats.norm.sf(abs(wb/wse)))
    if steiger_o is not None and n_e is not None:
        r2e=np.sum(r2_cont(bx,sx,n_e)); ncase,ncont,prev=steiger_o
        r2o=np.sum(r2_binary(by,x.eaf_o.values if 'eaf_o' in x else x.eaf_e.values,ncase,ncont,prev))
        # Steiger z-test on Fisher-transformed r
        ne=n_e; no=ncase+ncont
        z=(np.arctanh(np.sqrt(r2e))-np.arctanh(np.sqrt(r2o)))/np.sqrt(1/(ne-3)+1/(no-3))
        r.update(r2_exp=r2e,r2_out=r2o,steiger_correct=r2e>r2o,steiger_p=2*stats.norm.sf(abs(z)))
    return r
