# Two-sample MR of many exposures on ME/CFS: DecodeME GWAS-1 (discovery), DecodeME infection-onset, MVP-only (independent replication outcome)
import pandas as pd, numpy as np, glob, os, subprocess, sys
from scipy import stats
sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/genetics_mr/scripts'); import mrlib
W='/tmp/work_gmr/'; R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
I=pd.read_parquet(W+'instruments_catalog_pqtl.parquet')
# streamed full-sumstat exposures
meta={'RHR_Zhu2019':('Resting heart rate (Zhu 2019, UKB)',458969,'risk_factor'),'VitD_Revez2020':('25(OH)D (Revez 2020, UKB)',496946,'risk_factor'),
 'SleepDur_Dashti2019':('Sleep duration (Dashti 2019, UKB)',446118,'risk_factor'),'EBVload_Kamitaki2026':('EBV DNA load in blood (Kamitaki 2026)',654938,'infection'),
 'HHV6Bload_Kamitaki2026':('HHV-6B DNA load in blood (Kamitaki 2026)',645249,'infection'),'EBVreadpos_Schmidt2026':('EBV read positivity (Schmidt 2026, UKB WGS)',360284,'infection'),
 'EBNA1_BL2020':('EBV EBNA-1 IgG level (Butler-Laporte 2020)',7972,'infection'),'VCAp18_BL2020':('EBV VCA p18 IgG level',8518,'infection'),'ZEBRA_BL2020':('EBV ZEBRA IgG level',8191,'infection'),
 'EAD_BL2020':('EBV EA-D IgG level',7763,'infection'),'EBVsero_BL2020':('EBV seropositivity',8735,'infection'),'CMVsero_BL2020':('CMV seropositivity',8735,'infection'),
 'CMVpp150_BL2020':('CMV pp150 IgG level',5136,'infection'),'CMVpp52_BL2020':('CMV pp52 IgG level',5681,'infection'),
 'IM_FinnGen':('Infectious mononucleosis (FinnGen R12, 3,838 cases)',489925,'infection'),'ChronicPain_Johnston2019':('Multisite chronic pain (Johnston 2019, UKB)',387649,'risk_factor'),
 'MDD_GCST90726344':('Major depression (GCST90726344)',514455,'risk_factor')}
chen_names={'IPA':'Indolepropionate','ILA':'Indolelactate','IAA':'Indoleacetate','KYN':'Kynurenine','TRP':'Tryptophan','KYNA':'Kynurenate','QUIN':'Quinolinate','XANTH':'Xanthurenate','BHB':'3-hydroxybutyrate','ACAC':'Acetoacetate','CAR':'Carnitine','C2':'Acetylcarnitine','C16':'Palmitoylcarnitine','GCA':'Glycocholate','CA':'Cholate','DCA':'Deoxycholate','GCDCA':'Glycochenodeoxycholate','TCA':'Taurocholate','GLCA':'Glycolithocholate','SERO':'Serotonin','LACT':'Lactate'}
dme=pd.read_parquet(W+'dme_gwas_1.parquet',columns=['rs','chr','pos19','pos38','ea','oa','eaf'])
streams=[]
import re
done=set(re.findall(r'^(\S+) rows \d+ kept',''.join(open(f).read() for f in glob.glob(W+'exp_fetch*.log')),flags=re.M))
for f in sorted(glob.glob(W+'exp/*.tsv')):
    if f.endswith('.mhc.tsv'): continue
    eid=os.path.basename(f)[:-4]
    if eid not in done: continue
    try: d=pd.read_csv(f,sep='\t',dtype=str)
    except Exception: continue
    if len(d)==0: continue
    for c in ['eaf','beta','se','p','n']: d[c]=pd.to_numeric(d[c],errors='coerce')
    if eid.split('_')[0] in ('RetCount','IRF','Hb','EBVreadpos'): continue   # used for rg only (MR uses catalog instruments) / no usable file
    d=d[d.p<(5e-6 if eid.startswith('met_') else 5e-8)]
    if eid.startswith('RHR'):  # BOLT file: many IDs are chr:pos_A_B (GRCh37); map by position
        d['pos']=pd.to_numeric(d.pos); k=dme[['rs','chr','pos19']].copy(); k['chr']=k.chr.astype(str)
        d=d.drop(columns='rs').merge(k,left_on=['chr','pos'],right_on=['chr','pos19'])
        d['eaf']=np.nan
    if eid.startswith('met_'):
        nm=eid.split('_')[1]; lab=chen_names[nm]+' (Chen 2023 full sumstats)'; n=8000; cat='metabolite_full'
    else: lab,n,cat=meta.get(eid,(eid,np.nan,'other'))
    d['n']=d.n.fillna(n)
    d=d[d.rs.astype(str).str.startswith('rs')]
    mk=lambda z,e,c,l: pd.DataFrame({'rs':z.rs,'ea':z.ea.str.upper(),'oa':z.oa.str.upper(),'eaf':z.eaf,'beta':z.beta,'se':z.se,'p':z.p,'n':z.n,'exposure':e,'category':c,'source':'full sumstats','label':l})
    streams.append(mk(d[d.p<5e-8],eid,cat,lab))
    if eid.startswith('met_'): streams.append(mk(d,eid+'_p5e-6',cat+'_p5e-6',lab+' [p<5e-6 instruments]'))
S=pd.concat(streams) if streams else pd.DataFrame()
A=pd.concat([I,S],ignore_index=True)
A=A[np.isfinite(A.beta)&np.isfinite(A.se)&(A.se>0)]
# fill missing other allele from DecodeME biallelic record (catalog lead SNPs report only the risk allele)
A=A.merge(dme[['rs','ea','oa']].rename(columns={'ea':'d_ea','oa':'d_oa'}),on='rs',how='left')
comp=str.maketrans('ACGT','TGCA')
miss=A.oa.isna()|(A.oa=='NA')|(A.oa=='NAN')
A['ea']=A.ea.astype(str)
m1=miss&(A.ea==A.d_ea); m2=miss&(A.ea==A.d_oa)
A.loc[m1,'oa']=A.loc[m1,'d_oa']; A.loc[m2,'oa']=A.loc[m2,'d_ea']
m3=miss&~m1&~m2&A.d_ea.notna()
ec=A.ea.str.translate(comp)
m3a=m3&(ec==A.d_ea); m3b=m3&(ec==A.d_oa)
A.loc[m3a,'oa']=A.loc[m3a,'d_oa'].str.translate(comp); A.loc[m3b,'oa']=A.loc[m3b,'d_ea'].str.translate(comp)
A=A[A.oa.notna()&(A.oa!='NA')].drop(columns=['d_ea','d_oa'])
print(A.groupby('category').exposure.nunique(),flush=True)
# reduced LD reference for fast clumping
allrs=A.rs.unique(); open(W+'allrs.txt','w').write('\n'.join(allrs))
if not os.path.exists(W+'ref/inst.bed') or os.path.getmtime(W+'ref/inst.bed')<os.path.getmtime(W+'instruments_catalog_pqtl.parquet') or True:
    subprocess.run([mrlib.PLINK,'--bfile',mrlib.BFILE,'--extract',W+'allrs.txt','--make-bed','--out',W+'ref/inst','--threads','1','--memory','1500'],stdout=subprocess.DEVNULL)
mrlib.BFILE=W+'ref/inst'
inref=set(pd.read_csv(W+'ref/inst.bim',sep='\t',header=None)[1])
# outcomes restricted to instrument SNPs
def load_out(path,cols=['rs','ea','oa','eaf','beta','se']):
    d=pd.read_parquet(path,columns=cols); return d[d.rs.isin(set(allrs))]
OUT={'DME_GWAS1':(load_out(W+'dme_gwas_1.parquet'),(15579,259909,0.005)),
     'DME_infection':(load_out(W+'dme_gwas_1_infectious_onset.parquet'),(9738,259909,0.003)),
     'MVP_only':(load_out(W+'mvp_only_derived.parquet'),(3891,439202,0.005))}
res=[]; used=[]
for eid,g in A.groupby('exposure'):
    g=g.sort_values('p').drop_duplicates('rs')
    if len(g)>1:
        gi=g[g.rs.isin(inref)]
        keep=mrlib.clump(gi[['rs','p']],p1=(5e-6 if eid.endswith('_p5e-6') else 5e-8),r2=0.001,kb=10000) if len(gi)>1 else gi.rs.tolist()
        g=g[g.rs.isin(keep)]
    lab,cat,src=g.label.iloc[0],g.category.iloc[0],g.source.iloc[0]
    for oname,(o,steig) in OUT.items():
        x=mrlib.harmonise(g[['rs','ea','oa','eaf','beta','se','n']],o)
        if len(x)==0: continue
        r=mrlib.mr_all(x,n_e=float(np.nanmedian(g.n)) if g.n.notna().any() else None,steiger_o=steig,nboot=300)
        r.update(exposure=eid,label=lab,category=cat,source=src,outcome=oname); res.append(r)
        if oname=='DME_GWAS1': used.append(x.assign(exposure=eid,label=lab))
Rz=pd.DataFrame(res)
Rz['OR']=np.exp(Rz.b); Rz['OR_lo']=np.exp(Rz.b-1.96*Rz.se); Rz['OR_hi']=np.exp(Rz.b+1.96*Rz.se)
Rz['fdr_in_category']=np.nan
for (c,o),g in Rz.groupby(['category','outcome']):
    from statsmodels.stats.multitest import multipletests
    Rz.loc[g.index,'fdr_in_category']=multipletests(g.p,method='fdr_bh')[1]
Rz.to_csv(R+'mr_all_results.csv',index=False)
pd.concat(used).to_csv(R+'mr_harmonised_instruments_GWAS1.csv.gz',index=False)
print(Rz[Rz.outcome=='DME_GWAS1'].groupby('category').apply(lambda g: pd.Series({'n':len(g),'p<0.05':(g.p<0.05).sum(),'fdr<0.05':(g.fdr_in_category<0.05).sum()})))
