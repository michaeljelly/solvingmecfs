import pandas as pd, numpy as np, warnings
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
C=pd.read_pickle(OUT+'contrasts.pkl'); meta=pd.read_pickle(OUT+'meta.pkl'); fm=pd.read_pickle(OUT+'fm.pkl'); L=pd.read_pickle(OUT+'L.pkl')
sub=pd.Series(fm['SUB PATHWAY'].values,index=fm.BIOCHEMICAL.astype(str)).reindex(L.columns)
# standardize each metabolite by SD of D1-PRE levels, so set scores are in SD units
sd=C['LEVEL D1-PRE'].std(); sd=sd.where(sd>0.05)
def zs(Y): return Y/sd[Y.columns]
# ---- lead sets
def has(s): return [c for c in L.columns if s.lower() in c.lower()]
S={}
for sp in sub.dropna().unique():
    m=list(sub[sub==sp].index)
    if len(m)>=3 and not str(sp).startswith('Drug') and sp not in ('Chemical','Food Component/Plant','Tobacco Metabolite'): S['SP: '+sp]=m
lead={
 'IPA (indolepropionate)':['indolepropionate'],
 'lactate':['lactate'],'pyruvate':['pyruvate'],'alanine':['alanine'],
 'alpha-ketoglutarate':['alpha-ketoglutarate'],'succinate':['succinate'],'citrate':['citrate'],'malate':['malate'],
 '3-hydroxybutyrate (ketone)':['3-hydroxybutyrate (BHBA)'],
 'carnitine (free)':['carnitine'],'acetylcarnitine':['acetylcarnitine (C2)'],
 'kynurenine':['kynurenine'],'tryptophan':['tryptophan'],'kynurenate':['kynurenate'],'quinolinate':['quinolinate'],'xanthurenate':['xanthurenate'],
 'VMA':['vanillylmandelate (VMA)'],'cortisol':['cortisol'],
 'hypoxanthine':['hypoxanthine'],'inosine':['inosine'],'xanthine':['xanthine'],'urate':['urate'],'AMP':['AMP'],
 'sphingosine':['sphingosine'],'sphingosine 1-phosphate':['sphingosine 1-phosphate'],
 'glutamate':['glutamate'],'glutamine':['glutamine'],'creatine':['creatine'],'creatinine':['creatinine'],
 '1,5-anhydroglucitol':['1,5-anhydroglucitol (1,5-AG)'],'glucose':['glucose'],'N-acetylarginine':['N-acetylarginine'],
 'ammonia-proxy: urea':['urea'],'citrulline':['citrulline'],'ornithine':['ornithine'],'arginine':['arginine'],
 'hypotaurine':['hypotaurine'],'taurine':['taurine'],'DHEA-S':['dehydroepiandrosterone sulfate (DHEA-S)'],
}
lead={k:[x for x in v if x in L.columns] for k,v in lead.items()}
print('missing leads:',[k for k,v in lead.items() if not v])
lead={k:v for k,v in lead.items() if v}
comp={
 'SET acylcarnitines (all)':[c for c in L.columns if 'carnitine' in c and sub.get(c,'') and 'Acyl Carnitine' in str(sub[c])],
 'SET long-chain free fatty acids':[c for c in L.columns if str(sub.get(c)) in ('Long Chain Saturated Fatty Acid','Long Chain Monounsaturated Fatty Acid','Long Chain Polyunsaturated Fatty Acid (n3 and n6)')],
 'SET primary bile acids':list(sub[sub=='Primary Bile Acid Metabolism'].index),
 'SET secondary bile acids':list(sub[sub=='Secondary Bile Acid Metabolism'].index),
 'SET kynurenine pathway':[c for c in ['kynurenine','kynurenate','quinolinate','xanthurenate','picolinate','N-acetylkynurenine (2)'] if c in L.columns],
 'SET purine degradation':[c for c in ['hypoxanthine','xanthine','inosine','urate'] if c in L.columns],
 'SET sphingolipids':list(sub[sub.isin(['Sphingosines','Sphingolipid Synthesis','Ceramide PEs'])].index),
 'SET TCA':list(sub[sub=='TCA Cycle'].index),
 'SET BCAA metabolism':list(sub[sub=='Leucine, Isoleucine and Valine Metabolism'].index),
 'SET catecholamine metabolites (tyrosine subpath)':[c for c in L.columns if c in ('vanillylmandelate (VMA)','dopamine 3-O-sulfate','dopamine 4-sulfate','homovanillate (HVA)','3-methoxytyramine sulfate','normetanephrine','metanephrine','vanillactate')],
}
print({k:len(v) for k,v in comp.items()}); print(comp['SET catecholamine metabolites (tyrosine subpath)'], comp['SET sphingolipids'])
allsets={**{'LEAD '+k:v for k,v in lead.items()},**comp,**S}
pd.to_pickle(allsets,OUT+'sets.pkl')
def score(Y,mem):
    mem=[m for m in mem if pd.notna(sd.get(m))]
    return zs(Y[mem]).mean(axis=1) if mem else pd.Series(np.nan,index=Y.index)
cov='male + AgeBin + BMIBin'
rows=[]
for k in ['R1 (D1 post-pre)','R2 (D2 post-pre)','R2-R1','REC (D2pre-D1pre)','D1post->D2pre','CUM (D2post-D1pre)','LEVEL D1-PRE']:
    Y=C[k]
    for s,mem in allsets.items():
        d=meta.copy(); d['y']=score(Y,mem); d=d.dropna(subset=['y'])
        r=smf.ols('y ~ case + '+cov,d).fit(cov_type='HC3')
        mc=d.y[d.case==0].mean(); mm=d.y[d.case==1].mean()
        pc=stats.ttest_1samp(d.y[d.case==0],0).pvalue; pm=stats.ttest_1samp(d.y[d.case==1],0).pvalue
        # sex interaction
        r2=smf.ols('y ~ case*male + AgeBin + BMIBin',d).fit(cov_type='HC3')
        # male-only case effect
        dm=d[d.male==1]; rm=smf.ols('y ~ case + AgeBin + BMIBin',dm).fit(cov_type='HC3')
        # Bell within patients
        dp=d[d.case==1]; rb=smf.ols('y ~ Bell + '+cov,dp).fit(cov_type='HC3')
        rho=stats.spearmanr(dp.Bell,dp.y)
        rows.append(dict(contrast=k,set=s,n_members=len(mem),mean_ctrl=mc,p_ctrl0=pc,mean_me=mm,p_me0=pm,diff=r.params['case'],p_diff=r.pvalues['case'],
                         sexint=r2.params['case:male'],p_sexint=r2.pvalues['case:male'],male_diff=rm.params['case'],p_male=rm.pvalues['case'],
                         bell_beta10=rb.params['Bell']*10,p_bell=rb.pvalues['Bell'],bell_rho=rho[0],p_bell_rho=rho[1]))
R=pd.DataFrame(rows)
for c in ['p_diff','p_sexint','p_male','p_bell']:
    def bh(p):
        q=pd.Series(np.nan,index=p.index); ok=p.notna(); q[ok]=multipletests(p[ok],method='fdr_bh')[1]; return q
    R[c.replace('p_','q_')]=R.groupby('contrast')[c].transform(bh)
R.to_csv(OUT+'set_scores_results.csv',index=False)
pd.set_option('display.width',250)
cols=['set','n_members','mean_ctrl','mean_me','diff','p_diff','q_diff','p_sexint','p_male','bell_beta10','p_bell','q_bell']
for k in R.contrast.unique():
    x=R[R.contrast==k]
    print('\n=====',k,' sets q_diff<0.1:',(x.q_diff<0.1).sum(),' q_sexint<0.1:',(x.q_sexint<0.1).sum(),' q_bell<0.1:',(x.q_bell<0.1).sum(), ' p_bell<.05:',(x.p_bell<.05).sum(),'of',len(x))
    print(x[x.set.str.startswith(('LEAD','SET'))][cols].round(3).to_string(index=False))
    print('-- top subpathways by p_diff'); print(x[x.set.str.startswith('SP')].sort_values('p_diff').head(8)[cols].round(3).to_string(index=False))
