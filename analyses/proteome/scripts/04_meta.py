"""Cross-cohort replication and meta-analysis at gene level."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
# UK Biobank Olink (Beentjes 2025, summary): convert z to d
u=pd.read_csv(D+'ukb_olink_TE.csv')
u['z']=u.res_combined_te/u.se_combined_te
k=np.sqrt(1/u.cases_combined+1/u.controls_combined)
ukb=pd.DataFrame({'beta':(u.z*k).values,'se':k.values,'p':u.pval_combined_te.values,'n_case':u.cases_combined.values,'n_ctrl':u.controls_combined.values},index=u.name.astype(str).values)
ukb=ukb[~ukb.index.duplicated()]; ukb['q']=bh(ukb.p); ukb['t']=ukb.beta/ukb.se
ukb.index.name='gene'; ukb.to_csv(R+'cc_ukb_olink.csv')
print('UKB q<0.05:',(ukb.q<0.05).sum()); print(ukb.sort_values('p').head(25).round(4).to_string())
blood={'hoel':'hoel_serum_leakadj','germain':'germain_D1PRE','walitt':'walitt_serum','ukb':'ukb_olink'}
other={'hoel_raw':'hoel_serum_raw','walitt_csf':'walitt_csf','rydland_ev':'rydland_ev','giloteaux_ev':'giloteaux_ev_0h','igg':'igg_serum'}
T={k:pd.read_csv(R+f'cc_{v}.csv',index_col=0) for k,v in {**blood,**other}.items()}
# 1. global concordance of t statistics
names=list(T); rows=[]
for i,a in enumerate(names):
    for b in names[i+1:]:
        g=T[a].index.intersection(T[b].index)
        if len(g)<30: continue
        r=stats.spearmanr(T[a].loc[g,'t'],T[b].loc[g,'t'])
        # among genes nominal in a, sign agreement in b
        s=T[a].loc[g]; s=s[s.p<0.05]
        agree=(np.sign(s.beta)==np.sign(T[b].loc[s.index,'beta'])).mean() if len(s) else np.nan
        rows.append(dict(a=a,b=b,n_shared=len(g),rho_t=r.correlation,p=r.pvalue,n_nom_a=len(s),sign_agree_in_b=agree))
conc=pd.DataFrame(rows); conc.to_csv(R+'concordance_pairs.csv',index=False); print(conc.round(3).to_string())
# 2. meta-analysis over blood cohorts
genes=sorted(set().union(*[T[k].index for k in blood]))
out=[]
for g in genes:
    d=[];s=[];ks=[]
    for k in blood:
        if g in T[k].index and not np.isnan(T[k].loc[g,'beta']): d.append(T[k].loc[g,'beta']); s.append(T[k].loc[g,'se']); ks.append(k)
    if len(d)<2: continue
    r=meta_re(d,s); r['gene']=g; r['cohorts']=','.join(ks)
    for k in blood:
        r[f'd_{k}']=T[k].loc[g,'beta'] if g in T[k].index else np.nan
        r[f'p_{k}']=T[k].loc[g,'p'] if g in T[k].index else np.nan
    dd=np.array(d); r['n_nominal']=int(sum(T[k].loc[g,'p']<0.05 for k in ks)); r['same_sign']=int(max((dd>0).sum(),(dd<0).sum()))
    out.append(r)
M=pd.DataFrame(out).set_index('gene'); M['re_q']=bh(M.re_p); M['fe_q']=bh(M.fe_p)
# replication: nominal p<0.05 same direction in >=2 cohorts
def rep(row):
    ps=[(row[f'd_{k}'],row[f'p_{k}']) for k in blood if not np.isnan(row[f'p_{k}'])]
    up=sum(1 for d,p in ps if p<0.05 and d>0); dn=sum(1 for d,p in ps if p<0.05 and d<0)
    return max(up,dn) if min(up,dn)==0 else 0
M['n_rep_same_dir']=M.apply(rep,axis=1)
M.sort_values('re_p').to_csv(R+'meta_blood_genes.csv')
print('genes meta:',len(M),' RE q<0.05:',(M.re_q<0.05).sum(),' FE q<0.05:',(M.fe_q<0.05).sum(),' replicated (>=2 cohorts p<0.05 same dir):',(M.n_rep_same_dir>=2).sum())
cols=['k','re','re_p','re_q','I2','d_hoel','d_germain','d_walitt','d_ukb','p_hoel','p_germain','p_walitt','p_ukb','n_rep_same_dir']
print(M.sort_values('re_p').head(40)[cols].round(4).to_string())
print(M[M.n_rep_same_dir>=2].sort_values('re_p')[cols].round(4).to_string())
# expected number replicated by chance: permutation of signs across cohorts
rng=np.random.default_rng(1); exp=[]
P={k:T[k][['beta','p']] for k in blood}
for it in range(200):
    cnt=0
    sh={k:P[k].sample(frac=1,random_state=int(rng.integers(1e9))).set_axis(P[k].index) for k in blood}
    tab=pd.concat({k:sh[k] for k in blood},axis=1)
    up=sum(((tab[(k,'p')]<0.05)&(tab[(k,'beta')]>0)).astype(int) for k in blood)
    dn=sum(((tab[(k,'p')]<0.05)&(tab[(k,'beta')]<0)).astype(int) for k in blood)
    exp.append(int((((up>=2)&(dn==0))|((dn>=2)&(up==0))).sum()))
print('replicated-by-chance null: mean',np.mean(exp),'95th pct',np.percentile(exp,95))
pd.Series(exp).to_csv(R+'meta_blood_rep_null.csv',index=False)
