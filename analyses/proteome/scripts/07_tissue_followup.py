"""Follow-up of tissue-origin signals: which proteins drive them, secreted vs intracellular, activity/BMI confounding."""
import sys, os; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
import statsmodels.formula.api as smf
rd=lambda f: set(pd.read_csv(W+'raw/hpa/'+f,sep='\t').iloc[:,0].astype(str))
liver=rd('T_liver.tsv')|rd('C_hepatocytes.tsv'); muscle=rd('T_skeletal_muscle.tsv'); heart=rd('T_heart_muscle.tsv')|rd('C_cardiomyocytes.tsv')
plate=rd('C_platelets.tsv'); entero=rd('C_enterocytes.tsv')|rd('T_intestine.tsv'); kidney=rd('T_kidney.tsv')
loc=pd.read_csv(D+'hoel_gene_location.csv',index_col=0).iloc[:,0]
cc={k:pd.read_csv(R+f'cc_{v}.csv',index_col=0) for k,v in {'hoel':'hoel_serum_leakadj','germain':'germain_D1PRE','walitt':'walitt_serum','ukb':'ukb_olink','rydland_ev':'rydland_ev','walitt_csf':'walitt_csf'}.items()}
# 1. driver table
out=[]
for sname,S in [('liver',liver),('skeletal_muscle',muscle),('heart',heart),('platelet',plate),('enterocyte',entero),('kidney',kidney)]:
    genes=sorted(set().union(*[set(c.index)&S for c in cc.values()]))
    for g in genes:
        r=dict(set=sname,gene=g,location=loc.get(g,'?'))
        for k,c in cc.items():
            if g in c.index: r[f'd_{k}']=c.loc[g,'beta']; r[f'p_{k}']=c.loc[g,'p']
        out.append(r)
Dv=pd.DataFrame(out)
dcols=[c for c in Dv.columns if c.startswith('d_')]
Dv['mean_d_blood']=Dv[['d_hoel','d_germain','d_walitt','d_ukb']].mean(1); Dv['n_blood']=Dv[['d_hoel','d_germain','d_walitt','d_ukb']].notna().sum(1)
Dv.to_csv(R+'tissue_origin_drivers.csv',index=False)
for s,g in Dv.groupby('set'):
    g=g[g.n_blood>=2].sort_values('mean_d_blood')
    print('\n==',s,'n',len(g)); print(g.head(8)[['gene','location','mean_d_blood']+dcols].round(2).to_string(index=False)); print(g.tail(8)[['gene','location','mean_d_blood']+dcols].round(2).to_string(index=False))
# 2. secreted vs intracellular split for liver / muscle (competitive z)
rng=np.random.default_rng(3)
rows=[]
for k,c in cc.items():
    t=c.t.dropna()
    for sname,S in [('liver',liver),('skeletal_muscle',muscle),('platelet',plate),('enterocyte',entero)]:
        for L in ['Secreted','Intracellular','Membrane']:
            gg=[g for g in t.index if g in S and loc.get(g)==L]
            if len(gg)<4: continue
            obs=t[gg].mean(); null=np.array([t.values[rng.choice(len(t),len(gg),replace=False)].mean() for _ in range(3000)])
            rows.append(dict(dataset=k,set=sname,location=L,n=len(gg),mean_t=obs,z=(obs-null.mean())/null.std()))
S2=pd.DataFrame(rows); S2.to_csv(R+'tissue_origin_by_location.csv',index=False)
print(S2.pivot_table(index=['set','location'],columns='dataset',values='z').round(2).to_string())
