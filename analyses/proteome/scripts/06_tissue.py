"""Tissue/cell-of-origin analysis: are proteins enriched in a given tissue (HPA 'enriched' lists) shifted in patient blood?
Competitive test: mean gene t in the set vs 5000 random same-size gene sets from the same dataset. Replication across cohorts."""
import sys, glob, os; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
rng=np.random.default_rng(11)
sets={}
for f in glob.glob(W+'raw/hpa/*.tsv'):
    g=pd.read_csv(f,sep='\t').iloc[:,0].dropna().astype(str).tolist()
    if len(g)>=5: sets[os.path.basename(f)[:-4]]=set(g)
DS={'hoel_leakadj':'hoel_serum_leakadj','hoel_raw':'hoel_serum_raw','germain':'germain_D1PRE','walitt_serum':'walitt_serum','ukb':'ukb_olink',
    'walitt_csf':'walitt_csf','rydland_ev':'rydland_ev','igg':'igg_serum'}
rows=[]
for k,v in DS.items():
    t=pd.read_csv(R+f'cc_{v}.csv',index_col=0).t.dropna(); tv=t.values
    for s,g in sets.items():
        gg=[x for x in g if x in t.index]
        if len(gg)<5: continue
        obs=t[gg].mean()
        null=np.array([tv[rng.choice(len(tv),len(gg),replace=False)].mean() for _ in range(3000)])
        z=(obs-null.mean())/null.std(); p=(np.sum(np.abs(null-null.mean())>=abs(obs-null.mean()))+1)/3001
        rows.append(dict(dataset=k,set=s,n=len(gg),mean_t=obs,z=z,p=p))
T=pd.DataFrame(rows); T['q']=T.groupby('dataset').p.transform(lambda p: pd.Series(bh(p.values),index=p.index))
T.to_csv(R+'tissue_origin_by_dataset.csv',index=False)
Z=T.pivot(index='set',columns='dataset',values='z')
ind=['hoel_leakadj','germain','walitt_serum','ukb']
Z['stouffer_blood']=Z[ind].sum(1,min_count=2)/np.sqrt(Z[ind].notna().sum(1))
Z['p_blood']=2*stats.norm.sf(Z.stouffer_blood.abs()); Z['q_blood']=bh(Z.p_blood.values)
Z['n_same_sign_p05']=[max(sum((T[(T.set==s)&(T.dataset==d)].p.values<0.05)&(T[(T.set==s)&(T.dataset==d)].z.values>0)) for d in ind) for s in Z.index] if False else np.nan
Z.sort_values('p_blood').to_csv(R+'tissue_origin_meta.csv')
print(Z.sort_values('p_blood').round(2).to_string())
