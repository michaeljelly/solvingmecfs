"""Theme-level tests: self-contained (per-person set score) and competitive (mean gene t vs random same-size sets)."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
from themes import THEMES
POSTHOC={'liver_hepatocyte','nuclear_leak_artefact','ECM_fibroblast'}
rng=np.random.default_rng(7)
CFG={ # name: (ds, selector, covars, leak, test)
 'hoel_leakadj':('hoel_serum',None,('male','age','bmi'),True,'case'),
 'hoel_raw':('hoel_serum',None,('male','age','bmi'),False,'case'),
 'germain_D1PRE':('germain_plasma',lambda m:m.timepoint=='D1PRE',('male','age','bmi_bin','site'),False,'case'),
 'walitt_serum':('walitt_serum',None,('male','age'),False,'case'),
 'walitt_csf':('walitt_csf',None,('male','age'),False,'case'),
 'rydland_ev':('rydland_ev',None,(),False,'case'),
 'giloteaux_ev_0h':('giloteaux_ev',lambda m:m.time=='0h',('tmt',),False,'case'),
 'igg_serum':('igg_serum',None,('batch',),False,'case'),
 'bragee_csf_POTS':('bragee_csf',None,('severity',),False,'pots'),
}
def design(m,cols):
    des=m[list(cols)].copy()
    for c in list(des.columns):
        if des[c].dtype==object:
            des=des.drop(columns=c).join(pd.get_dummies(m[c],prefix=c,drop_first=True).astype(float))
    return des
rows=[]
for name,(ds,sel,cov,leak,test) in CFG.items():
    X,m,f=load(ds)
    if sel is not None: k=sel(m).values; X=X[k]; m=m[k]
    if leak: m=m.join(pd.read_csv(D+f'{ds}_leak.csv',index_col=0).leak); cov=tuple(cov)+('leak',)
    G=gene_collapse(X,f)
    des=design(m,(test,)+tuple(cov)); ok=des.notna().all(1).values; des=des[ok]; G=G[ok]
    cc=pd.read_csv(R+f"cc_{ {'hoel_leakadj':'hoel_serum_leakadj','hoel_raw':'hoel_serum_raw'}.get(name,name)}.csv",index_col=0)
    tvals=cc.t.dropna()
    for th,genes in THEMES.items():
        g=[x for x in genes if x in G.columns]
        if len(g)<2: continue
        sc=G[g].mean(1).to_frame('s')
        r=ols_many(sc,des,test).iloc[0]
        gt=[x for x in g if x in tvals.index]; obs=tvals[gt].mean()
        null=np.array([tvals.sample(len(gt),random_state=int(rng.integers(1e9))).mean() for _ in range(2000)])
        pc=(np.sum(np.abs(null)>=abs(obs))+1)/(len(null)+1)
        rows.append(dict(dataset=name,theme=th,posthoc=th in POSTHOC,n_genes=len(g),genes=';'.join(g),d=r.beta,se=r.se,p=r.p,mean_t=obs,p_competitive=pc,
                         n1=int(des[test].sum()),n0=int((des[test]==0).sum())))
# UKB summary-level
u=pd.read_csv(R+'cc_ukb_olink.csv',index_col=0); tv=u.t.dropna()
for th,genes in THEMES.items():
    g=[x for x in genes if x in tv.index]
    if len(g)<2: continue
    obs=tv[g].mean(); null=np.array([tv.sample(len(g),random_state=int(rng.integers(1e9))).mean() for _ in range(2000)])
    pc=(np.sum(np.abs(null)>=abs(obs))+1)/(len(null)+1)
    rows.append(dict(dataset='ukb_olink',theme=th,posthoc=th in POSTHOC,n_genes=len(g),genes=';'.join(g),d=u.loc[g,'beta'].mean(),se=np.nan,p=np.nan,mean_t=obs,p_competitive=pc,n1=int(u.n_case.median()),n0=int(u.n_ctrl.median())))
T=pd.DataFrame(rows)
T['q_selfcontained']=T.groupby('dataset').p.transform(lambda p: pd.Series(bh(p.values),index=p.index))
T['q_competitive']=T.groupby('dataset').p_competitive.transform(lambda p: pd.Series(bh(p.values),index=p.index))
T.to_csv(R+'themes_by_dataset.csv',index=False)
piv=T.pivot(index='theme',columns='dataset',values='d').round(2)
pv=T.pivot(index='theme',columns='dataset',values='p')
star=pv.map(lambda x:'*' if x<0.05 else '') if hasattr(pv,'map') else pv.applymap(lambda x:'*' if x<0.05 else '')
print((piv.astype(str)+star.reindex_like(piv).fillna('')).to_string())
print('\ncompetitive mean t:'); print(T.pivot(index='theme',columns='dataset',values='mean_t').round(2).to_string())
print(T.pivot(index='theme',columns='dataset',values='p_competitive').round(3).to_string())
# meta across blood cohorts (self-contained d): hoel_leakadj, germain, walitt_serum
M=[]
for th,s in T[T.dataset.isin(['hoel_leakadj','germain_D1PRE','walitt_serum'])].groupby('theme'):
    r=meta_re(s.d.values,s.se.values); r['theme']=th; M.append(r)
M=pd.DataFrame(M).set_index('theme'); M['re_q']=bh(M.re_p)
M.to_csv(R+'themes_meta_blood.csv'); print(M.round(4).to_string())
