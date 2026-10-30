"""Tissue-origin analysis redone with correlation-aware tests:
(a) self-contained: per-person mean z of set genes, OLS on case with covariates (Hoel, Germain, Walitt, EV);
(b) competitive CAMERA-style: z of mean t vs all genes, deflated by VIF = 1+(m-1)*mean inter-gene correlation (residual correlation from the same dataset;
    for UKB summary stats the correlation is borrowed from Germain D1PRE)."""
import sys, glob, os; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
loc=pd.read_csv(D+'hoel_gene_location.csv',index_col=0).iloc[:,0]
sets={}
for f in glob.glob(W+'raw/hpa/*.tsv'):
    g=set(pd.read_csv(f,sep='\t').iloc[:,0].dropna().astype(str))
    n=os.path.basename(f)[:-4]
    if len(g)>=5: sets[n]=g
    for L,ab in [('Intracellular','intra'),('Secreted','secr')]:
        gg={x for x in g if loc.get(x)==L}
        if len(gg)>=5: sets[f'{n}|{ab}']=gg
from themes import THEMES
for k,v in THEMES.items(): sets['theme_'+k]=set(v)
CFG={'hoel_leakadj':('hoel_serum',None,('male','age','bmi'),True,'hoel_serum_leakadj'),
     'germain':('germain_plasma',lambda m:m.timepoint=='D1PRE',('male','age','bmi_bin','site'),False,'germain_D1PRE'),
     'walitt_serum':('walitt_serum',None,('male','age'),False,'walitt_serum'),
     'walitt_csf':('walitt_csf',None,('male','age'),False,'walitt_csf'),
     'rydland_ev':('rydland_ev',None,(),False,'rydland_ev'),
     'lc_uk_explore':('lc_uk_explore',None,('male','age','bmi'),False,'lc_uk_explore'),
     'lc_iris_3mo':('lc_iris_3mo',None,('male',),False,'lc_iris_3mo'),
     'lc_bansal':('lc_bansal',None,('male','age'),False,'lc_bansal')}
def design(m,cols):
    des=m[list(cols)].copy()
    for c in list(des.columns):
        if des[c].dtype==object: des=des.drop(columns=c).join(pd.get_dummies(m[c],prefix=c,drop_first=True).astype(float))
    return des
rows=[]; RES={}
for name,(ds,sel,cov,leak,ccname) in CFG.items():
    X,m,f=load(ds)
    if sel is not None: k=sel(m).values; X=X[k]; m=m[k]
    if leak: m=m.join(pd.read_csv(D+f'{ds}_leak.csv',index_col=0).leak); cov=tuple(cov)+('leak',)
    G=gene_collapse(X,f); des=design(m,('case',)+tuple(cov)); ok=des.notna().all(1).values; des=des[ok]; G=G[ok].fillna(0)
    A=np.c_[np.ones(len(des)),des.values.astype(float)]; B=np.linalg.lstsq(A,G.values,rcond=None)[0]; Res=pd.DataFrame(G.values-A@B,columns=G.columns)
    RES[name]=Res
    t=pd.read_csv(R+f'cc_{ccname}.csv',index_col=0).t.dropna()
    for s,g in sets.items():
        gg=[x for x in g if x in G.columns and x in t.index]
        if len(gg)<5: continue
        sc=G[gg].mean(1).to_frame('s'); r=ols_many(sc,des,'case').iloc[0]
        C=np.corrcoef(Res[gg].values.T); mc=(C.sum()-len(gg))/(len(gg)*(len(gg)-1)); vif=1+(len(gg)-1)*mc
        others=t.drop(gg); diff=t[gg].mean()-others.mean()
        se=np.sqrt(vif*t.var()/len(gg)+t.var()/len(others)); zc=diff/se
        rows.append(dict(dataset=name,set=s,n=len(gg),d_selfcontained=r.beta,p_selfcontained=r.p,mean_t=t[gg].mean(),mean_corr=mc,vif=vif,z_camera=zc,p_camera=2*stats.norm.sf(abs(zc))))
# UKB with borrowed correlation from Germain residuals
u=pd.read_csv(R+'cc_ukb_olink.csv',index_col=0).t.dropna(); Rg=RES['germain']
for s,g in sets.items():
    gg=[x for x in g if x in u.index]
    if len(gg)<5: continue
    gc=[x for x in gg if x in Rg.columns]
    if len(gc)>=3:
        C=np.corrcoef(Rg[gc].values.T); mc=(C.sum()-len(gc))/(len(gc)*(len(gc)-1))
    else: mc=0.1
    vif=1+(len(gg)-1)*max(mc,0); others=u.drop(gg); diff=u[gg].mean()-others.mean(); se=np.sqrt(vif*u.var()/len(gg)+u.var()/len(others)); zc=diff/se
    rows.append(dict(dataset='ukb',set=s,n=len(gg),mean_t=u[gg].mean(),mean_corr=mc,vif=vif,z_camera=zc,p_camera=2*stats.norm.sf(abs(zc))))
T=pd.DataFrame(rows)
T['q_camera']=T.groupby('dataset').p_camera.transform(lambda p: pd.Series(bh(p.values),index=p.index))
T['q_self']=T.groupby('dataset').p_selfcontained.transform(lambda p: pd.Series(bh(p.values),index=p.index))
T.to_csv(R+'tissue_and_theme_camera.csv',index=False)
Z=T.pivot(index='set',columns='dataset',values='z_camera'); Dd=T.pivot(index='set',columns='dataset',values='d_selfcontained')
ind=['hoel_leakadj','germain','walitt_serum','ukb']
Z['stouffer_blood']=Z[ind].sum(1,min_count=2)/np.sqrt(Z[ind].notna().sum(1)); Z['p_blood']=2*stats.norm.sf(Z.stouffer_blood.abs()); Z['q_blood']=bh(Z.p_blood.values)
Z['n_coh_p05_same']=[max(((T[(T.set==s)&T.dataset.isin(ind)].p_camera<0.05)&(T[(T.set==s)&T.dataset.isin(ind)].z_camera>0)).sum(),((T[(T.set==s)&T.dataset.isin(ind)].p_camera<0.05)&(T[(T.set==s)&T.dataset.isin(ind)].z_camera<0)).sum()) for s in Z.index]
# self-contained meta over Hoel/Germain/Walitt
M=[]
for s in Z.index:
    x=T[(T.set==s)&T.dataset.isin(['hoel_leakadj','germain','walitt_serum'])].dropna(subset=['d_selfcontained'])
    se=(x.d_selfcontained/stats.norm.isf(x.p_selfcontained/2)).abs()
    r=meta_re(x.d_selfcontained.values,se.values) if len(x)>=2 else {}
    M.append(dict(set=s,sc_meta_d=r.get('re',np.nan),sc_meta_p=r.get('re_p',np.nan),sc_k=r.get('k',0)))
Z=Z.join(pd.DataFrame(M).set_index('set')); Z['sc_meta_q']=bh(Z.sc_meta_p.values)
Z.sort_values('p_blood').to_csv(R+'tissue_and_theme_camera_meta.csv')
pd.set_option('display.width',250)
print(Z.sort_values('p_blood').head(45).round(3).to_string())
