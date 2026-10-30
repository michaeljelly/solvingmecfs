import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
D='/tmp/work_micro/dl/esteban/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
tax=pd.read_csv(D+'taxonomy.csv').set_index('Feature_ID')
ab=pd.read_csv(D+'Relative_abudance_with_OTUID.csv',index_col=0)
met=pd.read_csv(D+'Metabolites_for_corr.csv',index_col=0)
grp=pd.read_csv(D+'metabolites_no_outliers.csv',encoding='utf-8-sig').rename(columns={'Sample':'Sample_ID'}).set_index('Sample_ID')['Group']
print('ab',ab.shape,'met',met.shape)
g=tax['Genus'].reindex(ab.columns).fillna('unassigned').str.replace('D_5__','',regex=False)
fam=tax['Family'].reindex(ab.columns).fillna('unassigned').str.replace('D_4__','',regex=False)
G=ab.T.groupby(g.values).sum().T
F=ab.T.groupby(fam.values).sum().T
common=sorted(set(G.index)&set(met.index))
print('paired samples',len(common))
G=G.loc[common]; F=F.loc[common]; met=met.loc[common]
z=lambda s:(s-s.mean())/s.std()
Mz=np.log2(met.astype(float).replace(0,np.nan))
targets={'IPA':'indole-3-propionate','ILA':'indole-3-lactate','PPA':'phenylpropionate','PLA':'phenyllactate',
         'OHPLA':'4-hydroxyphenyllactate','indole':'indole','butyrate':'butyrate'}
targets={k:v for k,v in targets.items() if v in Mz.columns}
def scan(T,lab,minprev=0.25):
    T=T.loc[:,(T>0).mean()>=minprev]
    Tl=np.log10(T.div(T.sum(axis=1),axis=0)+1e-5)
    out=[]
    for c in Tl.columns:
        for k,v in targets.items():
            ok=Mz[v].notna()
            if ok.sum()<25: continue
            r,p=stats.spearmanr(Tl.loc[ok,c],Mz.loc[ok,v])
            out.append(dict(taxon=c,target=k,rho=r,p=p,n=int(ok.sum())))
    O=pd.DataFrame(out)
    O['q']=O.groupby('target')['p'].transform(lambda s:multipletests(s,method='fdr_bh')[1])
    print(f'\n== {lab}: n taxa={Tl.shape[1]}; hits q<0.1 per target')
    print(O[O.q<0.1].groupby('target').size().to_dict())
    print(O[O.q<0.1].sort_values('p')[['taxon','target','rho','q','n']].head(20).to_string(index=False))
    O.to_csv(OUT+f'esteban_{lab}_metab_corr.csv',index=False)
    return O
scan(G,'genus'); scan(F,'family')
