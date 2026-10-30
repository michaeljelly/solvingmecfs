import pandas as pd, numpy as np, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
meta=pd.read_csv(I+'metadata.csv',index_col=0)
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
sp=pd.read_csv(B+'Specie Abundance.csv',index_col=0)
ko=pd.read_csv(B+'KEGG.csv',index_col=0)
kmeta=pd.read_csv(I+'feature_meta/kegg_feature_meta.csv',index_col=0)

cols=[c for c in met.columns if c in sp.columns and c in meta.index]
print('paired samples',len(cols))
md=meta.loc[cols]; subj=md['sample_id_tp1']
z=lambda s:(s-s.mean())/s.std()
L=np.log2(met[cols].astype(float).T)
def score(names): return z(sum(z(L[n]) for n in names))
arl=score(['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate'])
fld=score(['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine'])
tgt=pd.DataFrame({'arLDH':arl,'fld':fld,'IPA':z(L['indolepropionate']),'ILA':z(L['indolelactate'])},index=cols)
tgt['subj']=subj.values
tgt=tgt.groupby('subj').mean()

def prep(tab,minprev=0.1):
    X=tab[cols].astype(float).T
    X['subj']=subj.values
    X=X.groupby('subj').mean()
    prev=(X>0).mean()
    X=X.loc[:,prev>=minprev]
    # CLR-ish: log of relative abundance with pseudocount
    R=X.div(X.sum(axis=1),axis=0)
    Lg=np.log10(R+ (R[R>0].min().min()/2 if (R>0).any().any() else 1e-6))
    Lg=Lg.sub(Lg.mean(axis=1),axis=0)
    return Lg

S=prep(sp); K=prep(ko)
print('species',S.shape,'KO',K.shape)
mdf=md.groupby('sample_id_tp1').first().loc[S.index]
y=(mdf['study_ptorhc']=='MECFS').astype(int)
cv=pd.DataFrame({'age':mdf['age'].astype(float),'male':(mdf['gender']=='Male').astype(int),
                 'bmi':pd.to_numeric(mdf['bmi'],errors='coerce')})
cv['bmi']=cv['bmi'].fillna(cv['bmi'].median())
tgt=tgt.loc[S.index]

def scan(F,name):
    out=[]
    for c in F.columns:
        v=z(F[c])
        X=sm.add_constant(pd.concat([cv,pd.DataFrame({'pt':y})],axis=1))
        r=sm.OLS(v,X).fit(); pt_b,pt_p=r.params['pt'],r.pvalues['pt']
        d={}
        for t in ['arLDH','fld','IPA','ILA']:
            X2=sm.add_constant(pd.concat([cv,tgt[[t]]],axis=1))
            r2=sm.OLS(v,X2).fit(); d[t+'_b']=r2.params[t]; d[t+'_p']=r2.pvalues[t]
        out.append(dict(feature=c,pt_beta=pt_b,pt_p=pt_p,**d))
    D=pd.DataFrame(out)
    for col in ['pt_p','arLDH_p','fld_p','IPA_p','ILA_p']:
        D[col.replace('_p','_q')]=multipletests(D[col],method='fdr_bh')[1]
    D.to_csv(OUT+f'biomap_{name}_assoc.csv',index=False)
    return D

DS=scan(S,'species'); DK=scan(K,'kegg')
print('\n== SPECIES: patient vs control (q<0.1)')
print(DS[DS.pt_q<0.1].sort_values('pt_beta')[['feature','pt_beta','pt_q']].to_string(index=False))
for t in ['arLDH','fld','IPA','ILA']:
    s=DS[DS[t+'_q']<0.05].reindex(DS[t+'_b'].abs().sort_values(ascending=False).index).dropna(subset=['feature'])
    print(f'\n== SPECIES vs {t} (q<0.05, top 15)')
    print(s.head(15)[['feature',t+'_b',t+'_q','pt_beta','pt_q']].to_string(index=False))

DK['desc']=DK['feature'].map(kmeta.iloc[:,0] if kmeta.shape[1]>0 else {})
print('\n== KEGG vs arLDH/fld (q<0.05 counts)')
for t in ['arLDH','fld','IPA','ILA','pt']:
    print(t, (DK[t+'_q']<0.05).sum(), 'of', len(DK))
