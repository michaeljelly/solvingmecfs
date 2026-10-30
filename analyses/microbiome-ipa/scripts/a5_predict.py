import pandas as pd, numpy as np, statsmodels.api as sm
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold
from scipy import stats
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
meta=pd.read_csv(I+'metadata.csv',index_col=0)
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
sp=pd.read_csv(B+'Specie Abundance.csv',index_col=0)
ko=pd.read_csv(B+'KEGG.csv',index_col=0)
cols=[c for c in met.columns if c in sp.columns and c in meta.index]
subj=meta.loc[cols,'sample_id_tp1']
z=lambda s:(s-s.mean())/s.std()
L=np.log2(met[cols].astype(float).T)
sco=lambda n: z(sum(z(L[x]) for x in n))
T=pd.DataFrame({'IPA':z(L['indolepropionate']),'ILA':z(L['indolelactate']),
                'fld':sco(['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine']),
                'arLDH':sco(['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate']),
                'butyrate':z(L['butyrate/isobutyrate (4:0)']),
                'pcresol':z(L['p-cresol sulfate']),'hippurate':z(L['hippurate']),
                'PAG':z(L['phenylacetylglutamine']),'kyn':z(L['kynurenine'])},index=cols)
def prep(tab):
    X=tab[cols].astype(float).T
    R=X.div(X.sum(axis=1),axis=0)
    Lg=np.log10(R+1e-7); return Lg.sub(Lg.mean(axis=1),axis=0)
S=prep(sp); K=prep(ko)
F=pd.concat([S.add_prefix('sp_'),K.add_prefix('ko_')],axis=1)
F=F.loc[:,F.std()>0]
grp=subj.values
gkf=GroupKFold(n_splits=5)
res=[]
for t in T.columns:
    yv=T[t].values; pred=np.zeros(len(yv))
    for tr,te in gkf.split(F,yv,grp):
        m=RidgeCV(alphas=np.logspace(0,5,20)).fit(F.values[tr],yv[tr])
        pred[te]=m.predict(F.values[te])
    r,p=stats.spearmanr(pred,yv)
    res.append(dict(target=t,cv_spearman=round(r,3),p=p,cv_R2=round(1-((yv-pred)**2).sum()/((yv-yv.mean())**2).sum(),3)))
    T[t+'_pred']=pred
R=pd.DataFrame(res).sort_values('cv_spearman',ascending=False)
print('== 5-fold group(person)-split ridge: stool species+KO -> plasma metabolite (n=%d samples, %d people)'%(len(F),len(set(grp))))
print(R.to_string(index=False))
R.to_csv(OUT+'biomap_cv_prediction.csv',index=False)
# does microbiome-predicted IPA explain the case-control gap?
Ta=T.copy(); Ta['subj']=subj.values; A=Ta.groupby('subj').mean()
mdf=meta.groupby('sample_id_tp1').first().loc[A.index]
pt=(mdf['study_ptorhc']=='MECFS').astype(int)
cv=pd.DataFrame({'pt':pt,'age':mdf['age'].astype(float),'male':(mdf['gender']=='Male').astype(int)},index=A.index)
print('\n== case-control effect before/after adjusting for microbiome-predicted value')
for t in ['IPA','fld','ILA','arLDH','butyrate']:
    r1=sm.OLS(z(A[t]),sm.add_constant(cv)).fit()
    X2=sm.add_constant(cv.assign(pred=z(A[t+'_pred'])))
    r2=sm.OLS(z(A[t]),X2).fit()
    rp=sm.OLS(z(A[t+'_pred']),sm.add_constant(cv)).fit()
    print(f'{t:9s} raw pt beta={r1.params["pt"]:+.3f} (p={r1.pvalues["pt"]:.3g}) -> adj={r2.params["pt"]:+.3f} (p={r2.pvalues["pt"]:.3g}); predicted-value pt beta={rp.params["pt"]:+.3f} (p={rp.pvalues["pt"]:.3g})')
