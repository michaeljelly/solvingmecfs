import pandas as pd, numpy as np, statsmodels.api as sm
from scipy import stats
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'

meta=pd.read_csv(I+'metadata.csv',index_col=0)
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
fm=pd.read_csv(I+'feature_meta/metabolome_feature_meta.csv',index_col=0)

PANEL={
 # fld reductive pathway (Trp arm)
 'indolepropionate':'fld_Trp','indoleacrylate':'fld_Trp',
 # fld reductive pathway (Phe arm) - same operon, different substrate
 '3-phenylpropionate (hydrocinnamate)':'fld_Phe','cinnamoylglycine':'fld_Phe',
 # aromatic lactate step (ArAT + aromatic lactate dehydrogenase)
 'indolelactate':'arLDH','phenyllactate (PLA)':'arLDH','3-(4-hydroxyphenyl)lactate':'arLDH',
 # other indole routes
 'indoleacetate':'other_indole','3-formylindole':'other_indole','6-hydroxyindole sulfate':'other_indole',
 'tryptophan betaine':'other_indole',
 # tryptophanase / host
 'tryptophan':'host','kynurenine':'host','kynurenate':'host','N-acetyltryptophan':'host',
 'serotonin':'host',
 # other microbial aromatics
 'p-cresol sulfate':'tyr_micro','phenylacetylglutamine':'phe_micro','hippurate':'phe_micro',
 '4-hydroxyhippurate':'phe_micro','3-hydroxyhippurate':'phe_micro',
 'indolin-2-one':'other_indole',
}
panel=[m for m in PANEL if m in met.index]
print('panel found', len(panel), 'missing', [m for m in PANEL if m not in met.index])

cols=[c for c in met.columns if c in meta.index]
M=np.log2(met.loc[panel,cols].astype(float).T)
md=meta.loc[cols]
M['subj']=md['sample_id_tp1']; M['grp']=md['study_ptorhc']
# subject-level mean
agg=M.groupby('subj').mean(numeric_only=True)
g=M.groupby('subj')['grp'].first()
cov=md.groupby('sample_id_tp1').first()[['age','gender','bmi']].loc[agg.index]
cov['male']=(cov['gender']=='Male').astype(int)
cov['bmi']=pd.to_numeric(cov['bmi'],errors='coerce')
cov['bmi']=cov['bmi'].fillna(cov['bmi'].median())
y=(g.loc[agg.index]=='MECFS').astype(int)
print('subjects', len(agg), 'patients', y.sum())

rows=[]
X0=sm.add_constant(pd.DataFrame({'pt':y,'age':cov['age'].astype(float),'male':cov['male'],'bmi':cov['bmi']}))
for m in panel:
    v=agg[m]
    v=(v-v.mean())/v.std()
    r=sm.OLS(v,X0).fit()
    # unadjusted cohen d
    a=agg.loc[y==1,m]; b=agg.loc[y==0,m]
    sd=np.sqrt(((len(a)-1)*a.var()+(len(b)-1)*b.var())/(len(a)+len(b)-2))
    d=(a.mean()-b.mean())/sd
    rows.append(dict(metabolite=m,module=PANEL[m],beta_sd=r.params['pt'],p=r.pvalues['pt'],
                     cohen_d=d,n_pt=len(a),n_hc=len(b)))
res=pd.DataFrame(rows).sort_values('p')
from statsmodels.stats.multitest import multipletests
res['q']=multipletests(res['p'],method='fdr_bh')[1]
res.to_csv(OUT+'biomap_aromatic_panel.csv',index=False)
print(res.to_string(index=False))

# module means
res['sig']=res['q']<0.05
print('\nby module (mean beta_sd):')
print(res.groupby('module')['beta_sd'].agg(['mean','count']))

# correlation structure among controls (subject level)
cc=agg.loc[y==0,panel].corr(method='spearman')
cc.to_csv(OUT+'biomap_panel_corr_controls.csv')
key=['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine','indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate','indoleacetate','tryptophan']
key=[k for k in key if k in cc.index]
print('\nSpearman corr (controls):'); print(cc.loc[key,key].round(2).to_string())
agg.assign(pt=y).to_csv(OUT+'biomap_panel_subject_values.csv')
