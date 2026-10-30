import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
D='/tmp/work_micro/dl/esteban/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
m=pd.read_csv(D+'metabolites_no_outliers.csv',encoding='utf-8-sig')
m=m.rename(columns={'Sample':'Sample_ID'}).set_index('Sample_ID')
grp=m.pop('Group')
print('n',len(m),grp.value_counts().to_dict())
M=np.log2(m.astype(float).replace(0,np.nan))
z=lambda s:(s-s.mean())/s.std()

FLD=['indole-3-propionate','phenylpropionate','cinnamoylglycine','4-hydroxyphenylpropionate','3-hydroxyphenylpropionate']
ARL=['indole-3-lactate','phenyllactate','4-hydroxyphenyllactate']
AA=['tryptophan','phenylalanine','tyrosine']
keep=[c for c in FLD+ARL+AA+['indole','indoleacetate','kynurenine','kynurenate','butyrate','propionate','acetate','p-cresol','indoxyl sulfate','tryptamine','phenylpyruvate'] if c in M.columns]
print('missingness (frac NaN):'); print(M[keep].isna().mean().round(2).to_string())

rows=[]
for c in keep:
    a=M.loc[grp=='patient',c].dropna(); b=M.loc[grp=='control',c].dropna()
    if len(a)<8 or len(b)<8: continue
    sd=np.sqrt(((len(a)-1)*a.var()+(len(b)-1)*b.var())/(len(a)+len(b)-2))
    d=(a.mean()-b.mean())/sd
    u=stats.mannwhitneyu(a,b)
    rows.append(dict(metabolite=c,module=('fld' if c in FLD else 'arLDH' if c in ARL else 'aa' if c in AA else 'other'),
                     cohen_d=round(d,2),p=u.pvalue,n_pt=len(a),n_hc=len(b)))
R=pd.DataFrame(rows); R['q']=multipletests(R['p'],method='fdr_bh')[1]
R=R.sort_values('p')
print('\n== Esteban stool metabolome, patient vs control')
print(R.to_string(index=False))
R.to_csv(OUT+'esteban_stool_metabolites.csv',index=False)

# module correlation structure in controls
sub=[c for c in FLD+ARL if c in M.columns and M[c].isna().mean()<0.5]
cc=M.loc[grp=='control',sub].corr(method='spearman')
print('\n== Spearman among fld/arLDH metabolites (controls, stool)')
print(cc.round(2).to_string())
cc.to_csv(OUT+'esteban_module_corr.csv')
# cross-module mean r
f=[c for c in sub if c in FLD]; a=[c for c in sub if c in ARL]
print('\nmean within-fld r %.2f | within-arLDH r %.2f | cross r %.2f'%(
  cc.loc[f,f].values[np.triu_indices(len(f),1)].mean(),
  cc.loc[a,a].values[np.triu_indices(len(a),1)].mean(),
  cc.loc[f,a].values.mean()))
