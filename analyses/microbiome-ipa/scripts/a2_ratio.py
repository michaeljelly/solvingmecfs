import pandas as pd, numpy as np, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
meta=pd.read_csv(I+'metadata.csv',index_col=0)
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
need=['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate','tryptophan','phenylalanine','tyrosine',
      'indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine','lactate','pyruvate','phenylpyruvate',
      'indoleacetate','kynurenine','p-cresol sulfate','phenylacetylglutamine']
cols=[c for c in met.columns if c in meta.index]
M=np.log2(met.loc[need,cols].astype(float).T)
md=meta.loc[cols]
M['subj']=md['sample_id_tp1']
agg=M.groupby('subj').mean(numeric_only=True)
mdf=md.groupby('sample_id_tp1').first().loc[agg.index]
y=(mdf['study_ptorhc']=='MECFS').astype(int)
cov=pd.DataFrame({'pt':y,'age':mdf['age'].astype(float),
                  'male':(mdf['gender']=='Male').astype(int),
                  'bmi':pd.to_numeric(mdf['bmi'],errors='coerce')})
cov['bmi']=cov['bmi'].fillna(cov['bmi'].median())

def test(v,extra=None,label=''):
    v=(v-v.mean())/v.std()
    X=cov.copy()
    if extra is not None:
        for k,s in extra.items(): X[k]=(s-s.mean())/s.std()
    X=sm.add_constant(X)
    r=sm.OLS(v,X).fit()
    return dict(test=label,beta_sd=round(r.params['pt'],3),p=r.pvalues['pt'],n=len(v))

rows=[]
pairs=[('indolelactate','tryptophan'),('phenyllactate (PLA)','phenylalanine'),
       ('3-(4-hydroxyphenyl)lactate','tyrosine'),('indolepropionate','tryptophan'),
       ('3-phenylpropionate (hydrocinnamate)','phenylalanine')]
for lac,aa in pairs:
    rows.append(test(agg[lac],None,f'{lac} raw'))
    rows.append(test(agg[lac],{'aa':agg[aa]},f'{lac} | {aa}'))
    rows.append(test(agg[lac]-agg[aa],None,f'log ratio {lac}/{aa}'))
# parent AAs themselves
for aa in ['tryptophan','phenylalanine','tyrosine']:
    rows.append(test(agg[aa],None,f'{aa} raw'))
# lactate/pyruvate general redox control
rows.append(test(agg['lactate']-agg['pyruvate'],None,'log lactate/pyruvate'))
rows.append(test(agg['indolelactate'],{'lacpyr':agg['lactate']-agg['pyruvate'],'aa':agg['tryptophan']},'indolelactate | Trp + lac/pyr'))
# composite scores
z=lambda s:(s-s.mean())/s.std()
arLDH=z(z(agg['indolelactate'])+z(agg['phenyllactate (PLA)'])+z(agg['3-(4-hydroxyphenyl)lactate']))
fld=z(z(agg['indolepropionate'])+z(agg['3-phenylpropionate (hydrocinnamate)'])+z(agg['cinnamoylglycine']))
rows.append(test(arLDH,None,'COMPOSITE aromatic-lactate score'))
rows.append(test(fld,None,'COMPOSITE fld-reductive score'))
rows.append(test(arLDH,{'fld':fld},'arLDH score | fld score'))
rows.append(test(fld,{'ar':arLDH},'fld score | arLDH score'))
res=pd.DataFrame(rows)
res['q']=multipletests(res['p'],method='fdr_bh')[1]
res['p']=res['p'].map(lambda x:f'{x:.2e}')
res['q']=res['q'].map(lambda x:f'{x:.3f}')
print(res.to_string(index=False))
res.to_csv(OUT+'biomap_lactate_vs_parent_aa.csv',index=False)
print('\ncorr arLDH vs fld score (all): %.2f'%arLDH.corr(fld,method='spearman'))
pd.DataFrame({'arLDH_score':arLDH,'fld_score':fld,'pt':y}).to_csv(OUT+'biomap_scores.csv')
