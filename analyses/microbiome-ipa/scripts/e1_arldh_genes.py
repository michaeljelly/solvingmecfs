import pandas as pd, numpy as np, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/';B='/mnt/project-files/biomap/omics/Omics_Dataset/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
meta=pd.read_csv(I+'metadata.csv',index_col=0); met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
ko=pd.read_csv(B+'KEGG.csv',index_col=0); sp=pd.read_csv(B+'Specie Abundance.csv',index_col=0)
km=pd.read_csv(I+'feature_meta/kegg_feature_meta.csv',index_col=0)
cols=[c for c in met.columns if c in ko.columns and c in meta.index]
subj=meta.loc[cols,'sample_id_tp1']; z=lambda s:(s-s.mean())/s.std()
L=np.log2(met[cols].astype(float).T)
LAC=['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate','imidazole lactate']
FLD=['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine']
T=pd.DataFrame({m:z(L[m]) for m in LAC+FLD},index=cols)
T['arLDH4']=z(sum(z(L[m]) for m in LAC)); T['fld']=z(sum(z(L[m]) for m in FLD))
T['subj']=subj.values; A=T.groupby('subj').mean()
mdf=meta.groupby('sample_id_tp1').first().loc[A.index]
pt=(mdf['study_ptorhc']=='MECFS').astype(int)
cv=pd.DataFrame({'pt':pt,'age':mdf['age'].astype(float),'male':(mdf['gender']=='Male').astype(int),
 'bmi':pd.to_numeric(mdf['bmi'],errors='coerce')},index=A.index); cv['bmi']=cv['bmi'].fillna(cv['bmi'].median())

print('== does imidazole lactate join the aromatic-lactate module? (controls, subject level)')
print(A.loc[pt==0,LAC+FLD].corr(method='spearman').round(2).to_string())
print('\n== group differences (adj age/sex/BMI)')
for m in LAC+FLD+['arLDH4','fld']:
    r=sm.OLS(z(A[m]),sm.add_constant(cv)).fit()
    print(f'  {m:38s} beta={r.params["pt"]:+.3f} p={r.pvalues["pt"]:.3g}')

def prep(tab):
    X=tab[cols].astype(float).T; R=X.div(X.sum(axis=1),axis=0)
    Lg=np.log10(R+1e-7); Lg=Lg.sub(Lg.mean(axis=1),axis=0); Lg['subj']=subj.values
    return Lg.groupby('subj').mean()
K=prep(ko); S=prep(sp)

PANEL={'K00016':'ldh, L-lactate dehydrogenase (family of the Bifidobacterium aromatic LDH)',
 'K03778':'ldhA, D-lactate dehydrogenase','K00101':'lldD','K00102':'LDHD/dld','K03777':'dld (quinone)',
 'K18930':'dld','K00832':'tyrB, aromatic-amino-acid transaminase (ArAT, makes indolepyruvate)',
 'K01667':'tnaA, tryptophanase','K01745':'hutH, histidine ammonia-lyase','K11734':'aroP transporter',
 'K00812':'aspB/aspC aminotransferase','K00813':'aspB aminotransferase','K10797':'enr, 2-enoate reductase (acdA analogue)',
 'K00248':'bcd','K00074':'hbd'}
PANEL={k:v for k,v in PANEL.items() if k in K.columns}
print('\n== TARGETED GENE PANEL vs plasma targets (adj age/sex/BMI; beta per SD gene)')
rows=[]
for k,lab in PANEL.items():
    g=z(K[k])
    for t in ['arLDH4','indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate','imidazole lactate','fld','indolepropionate']:
        X=sm.add_constant(cv.drop(columns='pt').assign(g=g))
        r=sm.OLS(z(A[t]),X).fit()
        rows.append(dict(KO=k,gene=lab,target=t,beta=r.params['g'],p=r.pvalues['g']))
D=pd.DataFrame(rows); D['q']=multipletests(D['p'],method='fdr_bh')[1]
piv=D.pivot(index='KO',columns='target',values='beta').round(3)
pq=D.pivot(index='KO',columns='target',values='q')
ordc=['arLDH4','indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate','imidazole lactate','fld','indolepropionate']
print(piv[ordc].to_string())
print('\nhow many of the %d panel tests reach q<0.05: %d'%(len(D),(D.q<0.05).sum()))
print(D[D.q<0.05][['KO','gene','target','beta','q']].to_string(index=False) if (D.q<0.05).any() else '  none')
D.to_csv(OUT+'biomap_arldh_gene_panel.csv',index=False)

print('\n== Bifidobacterium / Lactobacillus species (BioMapAI 135-species table) vs the four lactates')
bif=[s for s in S.columns if s.startswith(('Bifidobacterium','Lactobacillus','Streptococcus'))]
for s in bif:
    out=[]
    for t in LAC+['arLDH4','indolepropionate']:
        r,p=stats.spearmanr(S[s],A[t]); out.append(f'{t.split()[0][:14]}={r:+.2f}({p:.2g})')
    print(f'  {s:32s} '+' '.join(out))
print('\n== whole-KO-table null for the 4-metabolite lactate score')
res=[]
for c in K.columns:
    r,p=stats.spearmanr(K[c],A['arLDH4']); res.append(p)
q=multipletests(res,method='fdr_bh')[1]
print(f'  KOs tested {len(res)}, q<0.05: {(q<0.05).sum()}, raw p<0.01: {(np.array(res)<0.01).sum()} (expected {0.01*len(res):.0f})')
