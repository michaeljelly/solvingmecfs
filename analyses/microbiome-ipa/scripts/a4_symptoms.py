import pandas as pd, numpy as np, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
meta=pd.read_csv(I+'metadata.csv',index_col=0)
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
sp=pd.read_csv(B+'Specie Abundance.csv',index_col=0)
ko=pd.read_csv(B+'KEGG.csv',index_col=0)
sc=pd.read_csv(I+'score.csv',index_col=0)
z=lambda s:(s-s.mean())/s.std()

cols=[c for c in met.columns if c in sp.columns and c in meta.index]
md=meta.loc[cols]; subj=md['sample_id_tp1']
L=np.log2(met[cols].astype(float).T)
def sco(n): return z(sum(z(L[x]) for x in n))
F=pd.DataFrame({'arLDH':sco(['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate']),
                'fld':sco(['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine']),
                'IPA':z(L['indolepropionate']),'ILA':z(L['indolelactate'])},index=cols)
# microbial reductive-guild gene score
R=ko[cols].astype(float).T; R=R.div(R.sum(axis=1),axis=0)
Rg=np.log10(R+1e-9); Rg=Rg.sub(Rg.mean(axis=1),axis=0)
guild=[k for k in ['K10797','K00248','K00074','K00209','K01715'] if k in Rg.columns]
F['guild_gene']=z(sum(z(Rg[k]) for k in guild))
Sp=sp[cols].astype(float).T; Sp=Sp.div(Sp.sum(axis=1),axis=0)
Spl=np.log10(Sp+1e-6); Spl=Spl.sub(Spl.mean(axis=1),axis=0)
for s in ['Ruminococcus_gnavus','Flavonifractor_plautii','Akkermansia_muciniphila','Bifidobacterium_longum',
          'Bifidobacterium_adolescentis','Dysosmobacter_welbionis','Faecalibacterium_prausnitzii','Escherichia_coli',
          'Streptococcus_thermophilus','Clostridium_sp_AT4']:
    if s in Spl.columns: F[s]=z(Spl[s])
F['subj']=subj.values
Fa=F.groupby('subj').mean()
S=sc.copy(); S['subj']=meta.loc[S.index,'sample_id_tp1']
Sa=S.groupby('subj').mean(numeric_only=True)
mdf=md.groupby('sample_id_tp1').first()
common=sorted(set(Fa.index)&set(Sa.index))
Fa=Fa.loc[common]; Sa=Sa.loc[common]; mdf=mdf.loc[common]
pt=(mdf['study_ptorhc']=='MECFS')
cv=pd.DataFrame({'age':mdf['age'].astype(float),'male':(mdf['gender']=='Male').astype(int),
                 'bmi':pd.to_numeric(mdf['bmi'],errors='coerce')},index=common)
cv['bmi']=cv['bmi'].fillna(cv['bmi'].median())
print('subjects with symptoms+omics:',len(common),'patients:',pt.sum())

rows=[]
for f in Fa.columns:
    for s in Sa.columns:
        sub=pt[pt].index  # patients only
        v=Sa.loc[sub,s]; x=Fa.loc[sub,f]
        ok=v.notna()&x.notna()
        if ok.sum()<40: continue
        X=sm.add_constant(pd.concat([cv.loc[sub][ok],pd.DataFrame({'f':z(x[ok])})],axis=1))
        r=sm.OLS(z(v[ok]),X).fit()
        rows.append(dict(feature=f,symptom=s,beta=r.params['f'],p=r.pvalues['f'],n=int(ok.sum())))
D=pd.DataFrame(rows); D['q']=multipletests(D['p'],method='fdr_bh')[1]
D.to_csv(OUT+'biomap_symptom_assoc_patients.csv',index=False)
print('\n== patient-only symptom associations, q<0.1')
print(D[D.q<0.1].sort_values('p')[['feature','symptom','beta','p','q','n']].to_string(index=False))
print('\n== all with raw p<0.01')
print(D[D.p<0.01].sort_values('p')[['feature','symptom','beta','p','q']].head(25).to_string(index=False))

# illness duration
dur=mdf.loc[pt[pt].index,'illness_duration']
print('\n== illness duration (patients):',dur.value_counts().to_dict())
for f in Fa.columns:
    x=Fa.loc[pt[pt].index,f]; d=(dur=='Long').astype(float)
    X=sm.add_constant(pd.concat([cv.loc[pt[pt].index],pd.DataFrame({'long':d})],axis=1))
    r=sm.OLS(z(x),X).fit()
    if r.pvalues['long']<0.1: print(f'{f:32s} long-vs-short beta={r.params["long"]:+.2f} p={r.pvalues["long"]:.3g}')
