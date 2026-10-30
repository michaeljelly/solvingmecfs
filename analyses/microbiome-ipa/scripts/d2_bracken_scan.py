import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
R='/mnt/project-files/mecfs-research/datasets/repos/Microbe2991_Myalgic-Encephalomyelitis-Chronic-Fatigue-Syndrome-Microbiome-Study/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
br=pd.read_csv(R+'ME&CFS_Microbiome_Results/combined_bracken_output.txt',sep='\t')
smeta=pd.read_excel(R+'SampleMetadata.xlsx')
frac=[c for c in br.columns if c.endswith('.bracken_frac')]
X=br.set_index('name')[frac]; X.columns=[c.split('.')[0] for c in X.columns]
X.columns=[dict(zip(smeta['Name'],smeta['Sample Name'])).get(c,c) for c in X.columns]
X=X.T
meta=pd.read_csv(I+'metadata.csv',index_col=0); sub=meta.groupby('sample_id_tp1').first()
X=X.loc[[s for s in X.index if s in sub.index]]; info=sub.loc[X.index]
y=(info['study_ptorhc']=='MECFS').astype(int)
print('kingdom check: any viruses/archaea/fungi names present?')
for k in ['Methanobrevibacter smithii','Candida albicans','Saccharomyces cerevisiae','crAssphage','Escherichia phage','Faecalibacterium phage']:
    print(' ',k, k in X.columns)
print(' n taxa containing "virus":',sum('virus' in c for c in X.columns),
      '| "phage":',sum('phage' in c for c in X.columns),
      '| Candida:',sum(c.startswith('Candida ') for c in X.columns),
      '| Methano:',sum(c.startswith('Methano') for c in X.columns))
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
cols=[c for c in met.columns if c in meta.index]
L=np.log2(met[cols].astype(float).T); L['subj']=meta.loc[cols,'sample_id_tp1'].values
La=L.groupby('subj').mean(); z=lambda s:(s-s.mean())/s.std()
sco=lambda n: z(sum(z(La[x]) for x in n))
T=pd.DataFrame({'IPA':z(La['indolepropionate']),'ILA':z(La['indolelactate']),
                'fld':sco(['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine']),
                'arLDH':sco(['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate'])})
sh=[s for s in X.index if s in T.index]
F=X.loc[sh]; F=F.loc[:,(F>0).mean()>=0.5]
Fl=np.log10(F+1e-7)
print('\ntaxa with prevalence>=50%%: %d, n=%d samples'%(Fl.shape[1],len(sh)))
out=[]
for c in Fl.columns:
    d={'taxon':c}
    for t in T.columns:
        r,p=stats.spearmanr(Fl[c],T.loc[sh,t]); d[t+'_rho']=r; d[t+'_p']=p
    a=Fl.loc[y.loc[sh]==1,c]; b=Fl.loc[y.loc[sh]==0,c]
    d['pt_p']=stats.mannwhitneyu(a,b).pvalue; d['pt_dir']=np.sign(a.median()-b.median())
    out.append(d)
D=pd.DataFrame(out)
for t in list(T.columns)+['pt']:
    D[t+'_q']=multipletests(D[t+'_p'],method='fdr_bh')[1]
    print(f'{t}: hits q<0.1 = {(D[t+"_q"]<0.1).sum()} / {len(D)}; raw p<0.01 = {(D[t+"_p"]<0.01).sum()}')
D.to_csv(OUT+'bracken_taxon_scan.csv',index=False)
for t in ['IPA','fld','ILA','arLDH']:
    print(f'\n== top 12 taxa by |rho| with {t} (q<0.1 only)')
    s=D[D[t+'_q']<0.1].reindex(D[t+'_rho'].abs().sort_values(ascending=False).index).dropna(subset=['taxon'])
    print(s.head(12)[['taxon',t+'_rho',t+'_q']].to_string(index=False) if len(s) else '  none')
