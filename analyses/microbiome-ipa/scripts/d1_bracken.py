import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
R='/mnt/project-files/mecfs-research/datasets/repos/Microbe2991_Myalgic-Encephalomyelitis-Chronic-Fatigue-Syndrome-Microbiome-Study/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
br=pd.read_csv(R+'ME&CFS_Microbiome_Results/combined_bracken_output.txt',sep='\t')
smeta=pd.read_excel(R+'SampleMetadata.xlsx')
frac=[c for c in br.columns if c.endswith('.bracken_frac')]
X=br.set_index('name')[frac]
X.columns=[c.split('.')[0] for c in X.columns]
srr2sam=dict(zip(smeta['Name'],smeta['Sample Name']))
X.columns=[srr2sam.get(c,c) for c in X.columns]
X=X.T  # samples x taxa
print('bracken',X.shape,'rowsum',X.sum(axis=1).describe().round(3).to_dict())
meta=pd.read_csv(I+'metadata.csv',index_col=0)
sub=meta.groupby('sample_id_tp1').first()
X=X.loc[[s for s in X.index if s in sub.index]]
info=sub.loc[X.index]
y=(info['study_ptorhc']=='MECFS').astype(int)
print('samples matched to BioMapAI metadata:',len(X),'patients',y.sum())

# plasma metabolome for the same subjects
met=pd.read_csv(B+'Metabolomics.csv',index_col=0)
cols=[c for c in met.columns if c in meta.index]
L=np.log2(met[cols].astype(float).T); L['subj']=meta.loc[cols,'sample_id_tp1'].values
La=L.groupby('subj').mean()
z=lambda s:(s-s.mean())/s.std()
sco=lambda n: z(sum(z(La[x]) for x in n))
T=pd.DataFrame({'IPA':z(La['indolepropionate']),'ILA':z(La['indolelactate']),
                'fld':sco(['indolepropionate','3-phenylpropionate (hydrocinnamate)','cinnamoylglycine']),
                'arLDH':sco(['indolelactate','phenyllactate (PLA)','3-(4-hydroxyphenyl)lactate'])})
sh=[s for s in X.index if s in T.index]
print('with paired plasma metabolome:',len(sh))

# --- canonical IPA / ILA producers
IPAprod=['Clostridium sporogenes','Clostridium botulinum','Clostridium cadaveris','Peptostreptococcus anaerobius',
         'Peptostreptococcus russellii','Peptostreptococcus stomatis','Clostridium paraputrificum',
         'Peptostreptococcus','Clostridium caloritolerans']
ILAprod=[c for c in X.columns if c.startswith('Bifidobacterium ') or c.startswith('Lactobacillus ') or c.startswith('Lacticaseibacillus ') or c.startswith('Limosilactobacillus ') or c.startswith('Ligilactobacillus ')]
ORAL=['Streptococcus salivarius','Streptococcus parasanguinis','Streptococcus mitis','Streptococcus oralis',
      'Streptococcus anginosus','Veillonella parvula','Veillonella atypica','Veillonella dispar',
      'Fusobacterium nucleatum','Fusobacterium periodonticum','Haemophilus parainfluenzae','Rothia mucilaginosa',
      'Rothia dentocariosa','Actinomyces oris','Actinomyces naeslundii','Gemella haemolysans','Gemella morbillorum',
      'Granulicatella adiacens','Atopobium parvulum','Solobacterium moorei','Porphyromonas gingivalis',
      'Prevotella melaninogenica','Prevotella intermedia','Campylobacter concisus','Leptotrichia buccalis',
      'Lactobacillus salivarius','Streptococcus sanguinis','Streptococcus gordonii','Eikenella corrodens',
      'Aggregatibacter aphrophilus','Neisseria subflava','Capnocytophaga gingivalis']
MUCIN=['Akkermansia muciniphila','Ruminococcus gnavus','Bacteroides fragilis','Barnesiella intestinihominis']
def present(l): return [x for x in l if x in X.columns]
groups={'IPA_producers':present(IPAprod),'Bifido_Lacto':ILAprod,'oral_taxa':present(ORAL),'mucin':present(MUCIN)}
for k,v in groups.items(): print(k,len(v))

S=pd.DataFrame(index=X.index)
for k,v in groups.items(): S[k]=X[v].sum(axis=1)
S['Akkermansia']=X.get('Akkermansia muciniphila',pd.Series(0,index=X.index))
S['R_gnavus']=X.get('Ruminococcus gnavus',pd.Series(0,index=X.index))
S['C_sporogenes']=X.get('Clostridium sporogenes',pd.Series(0,index=X.index))
S['fungi']=X[[c for c in X.columns if c.startswith(('Candida ','Saccharomyces ','Malassezia ','Aspergillus ','Debaryomyces '))]].sum(axis=1)
Sl=np.log10(S+1e-7)
rows=[]
for c in Sl.columns:
    a=Sl.loc[y==1,c]; b=Sl.loc[y==0,c]
    sd=np.sqrt(((len(a)-1)*a.var()+(len(b)-1)*b.var())/(len(a)+len(b)-2))
    d=(a.mean()-b.mean())/sd if sd>0 else np.nan
    p=stats.mannwhitneyu(a,b).pvalue
    r_ipa,p_ipa=stats.spearmanr(Sl.loc[sh,c],T.loc[sh,'IPA'])
    r_ila,p_ila=stats.spearmanr(Sl.loc[sh,c],T.loc[sh,'ILA'])
    r_fld,p_fld=stats.spearmanr(Sl.loc[sh,c],T.loc[sh,'fld'])
    r_arl,p_arl=stats.spearmanr(Sl.loc[sh,c],T.loc[sh,'arLDH'])
    rows.append(dict(group=c,med_pt=10**a.median(),med_hc=10**b.median(),cohen_d=round(d,2),p_grp=p,
                     rho_IPA=round(r_ipa,2),p_IPA=p_ipa,rho_fld=round(r_fld,2),p_fld=p_fld,
                     rho_ILA=round(r_ila,2),p_ILA=p_ila,rho_arLDH=round(r_arl,2),p_arLDH=p_arl))
D=pd.DataFrame(rows)
pd.set_option('display.width',250)
print('\n== guild-level (Bracken, n=%d; %d with plasma)'%(len(X),len(sh)))
print(D.to_string(index=False,float_format=lambda x:f'{x:.3g}'))
D.to_csv(OUT+'bracken_guilds.csv',index=False)
S.assign(pt=y).to_csv(OUT+'bracken_guild_abundances.csv')
