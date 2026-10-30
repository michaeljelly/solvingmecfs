import pandas as pd, numpy as np, os
B='/mnt/project-files/biomap/omics/Omics_Dataset/'
I='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'

def rd(p):
    d=pd.read_csv(p,index_col=0)
    return d

meta=rd(I+'metadata.csv')
met=rd(B+'Metabolomics.csv')
sp=rd(B+'Specie Abundance.csv')
ko=rd(B+'KEGG.csv')
print('meta',meta.shape,'met',met.shape,'sp',sp.shape,'ko',ko.shape)
for n,d in [('met',met),('sp',sp),('ko',ko)]:
    print(n,'overlap with meta idx:',len(set(d.columns)&set(meta.index)))
common=sorted(set(met.columns)&set(sp.columns)&set(meta.index))
print('paired met+sp+meta samples:',len(common))
subj=meta.loc[common,'sample_id_tp1']
print('subjects:',subj.nunique())
print(meta.loc[common].groupby('sample_id_tp1').first()['study_ptorhc'].value_counts())
