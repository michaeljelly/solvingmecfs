import pandas as pd, numpy as np
B='/mnt/project-files/biomap/BioMapAI/codes/AI/input/'
def load(layer):
    d=pd.read_csv(B+layer+'.csv',index_col=0).T
    return d
meta=pd.read_csv(B+'metadata.csv',index_col=0)
score=pd.read_csv(B+'score.csv',index_col=0)
meta['subject']=meta['sample_id_tp1']
meta['y']=(meta.study_ptorhc=='MECFS').astype(int)
LAYERS=['immune','specie','kegg','metabolome','quest','omics']
