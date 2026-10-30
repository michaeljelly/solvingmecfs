import pandas as pd, numpy as np
m=pd.read_csv('/tmp/work_geo/GSE227375_meta_raw.csv',index_col=0).sort_values('num'); C=pd.read_csv('/tmp/work_geo/GSE227375_idcorr.csv',index_col=0)
m['prefix']=m['description'].str[:2]; subj=[]; s=-1; prev=None
for i,r in m.iterrows():
    if r.num>5000: subj.append(None); continue
    if prev is None or r.prefix!=prev.prefix or r.time<=prev.time or r.num-prev.num>2: s+=1
    subj.append(f'S{s}'); prev=r
m['subject']=subj
d2s=dict(zip(m.description,m.subject))
m.loc[m.description=='HF2_5264','subject']=d2s['HF0_2956']; m.loc[m.description=='CF0_5266','subject']=d2s['CF1_2690']; m.loc[m.description=='HF1_5271','subject']=d2s['HF0_2836']
# validate
w=[];b=[]
for i in m.index:
    for j in m.index:
        if i<j: (w if m.subject[i]==m.subject[j] else b).append(C.loc[i,j])
print('within-subject r median',np.median(w),'n',len(w),' min',np.min(w),'| between median',np.median(b),'99th pct',np.percentile(b,99))
print(m.groupby(['case','subject']).size().groupby('case').describe())
dup=m.groupby(['subject','time']).size(); print('dup subject-time',(dup>1).sum())
m.to_csv('/tmp/work_geo/GSE227375_meta.csv')
