# Sample-level mitochondrial read fraction and other QC by group (is the pct_mt signal global or cell-type specific?)
import pandas as pd, numpy as np, statsmodels.formula.api as smf
R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
CE=pd.read_parquet('/tmp/work_scm/cells_all.parquet',columns=['gsm','type','pct_mt','pct_ribo','lib','ngenes'])
meta=pd.read_csv(f'{R}/sc_sample_meta.csv').set_index('geo_accession'); meta['case']=(meta['disease state']=='case').astype(int); meta['age']=meta['age bin'].astype(float)
S=CE.groupby('gsm').agg(pct_mt=('pct_mt','mean'),pct_mt_med=('pct_mt','median'),ribo=('pct_ribo','mean'),lib=('lib','median'),ngenes=('ngenes','median'),ncell=('lib','size')).join(meta)
S['hi_mt']=CE.assign(h=CE.pct_mt>0.1).groupby('gsm').h.mean()
rows=[]
for v in ['pct_mt','pct_mt_med','hi_mt','ribo','lib','ngenes','ncell']:
    for day in ['D1','D2']:
        d=S[S.day==day].copy(); d['z']=(d[v]-d[v].mean())/d[v].std(); m=smf.ols('z~case+C(Sex)+age',d).fit()
        rows.append(dict(var=v,day=day,mean_case=d.loc[d.case==1,v].mean(),mean_ctrl=d.loc[d.case==0,v].mean(),d=m.params['case'],p=m.pvalues['case']))
QT=pd.DataFrame(rows); print(QT.round(4).to_string())
d=S[S.day=='D1'].copy()
print('adj log lib p=%.4f'%smf.ols('pct_mt~case+C(Sex)+age+np.log(lib)',d).fit().pvalues['case'])
m=smf.ols('pct_mt~case+C(Sex)+age+C(Q("chromium batch"))',d).fit(); print('adj batch p=%.4f'%m.pvalues['case'])
w=S.pivot_table(index='person',columns='day',values='pct_mt'); print('r(D1,D2) per person=%.2f'%w.corr().iloc[0,1])
# per type pct_mt case effect using both days (mixed-free: average D1,D2 per person)
T=CE.groupby(['gsm','type']).pct_mt.mean().unstack(); T=T.join(meta[['person','case','Sex','age']])
P=T.groupby('person').mean(numeric_only=True); P['Sex']=T.groupby('person').Sex.first()
out=[]
for ct in [c for c in CE.type.unique() if c in P.columns and P[c].notna().sum()>50]:
    d=P.dropna(subset=[ct]).copy(); d['z']=(d[ct]-d[ct].mean())/d[ct].std(); m=smf.ols('z~case+C(Sex)+age',d).fit()
    out.append(dict(celltype=ct,d=m.params['case'],p=m.pvalues['case'],n=len(d)))
out=pd.DataFrame(out).sort_values('p'); print(out.round(4).to_string())
S.to_csv(f'{R}/sc_sample_quality.csv'); QT.to_csv(f'{R}/sc_sample_quality_tests.csv',index=False); out.to_csv(f'{R}/sc_pctmt_by_celltype_personmean.csv',index=False)
