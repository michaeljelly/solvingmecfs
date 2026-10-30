import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='/tmp/work_geo/results/'
T=pd.read_csv(O+'SUMMARY_sets_baseline_by_dataset.csv',index_col=0,dtype=str,keep_default_na=False)
v=T.map(lambda s: np.nan if s in('na','nd') else float(s.replace('*','')))
st=T.map(lambda s: '' if s in('na','nd') else ('**' if s.endswith('**') else ('*' if s.endswith('*') else '')))
cols=[c.replace(' (','\n(') for c in v.columns]
fig,ax=plt.subplots(figsize=(13,7)); im=ax.imshow(v.values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
ax.set_xticks(range(v.shape[1])); ax.set_xticklabels(cols,rotation=45,ha='right',fontsize=7); ax.set_yticks(range(v.shape[0])); ax.set_yticklabels(v.index,fontsize=8)
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        if not np.isnan(v.values[i,j]): ax.text(j,i,f'{v.values[i,j]:+.2f}{st.values[i,j]}',ha='center',va='center',fontsize=6)
plt.colorbar(im,label='effect size (SD units), ME/CFS minus control'); ax.set_title('Gene-set scores, ME/CFS vs controls at baseline (* p<0.05, ** q<0.05 within dataset)',fontsize=9)
plt.tight_layout(); plt.savefig(O+'fig_sets_by_dataset.png',dpi=130)
