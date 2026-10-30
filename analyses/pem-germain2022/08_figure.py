import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
C=pd.read_pickle('contrasts.pkl'); meta=pd.read_pickle('meta.pkl'); perm=pd.read_csv('permutation_global.csv')
BL,OR,TX,MU='#2a78d6','#eb6834','#0b0b0b','#52514e'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':MU,'axes.labelcolor':TX,'xtick.color':MU,'ytick.color':MU})
fig,ax=plt.subplots(2,2,figsize=(10,7.5),facecolor='#fcfcfb')
# A: permutation
a=ax[0,0]; p=perm.set_index('contrast'); lab=['Rest\nday 1','Rest\nday 2','Response\nday 1','Response\nday 2','Day 2 minus\nday 1','24 h\nrecovery']
keys=['LEVEL D1-PRE','LEVEL D2-PRE','R1 (D1 post-pre)','R2 (D2 post-pre)','R2-R1','REC (D2pre-D1pre)']
x=np.arange(len(keys)); a.bar(x,p.loc[keys,'obs_p05'],color=OR,width=0.6,label='observed')
a.errorbar(x,p.loc[keys,'null_median'],yerr=[np.zeros(len(x)),p.loc[keys,'null_95']-p.loc[keys,'null_median']],fmt='o',color=MU,capsize=4,label='chance (median, 95th pct)')
a.set_xticks(x,lab,fontsize=8); a.set_ylabel('metabolites with p<0.05 (of ~780)'); a.set_title('Patients differ at rest, not in their response',loc='left',color=TX,fontsize=10); a.legend(frameon=False,fontsize=8)
def traj(a,f,title,ylab):
    tp=['LEVEL D1-PRE','LEVEL D1-POST','LEVEL D2-PRE','LEVEL D2-POST']
    for g,c,n in [(0,BL,'Controls'),(1,OR,'ME/CFS')]:
        ids=meta.index[meta.case==g]; M=[C[t][f].reindex(ids) for t in tp]
        m=[v.mean() for v in M]; se=[v.std()/np.sqrt(v.notna().sum()) for v in M]
        a.errorbar([0,1,2.2,3.2],m,yerr=se,color=c,lw=2,marker='o',ms=6,capsize=3,label=n)
    a.set_xticks([0,1,2.2,3.2],['D1 pre','D1 post','D2 pre','D2 post']); a.set_title(title,loc='left',color=TX,fontsize=10); a.set_ylabel(ylab)
traj(ax[0,1],'lactate','Lactate: smaller rise in patients (tracks effort)','log2 scaled level'); ax[0,1].legend(frameon=False)
traj(ax[1,0],'glucose','Glucose: rises after exercise in controls only','log2 scaled level')
traj(ax[1,1],'cys-gly, oxidized','Oxidised cysteinylglycine: rises in patients on day 1','log2 scaled level')
fig.text(0.01,0.005,'Germain et al. 2022 (JCI Insight) data, 60 ME/CFS, 45 controls; mean +/- SE. Re-analysis: mecfs-research/pem.',color=MU,fontsize=7)
plt.tight_layout(rect=(0,0.02,1,1)); plt.savefig('pem_figure.png',dpi=150)
