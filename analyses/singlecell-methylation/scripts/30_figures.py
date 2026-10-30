import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'; F='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/figures'
BLUE,ORANGE,INK,MUTED='#2a78d6','#eb6834','#0b0b0b','#8a8984'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':MUTED,'axes.labelcolor':INK})
# Fig 1: scRNA pre-specified features, baseline effect sizes with 95% CI
O=pd.read_csv(f'{R}/sc_composition_scores_tests.csv'); O=O[O.contrast=='baseline_D1']
sel=['score:CD8_Tem:pct_mt','score:CD4_naive_CM:pct_mt','score:MAIT:pct_mt','score:NK_CD56dim:pct_mt','score:Mono_classical:pct_mt',
     'sub:KLRG1+CX3CR1+_of_CD8mem','sub:NKadaptive_KLRC2+FCER1G-_of_NK','sub:ABC_celltypist_of_B','sub:ABC_TBX21+orITGAX+_of_B','sub:Exhausted_PDCD1+TOX+_of_CD8mem',
     'score:pDC:s_ISG','score:Mono_classical:s_ISG','score:Mono_classical:s_hypoxia','sub:PlateletRNA+_of_Mono','clr:Platelet','sub:Th17like_CCR6+KLRB1+_of_CD4mem','clr:MAIT','sub:gdT_of_T','score:Mono_classical:s_adrenergic']
d=O.set_index('feature').loc[[s for s in sel if s in set(O.feature)]]
se=np.abs(d.d)/np.maximum(stats.norm.isf(d.p/2),1e-6)
fig,ax=plt.subplots(figsize=(6.2,5.6)); y=np.arange(len(d))[::-1]
ax.axvline(0,color=MUTED,lw=1)
ax.errorbar(d.d,y,xerr=1.96*se,fmt='o',color=BLUE,ms=5,lw=1.5,capsize=0)
for yi,(f,r) in zip(y,d.iterrows()): ax.text(1.45,yi,f'q={r.q:.2f}',va='center',fontsize=7,color=INK)
ax.set_yticks(y); ax.set_yticklabels([s.replace('score:','').replace('sub:','').replace('clr:','prop ') for s in d.index],fontsize=7)
ax.set_xlim(-1.6,1.7); ax.set_xlabel('ME/CFS minus control, SD units (sex, age bin adjusted), 95% CI')
ax.set_title('GSE214283 scRNA, baseline (30 ME/CFS, 28 controls)',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig(f'{F}/fig1_scRNA_prespecified_baseline.png',dpi=160)
# Fig 2: GSE93266 granulocyte contamination + clocks across datasets
S=pd.read_csv(f'{R}/meth_GSE93266_sample_scores.csv')
T=pd.read_csv(f'{R}/meth_clocks_composition_group_tests.csv')
fig,axs=plt.subplots(1,2,figsize=(9,4),gridspec_kw={'width_ratios':[1,1.6]})
ax=axs[0]
for i,(g,c) in enumerate([('HC',BLUE),('ME',ORANGE)]):
    v=S.loc[S.group==g,'rein_neutrophil']; ax.scatter(np.random.default_rng(i).normal(i,0.06,len(v)),v,s=14,color=c,alpha=.8)
    ax.plot([i-.2,i+.2],[v.median()]*2,color=INK,lw=2)
ax.set_xticks([0,1]); ax.set_xticklabels(['Controls (26)','ME/CFS (49)']); ax.set_ylabel('Estimated granulocyte fraction (Reinius ref.)')
ax.set_title('GSE93266 "PBMC": granulocyte DNA\nmostly in control samples',fontsize=9,loc='left')
ax=axs[1]; cl=['Horvath','HorvathSkinBlood','Hannum','PhenoAge','DunedinPACE','DNAmTL','DNAmGDF15','DNAmPAI1']
ds=[('GSE156792','T cells 61/48',BLUE),('GSE93266','PBMC 49/26',ORANGE),('GSE111183','PBMC 13/12','#1baf7a'),('GSE327255','saliva 54/21 (age adj)','#4a3aa7')]
for k,(g,lab,c) in enumerate(ds):
    t=T[(T.dataset==g)].set_index('score').reindex(cl); yy=np.arange(len(cl))[::-1]+ (k-1.5)*0.18
    se=np.abs(t.d_adj)/np.maximum(stats.norm.isf(t.p/2),1e-6)
    ax.errorbar(t.d_adj,yy,xerr=1.96*se,fmt='o',ms=4,color=c,lw=1,label=lab,capsize=0)
ax.axvline(0,color=MUTED,lw=1); ax.set_yticks(np.arange(len(cl))[::-1]); ax.set_yticklabels(cl); ax.set_xlabel('ME/CFS minus control, SD units, 95% CI')
ax.legend(fontsize=7,frameon=False,loc='lower right'); ax.set_title('Epigenetic clocks and DNAm surrogates',fontsize=9,loc='left')
fig.tight_layout(); fig.savefig(f'{F}/fig2_methylation_confound_clocks.png',dpi=160)
print('ok')
