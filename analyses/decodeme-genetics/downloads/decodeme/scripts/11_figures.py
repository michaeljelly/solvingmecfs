import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W='/tmp/work_decodeme/'
fig,axes=plt.subplots(2,1,figsize=(13,8))
for ax,g,t in zip(axes,['gwas_1','gwas_1_infectious_onset'],['GWAS-1: 15,579 cases vs 259,909 controls','GWAS-Infection: 9,738 infection-onset cases vs 259,909 controls']):
    d=pd.read_parquet(W+f'pq/{g}.parquet',columns=['CHROM','GENPOS','LOG10P'])
    d=d[(d.LOG10P>2)|(np.random.rand(len(d))<0.02)].sort_values(['CHROM','GENPOS'])
    off=d.groupby('CHROM').GENPOS.max().cumsum().shift(fill_value=0); d['x']=d.GENPOS+d.CHROM.map(off)
    for c,dd in d.groupby('CHROM'): ax.scatter(dd.x,dd.LOG10P,s=2,c='#1f4e79' if c%2 else '#7f9fbf',rasterized=True)
    ax.axhline(-np.log10(5e-8),c='r',lw=.7,ls='--'); ax.set_ylabel('-log10 p'); ax.set_title(t,fontsize=10)
    ax.set_xticks((d.groupby('CHROM').x.min()+d.groupby('CHROM').x.max())/2); ax.set_xticklabels(range(1,23),fontsize=7)
plt.tight_layout(); plt.savefig(W+'results/manhattan_gwas1_infection.png',dpi=130)
P=pd.read_csv(W+'results/magma_gwas_1_tissue_celltype.csv')
fig,axes=plt.subplots(1,3,figsize=(15,7))
for ax,p,t in zip(axes,['gtex54','monaco_immune','brain_sn'],['GTEx v8 tissues','Monaco immune cells (HPA)','Human brain single-nuclei clusters (HPA/Siletti)']):
    x=P[P.panel==p].sort_values('P').head(25).iloc[::-1]
    ax.barh(x.VARIABLE.str.slice(0,38),-np.log10(x.P),color=np.where(x.FDR<0.05,'#b03a2e','#999999'))
    ax.axvline(-np.log10(0.05/len(P[P.panel==p])),ls='--',c='k',lw=.7); ax.set_title(t,fontsize=10); ax.set_xlabel('-log10 p (MAGMA, conditioned on average expression)'); ax.tick_params(axis='y',labelsize=7)
plt.tight_layout(); plt.savefig(W+'results/magma_tissue_celltype.png',dpi=130)
