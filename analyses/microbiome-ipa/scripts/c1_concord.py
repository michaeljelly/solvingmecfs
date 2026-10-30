import pandas as pd, numpy as np
from scipy import stats
OUT='/mnt/project-files/mecfs-research/deep2/microbiome_function/'
g=pd.read_csv(OUT+'guo2023_species_maaslin.csv')
b=pd.read_csv(OUT+'biomap_species_assoc.csv')
# harmonise names (BioMapAI uses newer taxonomy: Phocaeicola<-Bacteroides etc.)
ren={'Bacteroides_vulgatus':'Phocaeicola_vulgatus','Bacteroides_dorei':'Phocaeicola_dorei',
     'Eubacterium_rectale':'Eubacterium_rectale','Bacteroides_plebeius':'Phocaeicola_plebeius',
     'Eubacterium_eligens':'Lachnospira_eligens','Eubacterium_hallii':'Anaerobutyricum_hallii',
     'Prevotella_copri':'Prevotella_copri_clade_A'}
g['sp2']=g['species'].replace(ren)
m=b.merge(g,left_on='feature',right_on='sp2',how='inner')
print('overlapping species tested in both:',len(m))
m['biomap_up']=np.sign(m['pt_beta'])
print(m[['feature','pt_beta','pt_q','coefficient','q-value (BH)','mecfs_up','biomap_up']].sort_values('pt_beta',ascending=False).to_string(index=False))
c=(m['biomap_up']==m['mecfs_up']).sum()
print(f'\nconcordant direction: {c}/{len(m)}  binomial p={stats.binomtest(c,len(m),0.5).pvalue:.4f}')
r,p=stats.spearmanr(m['pt_beta'],-m['coefficient'])
print(f'Spearman of effect sizes: rho={r:.2f} p={p:.4f}')
# restricted to those significant in BioMapAI
s=m[m.pt_q<0.1]
c2=(s['biomap_up']==s['mecfs_up']).sum()
print(f'among BioMapAI q<0.1 species ({len(s)}): concordant {c2}/{len(s)}')
m.to_csv(OUT+'biomap_vs_guo_species_concordance.csv',index=False)
