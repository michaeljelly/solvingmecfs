import pandas as pd, numpy as np, sys, gzip
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
R='/tmp/work_geo/raw/'; O='/tmp/work_geo/results/'
gi=pd.read_csv(R+'Homo_sapiens.gene_info.gz',sep='\t',usecols=['Symbol','dbXrefs'])
gi['ens']=gi.dbXrefs.str.extract(r'Ensembl:(ENSG\d+)'); e2s=gi.dropna(subset=['ens']).drop_duplicates('ens').set_index('ens')['Symbol']
c=pd.read_csv(R+'GSE293840_raw_counts_all.csv.gz',index_col=0)
c.index=c.index.str.split('.').str[0]; c=c[c.index.isin(e2s.index)]; c.index=e2s.reindex(c.index).values; c=c.groupby(level=0).sum()
d={}
for line in gzip.open(R+'GSE293840_series_matrix.txt.gz','rt'):
    if line.startswith('!Sample_characteristics_ch1') or line.startswith('!Sample_description') or line.startswith('!Sample_geo_accession'):
        k,*v=line.rstrip('\n').split('\t'); v=[x.strip('"') for x in v]
        if k=='!Sample_geo_accession': d['gsm']=v
        elif 'Column name' in v[0]: d['col']=[x.split(':',1)[1].strip() for x in v]
        elif k=='!Sample_characteristics_ch1': d[v[0].split(':')[0]]=[x.split(':',1)[1].strip() for x in v]
m=pd.DataFrame(d).set_index('col'); m['case']=(m.phenotype=='case').astype(int); m['male']=(m.Sex=='male').astype(int); m['age']=pd.to_numeric(m['age bin'],errors='coerce'); m['age']=m['age'].fillna(m['age'].median())
m['site']=m['test site'].astype('category').cat.codes
m=pd.get_dummies(m,columns=['batch','test site'],drop_first=True,dtype=int); m.columns=[x.replace(' ','_').replace('.','_').replace(',','') for x in m.columns]
cov=['male','age']+[x for x in m.columns if x.startswith('batch_') or x.startswith('test_site_')]
m=m.loc[m.index.intersection(c.columns)]; c=c[m.index]
print(m.case.value_counts().to_dict(), c.shape, 'median lib',c.sum(0).median())
lm=filter_expr(logcpm(c),c,min_cpm=1,min_frac=0.25)
g,s,sc,miss=run_case_control('GSE293840_cfRNA',lm,m,covars=cov)
g.to_csv(O+'GSE293840_genes_case_control.csv'); s.to_csv(O+'GSE293840_sets_case_control.csv')
print('missing',miss); print(s.round(3)[['beta','d','p','q','n_genes']].to_string()); print(g[g.p<0.05].sort_values('p').round(3)[['beta','d','p','q']].to_string())
m.to_csv('/tmp/work_geo/GSE293840_meta.csv')
m['glob']=zscore_rows(lm).mean(0).reindex(m.index); print('global shift d:',round(ols_test(m.glob,m,'case',cov)['d'],2))
g2,s2,_,_=run_case_control('GSE293840_cfRNA_globadj',lm,m,covars=cov+['glob'],extra_label='_global_adjusted')
g2.to_csv(O+'GSE293840_genes_case_control_globadj.csv'); s2.to_csv(O+'GSE293840_sets_case_control_globadj.csv')
print(s2.round(3)[['beta','d','p','q']].to_string()); print(g2[g2.p<0.01].sort_values('p').round(3)[['beta','d','p','q']].to_string())
