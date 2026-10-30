import pandas as pd, numpy as np, gzip, re, sys
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
R='/tmp/work_geo/raw/'
gi=pd.read_csv(R+'Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','Symbol'],dtype={'GeneID':str}).set_index('GeneID')['Symbol']
c=pd.read_csv(R+'GSE227375_ncbi_counts.tsv.gz',sep='\t',dtype={'GeneID':str}).set_index('GeneID')
c.index=gi.reindex(c.index).fillna(pd.Series(c.index,index=c.index)).values; c=c.groupby(level=0).sum()
rows=[]
import glob
for f in glob.glob(R+'GSE227375-GPL*_series_matrix.txt.gz'):
    d={}
    for line in gzip.open(f,'rt'):
        if line.startswith('!Sample_'):
            k,*v=line.rstrip('\n').split('\t'); v=[x.strip('"') for x in v]
            if k=='!Sample_characteristics_ch1': d[v[0].split(':')[0]]=[x.split(':',1)[1].strip() for x in v]
            elif k in('!Sample_geo_accession','!Sample_description','!Sample_title'): d[k[8:]]=v
    df=pd.DataFrame(d); df['platform']=f.split('-')[1].split('_')[0]; rows.append(df)
m=pd.concat(rows).set_index('geo_accession')
m['case']=m['disease state'].str.startswith('ME').astype(int); m['sex']=m['Sex']
m['time']=m['description'].str[2]; m['num']=m['description'].str.split('_').str[1].astype(int)
print(m.groupby(['case','sex','time']).size())
print('counts samples',c.shape, 'in meta', len(set(c.columns)&set(m.index)))
m=m.loc[m.index.intersection(c.columns)]; print(m.groupby(['case','time']).size())
lc=logcpm(c[m.index])
idg=['GSTM1','ERAP2','HLA-DRB4','HLA-DRB5','HLA-DQA2','HLA-DQB2','UGT2B17','CHI3L2','PSPHP1','LILRA3','NLRP2','GSTT2B','THNSL2','ZNF890P','KANSL1-AS1','LRRC37A4P','CCL3L3','CCL4L2','SIGLEC14','APOBEC3B','XIST','RPS4Y1','HLA-DRB1','HLA-DQA1','HLA-DQB1','MICA','GSTT1','CYP2E1','PILRB','FAM118A','ZNF37BP','POLR2J2','WBSCR27','ACCS','IFITM3']
idg=[g for g in idg if g in lc.index and lc.loc[g].std()>0.5]
print('id genes',len(idg))
z=zscore_rows(lc.loc[idg]); C=np.corrcoef(z.T.values); np.fill_diagonal(C,-1)
C=pd.DataFrame(C,index=lc.columns,columns=lc.columns)
# consecutive-ID hypothesis check
m=m.sort_values('num'); desc=m['description']
top=C.idxmax(1); m['top_match']=desc.reindex(top.values).values; m['top_r']=C.max(1).reindex(m.index).values
print(m[['description','top_match','top_r']].to_string())
m.to_csv('/tmp/work_geo/GSE227375_meta_raw.csv'); C.to_csv('/tmp/work_geo/GSE227375_idcorr.csv')
lc.to_parquet('/tmp/work_geo/GSE227375_logcpm.parquet'); c.to_parquet('/tmp/work_geo/GSE227375_counts.parquet')
