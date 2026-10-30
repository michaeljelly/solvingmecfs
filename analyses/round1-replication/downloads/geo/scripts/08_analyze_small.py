import pandas as pd, numpy as np, sys, gzip
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
R='/tmp/work_geo/raw/'; O='/tmp/work_geo/results/'
def smeta(f):
    d={}
    for line in gzip.open(f,'rt'):
        if line.startswith('!Sample_'):
            k,*v=line.rstrip('\n').split('\t'); v=[x.strip('"') for x in v]
            if k=='!Sample_characteristics_ch1': d[v[0].split(':')[0].strip()]=[x.split(':',1)[1].strip() for x in v]
            elif k in('!Sample_geo_accession','!Sample_title','!Sample_source_name_ch1'): d[k[8:]]=v
    return pd.DataFrame(d)
gi=pd.read_csv(R+'Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','Symbol','dbXrefs'],dtype={'GeneID':str})
id2s=gi.set_index('GeneID')['Symbol']; gi['ens']=gi.dbXrefs.str.extract(r'Ensembl:(ENSG\d+)'); e2s=gi.dropna(subset=['ens']).drop_duplicates('ens').set_index('ens')['Symbol']
ALLG=[];ALLS=[];ALLP=[]
# ---- GSE98139 adolescent whole blood (Nguyen 2017), normalized counts
c=pd.read_csv(R+'GSE98139_All_normalized_counts.csv.gz',sep=';',index_col=0)
c=c[c.index.isin(e2s.index)]; c.index=e2s.reindex(c.index).values; c=c.groupby(level=0).sum()
m=smeta(R+'GSE98139_series_matrix.txt.gz'); m['id']=m.title.str.extract(r'(\d+)$')[0]; m=m.set_index('id'); m['case']=(m.condition=='CFS').astype(int)
m=m.loc[m.index.intersection(c.columns)]; c=c[m.index].round()
lm=filter_expr(logcpm(c),c,1,0.25); g,s,_,_=run_case_control('GSE98139_adolescent_WB',lm,m); ALLG.append(g);ALLS.append(s); print('GSE98139',m.case.value_counts().to_dict()); print(s.round(3)[['d','p','q']].to_string())
# ---- GSE128078 Bouquet 2019 whole blood, days 1,2,3,7 around CPET
c=pd.read_csv(R+'GSE128078_ncbi_counts.tsv.gz',sep='\t',dtype={'GeneID':str}).set_index('GeneID'); c.index=id2s.reindex(c.index).fillna('NA').values; c=c.groupby(level=0).sum()
m=smeta(R+'GSE128078_series_matrix.txt.gz').set_index('geo_accession'); m['case']=(m['disease state']=='ME/CFS').astype(int); m['day']=m['timepoint (day)'].astype(int); m['subject']=m['individual identifier']
m=m.loc[m.index.intersection(c.columns)]; c=c[m.index]
lm=filter_expr(logcpm(c),c,1,0.25)
g,s,_,_=run_case_control('GSE128078_WB',lm,m[m.day==1],extra_label='_day1'); ALLG.append(g);ALLS.append(s)
print('GSE128078 day1',m[m.day==1].case.value_counts().to_dict()); print(s.round(3)[['d','p','q']].to_string())
for d2 in [2,3]:
    r,w=run_paired('GSE128078_WB',lm,m,'subject','day',1,d2,label=f'day1->day{d2}'); ALLP.append(r)
    print(f'change day1->day{d2} sets p<0.1'); print(r[r.is_set&(r.p<0.1)].round(3)[['beta','d','p','q','n_case','n_ctrl']].to_string())
# ---- GSE16059 twins (Byrnes 2009), GPL570 RMA, paired within twin pair
ann=pd.read_csv(R+'GPL570.txt',sep='\t',skiprows=18,usecols=['ID','Gene Symbol'],dtype=str).dropna(); ann=ann[~ann['Gene Symbol'].str.contains('///')].set_index('ID')['Gene Symbol']
x=pd.read_csv(R+'GSE16059_series_matrix.txt.gz',sep='\t',comment='!',index_col=0)
x=x[x.index.isin(ann.index)]; x['sym']=ann.reindex(x.index).values
x['mean']=x.drop(columns='sym').mean(1); x=x.sort_values('mean',ascending=False).drop_duplicates('sym').set_index('sym').drop(columns='mean')  # most-expressed probe per gene
m=smeta(R+'GSE16059_series_matrix.txt.gz').set_index('geo_accession'); m['case']=(m['diagnonsis']!='unaffected').astype(int); m['pair']=m['twin pair']
m['pairdx']=m.groupby('pair')['diagnonsis'].transform(lambda s:[v for v in s if v!='unaffected'][0])
lm=x[m.index]; lm=lm[lm.mean(1)>np.percentile(lm.mean(1),25)]
for lab,mm in [('all_44pairs',m),('CFS_only_32pairs',m[m.pairdx=='CFS'])]:
    r,w=run_paired('GSE16059_twins_'+lab,lm,mm.assign(case2=mm.case),'pair','case',0,1)  # change = affected - unaffected
    # within 'group' here are all pairs; take within results from group of unaffected rows (case==0 at t0)
    ww=w[w.group=='control'].set_index('feature').rename(columns={'mean_change':'beta'}); ww['contrast']='affected_minus_cotwin'; ww['dataset']='GSE16059_twins_'+lab
    ALLP.append(ww.assign(is_set=ww.index.isin(GENESETS.keys())))
    print('twins',lab,len(mm)//2); print(ww[ww.index.isin(GENESETS.keys())].round(3)[['beta','d','p','q','n']].to_string())
    print(ww[(~ww.index.isin(GENESETS.keys()))&(ww.p<0.01)].round(3)[['beta','p','q']].to_string())
pd.concat(ALLG).to_csv(O+'small_genes_case_control.csv'); pd.concat(ALLS).to_csv(O+'small_sets_case_control.csv'); pd.concat(ALLP).to_csv(O+'small_paired.csv')
