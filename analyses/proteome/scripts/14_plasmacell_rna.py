"""Self-contained test of a plasma-cell / immunoglobulin transcript score per blood RNA cohort (correlation-safe)."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
exec(open('/tmp/work_proteome/scripts/12_rna_genomewide.py').read().split('out={}')[0])
from common import *
core=['JCHAIN','MZB1','TNFRSF17','IGHG1','IGHG3','IGHA1','IGKC','TXNDC5','POU2AF1']
rows=[]
def score(name,c,m,cov,glob=False):
    lm=filt(c,logcpm(c)).T; z=(lm-lm.mean())/lm.std(); g=[x for x in core if x in z.columns]
    des=m[['case']+list(cov)].astype(float)
    if glob: des['glob']=z.mean(1).reindex(des.index)
    s=z[g].mean(1).reindex(des.index).to_frame('pc'); r=ols_many(s,des,'case').iloc[0]
    rows.append(dict(cohort=name,n_genes=len(g),d=r.beta,p=r.p,n_case=int(des.case.sum()),n_ctrl=int((des.case==0).sum())))
m=pd.read_csv(G+'processed/GSE227375_meta.csv',index_col=0); m=m[m.time==0]; m['male']=(m.sex=='male').astype(int); m['plat']=(m.platform=='GPL21290').astype(int)
c=pd.read_parquet(G+'processed/GSE227375_counts_symbol.parquet'); ix=pd.Series(c.index.astype(str)); c.index=np.where(ix.str.isdigit(),e2s.reindex(ix).values,ix.values); c=c[pd.notna(c.index)].groupby(level=0).sum()
score('GSE227375 PBMC',c[m.index],m,('male','plat'))
m=pd.read_csv(G+'processed/GSE214283_meta.csv').set_index('geo_accession'); m=m[m['collection day']=='D1']; m['case']=(m['disease state']=='case').astype(int); m['male']=(m.Sex=='M').astype(int); m['age']=m['age bin'].astype(float)
c=pd.read_parquet(G+'processed/GSE214283_pseudobulk_ALL.parquet'); m=m.loc[m.index.intersection(c.columns)]; score('GSE214283 PBMC scRNA',c[m.index],m,('male','age'))
c=pd.read_csv(G+'raw_small/GSE98139_All_normalized_counts.csv.gz',sep=';',index_col=0); c.index=ens.reindex(c.index.str.split('.').str[0]).values; c=c[c.index.notna()].groupby(level=0).sum().round()
m=smeta(G+'raw_small/GSE98139_series_matrix.txt.gz'); m['id']=m.title.str.extract(r'(\d+)$')[0]; m=m.set_index('id'); m['case']=(m['condition']=='CFS').astype(int); m=m.loc[m.index.intersection(c.columns)]
score('GSE98139 whole blood',c[m.index],m,())
c=pd.read_csv(G+'raw_small/GSE293840_raw_counts_all.csv.gz',index_col=0); c.index=ens.reindex(c.index.str.split('.').str[0]).values; c=c[c.index.notna()].groupby(level=0).sum()
m=pd.read_csv(G+'processed/GSE293840_meta.csv',index_col=0); m=m.loc[m.index.intersection(c.columns)]
score('GSE293840 cfRNA',c[m.index],m,['male','age']+[k for k in m.columns if k.startswith('batch_') or k.startswith('test_site')],glob=True)
T=pd.DataFrame(rows); T.to_csv(R+'plasmacell_rna_score.csv',index=False); print(T.round(3).to_string())
se=(T.d/stats.norm.isf(T.p/2)).abs(); print(meta_re(T.d.values,se.values))
