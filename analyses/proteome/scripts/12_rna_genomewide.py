"""Genome-wide blood transcriptome case-control (4 cohorts, baseline) for the convergence analysis. Uses locally processed GEO data."""
import sys, gzip, re; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
G='/mnt/project-files/mecfs-research/downloads/geo/data/'
gi=pd.read_csv(W+'raw/Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','Symbol','dbXrefs','type_of_gene'])
e2s=gi.set_index(gi.GeneID.astype(str)).Symbol
ens=gi.assign(ens=gi.dbXrefs.str.extract(r'Ensembl:(ENSG\d+)')[0]).dropna(subset=['ens']).drop_duplicates('ens').set_index('ens').Symbol
def logcpm(c): return np.log2((c+0.5)/(c.sum(0)+1)*1e6)
def filt(c,lm,min_cpm=1,frac=0.25): return lm[((c/c.sum(0)*1e6)>=min_cpm).mean(1)>=frac]
def smeta(path):
    d={}
    for l in gzip.open(path,'rt'):
        if l.startswith('!Sample_'):
            p=l.rstrip('\n').split('\t'); k=p[0][8:]; v=[x.strip('"') for x in p[1:]]
            if k=='characteristics_ch1':
                kk=v[0].split(':')[0].strip(); d[kk]=[x.split(':',1)[1].strip() if ':' in x else x for x in v]
            else: d[k]=v
    return pd.DataFrame(d)
out={}
def run(name,c,m,cov,glob=False):
    lm=filt(c,logcpm(c)).T; des=m[['case']+list(cov)].astype(float)
    if glob:
        zz=(lm-lm.mean())/lm.std(); des['glob']=zz.mean(1).reindex(des.index)
    r=ols_many(lm.loc[des.index],des,'case'); r['n_case']=int(des.case.sum()); r['n_ctrl']=int((des.case==0).sum()); r.index.name='gene'
    r.to_csv(R+f'rna_{name}.csv'); out[name]=r
    print(name,des.case.value_counts().to_dict(),'genes',len(r),'q<0.05',(r.q<0.05).sum(),'p<0.05',(r.p<0.05).sum())
# GSE227375 PBMC baseline
m=pd.read_csv(G+'processed/GSE227375_meta.csv',index_col=0); m=m[m.time==0]; m['male']=(m.sex=='male').astype(int); m['plat']=(m.platform=='GPL21290').astype(int)
c=pd.read_parquet(G+'processed/GSE227375_counts_symbol.parquet'); ix=pd.Series(c.index.astype(str)); c.index=np.where(ix.str.isdigit(),e2s.reindex(ix).values,ix.values); c=c[c.index.notna()].groupby(level=0).sum()
run('GSE227375_PBMC',c[m.index],m,('male','plat'))
# GSE214283 scRNA pseudobulk all cells, D1
m=pd.read_csv(G+'processed/GSE214283_meta.csv').set_index('geo_accession'); m=m[m['collection day']=='D1']
m['case']=(m['disease state']=='case').astype(int); m['male']=(m.Sex=='M').astype(int); m['age']=m['age bin'].astype(float)
c=pd.read_parquet(G+'processed/GSE214283_pseudobulk_ALL.parquet'); m=m.loc[m.index.intersection(c.columns)]
run('GSE214283_PBMCsc',c[m.index],m,('male','age'))
# GSE98139 adolescent whole blood
c=pd.read_csv(G+'raw_small/GSE98139_All_normalized_counts.csv.gz',sep=';',index_col=0); c.index=ens.reindex(c.index.str.split('.').str[0]).values; c=c[c.index.notna()].groupby(level=0).sum().round()
m=smeta(G+'raw_small/GSE98139_series_matrix.txt.gz'); print(m.columns.tolist()[:15])
m['id']=m.title.str.extract(r'(\d+)$')[0]; m=m.set_index('id'); cond=[k for k in m.columns if m[k].astype(str).str.contains('CFS').any()][0]; m['case']=(m[cond]=='CFS').astype(int)
m=m.loc[m.index.intersection(c.columns)]
run('GSE98139_WB',c[m.index],m,())
# GSE293840 cfRNA
c=pd.read_csv(G+'raw_small/GSE293840_raw_counts_all.csv.gz',index_col=0); c.index=ens.reindex(c.index.str.split('.').str[0]).values; c=c[c.index.notna()].groupby(level=0).sum()
m=pd.read_csv(G+'processed/GSE293840_meta.csv',index_col=0); m=m.loc[m.index.intersection(c.columns)]
cov=['male','age']+[k for k in m.columns if k.startswith('batch_') or k.startswith('test_site')]
run('GSE293840_cfRNA',c[m.index],m,cov,glob=True)
# meta (Stouffer weighted by sqrt n) across the 4
genes=sorted(set().union(*[set(r.index) for r in out.values()]))
Zm=pd.DataFrame({k:(r.t.apply(lambda t: t)) for k,r in out.items()}).reindex(genes)
# convert t to z via p
Zz=pd.DataFrame({k:np.sign(r.t)*stats.norm.isf(r.p/2) for k,r in out.items() if len(r)}).reindex(genes)
w=pd.Series({k:np.sqrt(r.n_case.dropna().iloc[0]+r.n_ctrl.dropna().iloc[0]) for k,r in out.items() if len(r)})
num=(Zz*w).sum(1,min_count=2); den=np.sqrt((Zz.notna()*w**2).sum(1)); z=num/den
M=pd.DataFrame({'z':z,'k':Zz.notna().sum(1)}).dropna(); M['p']=2*stats.norm.sf(M.z.abs()); M['q']=bh(M.p.values)
M=M.join(Zz.add_prefix('z_')); M.to_csv(R+'rna_meta_genomewide.csv')
print('RNA meta genes',len(M),'q<0.05',(M.q<0.05).sum()); print(M.sort_values('p').head(20).round(4).to_string())
