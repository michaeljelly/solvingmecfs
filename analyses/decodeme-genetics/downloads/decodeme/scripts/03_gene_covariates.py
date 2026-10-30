# Build MAGMA gene-property covariate files (Entrez IDs): GTEx v8 54 tissues, HPA single-cell types,
# HPA/Monaco immune cells (incl. Th17), HPA brain single-nuclei clusters. log2(x+1) per cell/tissue + Average column (FUMA style).
import pandas as pd, numpy as np, gzip
W='/tmp/work_decodeme/'
gi=pd.read_csv(W+'ref/Homo_sapiens.gene_info.gz',sep='\t',usecols=['GeneID','Symbol','dbXrefs','type_of_gene'])
gi=gi[gi.type_of_gene=='protein-coding']
ens={}
for g,x in zip(gi.GeneID,gi.dbXrefs):
    for t in str(x).split('|'):
        if t.startswith('Ensembl:'): ens.setdefault(t[8:],g)
loc=pd.read_csv(W+'magma/NCBI37.3.gene.loc',sep='\t',header=None,names=['id','chr','s','e','str','sym'])
keep=set(loc.id)
def write(mat,name):
    mat=mat.copy(); mat.index=[ens.get(e.split('.')[0]) for e in mat.index]
    mat=mat[pd.notna(mat.index)]; mat.index=mat.index.astype(int)
    mat=mat.groupby(level=0).mean(); mat=mat[mat.index.isin(keep)]
    mat=np.log2(mat+1); mat.columns=[c.replace(' ','_').replace('/','_').replace('(','').replace(')','').replace('-','_').replace(',','') for c in mat.columns]
    mat['Average']=mat.mean(axis=1)
    mat.index.name='GENE'; mat.to_csv(W+f'ref/covar_{name}.txt',sep='\t')
    print(name,mat.shape)
g=pd.read_csv(W+'ref/gtex_v8_median_tpm.gct.gz',sep='\t',skiprows=2).set_index('Name').drop(columns='Description')
write(g,'gtex54')
for f,col,val,name in [('rna_single_cell_type','Cell type','nCPM','hpa_sc'),
                       ('rna_immune_cell_monaco','Immune cell','nTPM','monaco_immune'),
                       ('rna_single_nuclei_cluster_type','Cluster type','nCPM','brain_sn')]:
    d=pd.read_csv(W+f'ref/{f}.tsv',sep='\t')
    m=d.pivot_table(index='Gene',columns=col,values=val,aggfunc='mean')
    if 'total PBMC' in m.columns: m=m.drop(columns='total PBMC')
    if m.isna().any().any(): m=m.dropna()  # HPA sc: 1031 genes lack some cell types; dropped
    write(m,name)
