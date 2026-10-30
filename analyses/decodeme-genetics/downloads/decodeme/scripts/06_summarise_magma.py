# Summarise MAGMA outputs: gene-level (BH-FDR), gene-set per collection (BH-FDR within collection), gene-property (tissue/cell), project candidate genes
import pandas as pd, numpy as np, sys, os
from statsmodels.stats.multitest import multipletests as mt
W='/tmp/work_decodeme/'; O=W+'magma_out/'; R=W+'results/'
G=sys.argv[1] if len(sys.argv)>1 else 'gwas_1'
loc=pd.read_csv(W+'magma/NCBI37.3.gene.loc',sep='\t',header=None,names=['GENE','c','s','e','str','SYMBOL'])
g=pd.read_csv(O+f'{G}.genes.out',sep=r'\s+').merge(loc[['GENE','SYMBOL']],on='GENE')
g['FDR']=mt(g.P,method='fdr_bh')[1]; g['bonf_sig']=g.P<0.05/len(g)
g.sort_values('P').to_csv(R+f'magma_{G}_genes.csv',index=False)
def rd(f):
    x=pd.read_csv(f,sep=r'\s+',comment='#'); 
    if 'FULL_NAME' in x.columns: x['VARIABLE']=x.FULL_NAME
    return x
rows=[]
for s in ['h.all','c2.cp.reactome','c2.cp.kegg_medicus','c5.go.bp','c5.go.cc','c5.go.mf','c8.all','c7.immunesigdb','project_sets']:
    f=O+f'{G}.gs.{s}.gsa.out'
    if not os.path.exists(f): continue
    x=rd(f); x['collection']=s; x['FDR']=mt(x.P,method='fdr_bh')[1]; x['bonf']=np.minimum(x.P*len(x),1); rows.append(x)
S=pd.concat(rows); S.sort_values('P').to_csv(R+f'magma_{G}_genesets.csv',index=False)
rows=[]
for c in ['gtex54','hpa_sc','monaco_immune','brain_sn']:
    f=O+f'{G}.gp.{c}.gsa.out'
    if not os.path.exists(f): continue
    x=rd(f); x['panel']=c; x['FDR']=mt(x.P,method='fdr_bh')[1]; x['bonf']=np.minimum(x.P*len(x),1); rows.append(x)
P=pd.concat(rows); P.sort_values(['panel','P']).to_csv(R+f'magma_{G}_tissue_celltype.csv',index=False)
pd.set_option('display.width',250)
print('genes',len(g),'bonf',g.bonf_sig.sum(),'FDR<.05',(g.FDR<.05).sum())
print(S.sort_values('P').groupby('collection').head(5)[['collection','VARIABLE','NGENES','BETA','P','FDR']].to_string())
print(P.sort_values('P').groupby('panel').head(8)[['panel','VARIABLE','BETA','P','FDR']].to_string())
# project candidate genes
sym=pd.read_csv(W+'gs/project_sets_symbols.tsv',sep='\t',header=None,index_col=0)[1]
C=[]
for k,v in sym.items():
    if k=='ALL_PROJECT_CANDIDATES': continue
    for s_ in v.split():
        r=g[g.SYMBOL==s_]
        C.append(dict(set=k,gene=s_,P=r.P.iloc[0] if len(r) else np.nan,Z=r.ZSTAT.iloc[0] if len(r) else np.nan,genome_FDR=r.FDR.iloc[0] if len(r) else np.nan))
C=pd.DataFrame(C); u=C.dropna().drop_duplicates('gene').copy(); u['FDR_within_candidates']=mt(u.P,method='fdr_bh')[1]
C=C.merge(u[['gene','FDR_within_candidates']],on='gene',how='left')
C.sort_values('P').to_csv(R+f'magma_{G}_project_candidate_genes.csv',index=False)
print(C.sort_values('P').head(15).to_string())
print('candidates tested',u.shape[0],'nominal p<0.05:',(u.P<0.05).sum(),'expected',round(0.05*len(u),1))
