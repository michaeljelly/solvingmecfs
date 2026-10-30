# Targeted CpG tests: DecodeME locus genes (+-500kb gene lists), project pathway genes, published single-CpG claims.
import pandas as pd, numpy as np, sys
sys.path.insert(0,'/mnt/project-files/mecfs-research/downloads/geo/scripts')
from common import GENESETS
from statsmodels.stats.multitest import multipletests
from scipy import stats
W='/tmp/work_scm/meth'; R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
ann450=pd.read_parquet(f'{W}/ann450.parquet').set_index('IlmnID'); annE=pd.read_parquet(f'{W}/annEPIC.parquet').set_index('IlmnID')
loc=pd.read_csv('/mnt/project-files/mecfs-research/genetics/locus_genes_500kb.csv')
gcol=[c for c in loc.columns if 'gene' in c.lower()][0]; bcol=[c for c in loc.columns if c.lower() in('band','locus')]
dec_genes=set(loc[gcol].astype(str))|{'RABGAP1L','DARS2','BTN2A2','TRIM38','FBXL4','SUDS3','OLFM4','CCPG1','CA10','ARFGEF2','CSE1L','ZNFX1','RC3H1','TNFSF4'}
path={k:v for k,v in GENESETS.items() if k!='DecodeME_named'}
path['EBV_host']=['CR2','CD21','EBI3','BHLHE40','IRF4','TNFAIP3','ZEB2','TBX21','ITGAX','FCRL5','CD58','CD48','SLAMF1','SLAMF6','SH2D1A','HLA-DQB1','HLA-DRB1']
path['coagulation']=['F3','SERPINE1','PLAUR','THBD','F2R','TFPI','VWF','GP1BA','ITGA2B','PF4','PPBP']
path['PTPRN2_saliva']=['PTPRN2']
def genes_of(ann):
    return ann.UCSC_RefGene_Name.fillna('').str.split(';').apply(lambda x:set(x)-{''})
MOD=sys.argv[1] if len(sys.argv)>1 else 'C'
rows=[]
for g in ['GSE156792','GSE93266','GSE111183','GSE327255']:
    E=pd.read_parquet(f'{W}/ewas_{g}.parquet'); E['t_B']=E[f't_{MOD}']; E['p_B']=E[f'p_{MOD}']; ann=(ann450 if g in('GSE156792','GSE93266') else annE).reindex(E.index)
    gs=genes_of(ann); grp=ann.UCSC_RefGene_Group.fillna('')
    for setname,glist in [('DecodeME_500kb',dec_genes)]+[(k,set(v)) for k,v in path.items()]:
        m=gs.apply(lambda s:len(s&glist)>0)
        sub=E[m].copy(); sub['gene']=gs[m].apply(lambda s:';'.join(sorted(s&glist))); sub['region']=grp[m]
        sub['q_in_set']=multipletests(sub.p_B,method='fdr_bh')[1] if len(sub) else []
        sub['dataset']=g; sub['set']=setname; sub['cpg']=sub.index
        rows.append(sub[['dataset','set','cpg','gene','region','dbeta','t_B','p_B','q_in_set','n']])
T=pd.concat(rows); T.to_csv(f'{R}/meth_targeted_cpgs_model{MOD}.csv',index=False)
summ=T.groupby(['set','dataset']).agg(n_cpg=('p_B','size'),n_p05=('p_B',lambda p:(p<0.05).sum()),n_q05=('q_in_set',lambda q:(q<0.05).sum()),min_q=('q_in_set','min')).reset_index()
summ['frac_p05']=summ.n_p05/summ.n_cpg
print(summ.to_string())
# cross-dataset: CpGs with q_in_set<0.1 in any dataset -> look at others
hit=T[T.q_in_set<0.1]
for _,r in hit.iterrows():
    o=T[(T.cpg==r.cpg)&(T.set==r.set)&(T.dataset!=r.dataset)]
    print(r.set,r.dataset,r.cpg,r.gene,r.region,'dbeta %.3f t %.2f q %.3f'%(r.dbeta,r.t_B,r.q_in_set),'| others:','; '.join(f'{x.dataset} t={x.t_B:.2f} p={x.p_B:.3f}' for x in o.itertuples()))
print(T[T.cpg=='cg19803194'].to_string())
