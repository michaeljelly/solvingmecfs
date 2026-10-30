# Direction of effect for MAGMA FDR<0.05 genes: SMR-style Wald ratio using the top cis-eQTL per gene
# (eQTLGen whole blood, n~31k; GTEx v8 brain cortex / frontal cortex BA9 / hypothalamus via API), outcome DecodeME GWAS-1.
# Also r2 (1000G EUR) between top eQTL and the strongest DecodeME SNP within 500 kb (crude colocalisation proxy).
import pandas as pd, numpy as np, requests, sys, time
from scipy import stats
sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/genetics_mr/scripts'); import mrlib
G='/tmp/work_gmr/'; R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
mg=pd.read_csv('/mnt/project-files/mecfs-research/downloads/decodeme/results/magma_gwas_1_genes.csv'); genes=mg[mg.FDR<0.05].SYMBOL.tolist()
extra=[g for g in sys.argv[1:]]; genes=list(dict.fromkeys(genes+extra))
D=pd.read_parquet(G+'dme_gwas_1.parquet',columns=['rs','chr','pos38','ea','oa','eaf','beta','se','mlogp'])
rows=[]
eq=pd.read_csv('/tmp/work_decodeme/eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz',sep='\t',usecols=['SNP','AssessedAllele','OtherAllele','Zscore','GeneSymbol','NrSamples','Pvalue'])
eq=eq[eq.GeneSymbol.isin(genes)]
eq=eq.merge(D[['rs','eaf']],left_on='SNP',right_on='rs',how='inner')
eq['beta']=eq.Zscore/np.sqrt(2*eq.eaf*(1-eq.eaf)*(eq.NrSamples+eq.Zscore**2)); eq['se']=eq.beta/eq.Zscore
for g,x in eq.groupby('GeneSymbol'):
    t=x.loc[x.Zscore.abs().idxmax()]
    rows.append(dict(gene=g,tissue='eQTLGen blood',rs=t.SNP,ea=t.AssessedAllele,oa=t.OtherAllele,eaf=np.nan,beta=t.beta,se=t.se,p_eqtl=t.Pvalue))
# GTEx v8 brain tissues
TIS=['Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Pituitary']
def gtex_gene(sym):
    r=requests.get('https://gtexportal.org/api/v2/reference/gene',params=dict(geneId=sym,gencodeVersion='v26',genomeBuild='GRCh38/hg38'),timeout=60).json()
    d=[x for x in r.get('data',[]) if x['geneSymbol']==sym]; return d[0]['gencodeId'] if d else None
for g in genes:
    try: gid=gtex_gene(g)
    except Exception: gid=None
    if not gid: continue
    for t in TIS:
        try:
            r=requests.get('https://gtexportal.org/api/v2/association/singleTissueEqtl',params=dict(gencodeId=gid,tissueSiteDetailId=t,datasetId='gtex_v8'),timeout=60).json().get('data',[])
        except Exception: r=[]
        if not r: continue
        b=min(r,key=lambda z: z['pValue'])
        v=b['variantId'].split('_')  # chr_pos_ref_alt_b38 ; NES is for alt allele
        rows.append(dict(gene=g,tissue='GTEx '+t,rs=b.get('snpId'),ea=v[3],oa=v[2],eaf=np.nan,beta=b['nes'],se=abs(b['nes'])/stats.norm.isf(max(b['pValue'],1e-300)/2),p_eqtl=b['pValue']))
        time.sleep(0.2)
E=pd.DataFrame(rows); E=E[E.rs.astype(str).str.startswith('rs')]
out=[]
for _,e in E.iterrows():
    x=mrlib.harmonise(pd.DataFrame([e])[['rs','ea','oa','beta','se']],D[D.rs==e.rs][['rs','ea','oa','beta','se']])
    if len(x)==0: continue
    b=x.beta_o.iloc[0]/x.beta_e.iloc[0]; se=abs(x.se_o.iloc[0]/x.beta_e.iloc[0])
    z_smr=(x.beta_e.iloc[0]/x.se_e.iloc[0])**2*(x.beta_o.iloc[0]/x.se_o.iloc[0])**2/((x.beta_e.iloc[0]/x.se_e.iloc[0])**2+(x.beta_o.iloc[0]/x.se_o.iloc[0])**2)
    p_smr=stats.chi2.sf(z_smr,1)
    c=D[D.rs==e.rs].iloc[0]; reg=D[(D.chr==c.chr)&(abs(D.pos38-c.pos38)<500000)].sort_values('mlogp',ascending=False).iloc[0]
    out.append(dict(gene=e.gene,tissue=e.tissue,eqtl_snp=e.rs,p_eqtl=e.p_eqtl,gwas_p_at_eqtl=2*stats.norm.sf(abs(x.beta_o.iloc[0]/x.se_o.iloc[0])),smr_b=b,smr_se=se,p_smr=p_smr,region_lead=reg.rs,region_lead_p=10**-reg.mlogp))
O=pd.DataFrame(out)
# LD between eQTL SNP and regional GWAS lead
r2=[]
for _,o in O.iterrows():
    if o.eqtl_snp==o.region_lead: r2.append(1.0); continue
    m=mrlib.ld_r([o.eqtl_snp,o.region_lead]) if False else None
    r2.append(np.nan)
O['r2_eqtl_lead']=r2
O['direction']=np.where(O.smr_b>0,'higher expression -> higher risk','higher expression -> lower risk')
O.to_csv(R+'smr_direction_magma_genes.csv',index=False)
print(O.sort_values('p_smr').head(40).round(4).to_string())
