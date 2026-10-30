# Map MAGMA-significant genes and cis-MR protein hits to drugs (Open Targets Platform GraphQL, DGIdb GraphQL) and rank with direction of effect.
import pandas as pd, numpy as np, requests, json, time, sys
R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
OT='https://api.platform.opentargets.org/api/v4/graphql'; DG='https://dgidb.org/api/graphql'
mg=pd.read_csv('/mnt/project-files/mecfs-research/downloads/decodeme/results/magma_gwas_1_genes.csv')
mag=mg[mg.FDR<0.05][['SYMBOL','P','FDR']].rename(columns={'SYMBOL':'gene','P':'magma_p','FDR':'magma_fdr'})
mr=pd.read_csv(R+'mr_all_results.csv')
cis=mr[mr.category.isin(['protein_cis','protein_cis_rep'])&(mr.outcome=='DME_GWAS1')&(mr.p<0.01)][['label','category','b','se','p','fdr_in_category']].rename(columns={'label':'gene'})
# replication info for cis hits
rep=mr[mr.category.isin(['protein_cis','protein_cis_rep'])&(mr.outcome=='MVP_only')][['label','category','b','p']].rename(columns={'label':'gene','b':'b_mvp','p':'p_mvp'})
cis=cis.merge(rep,on=['gene','category'],how='left')
smr=pd.read_csv(R+'smr_direction_magma_genes.csv') if __import__('os').path.exists(R+'smr_direction_magma_genes.csv') else pd.DataFrame(columns=['gene','tissue','smr_b','p_smr'])
genes=sorted(set(mag.gene)|set(cis.gene.str.split(' ').str[0]))
def ot_q(q,v):
    for a in range(4):
        try: return requests.post(OT,json={'query':q,'variables':v},timeout=60).json()
        except Exception: time.sleep(3)
    return {}
Qs='query($q:String!){ search(queryString:$q, entityNames:["target"], page:{index:0,size:3}){ hits{ id name } } }'
Qt='''query($id:String!){ target(ensemblId:$id){ approvedSymbol approvedName tractability{label modality value}
  drugAndClinicalCandidates{ count rows{ maxClinicalStage drug{ name maximumClinicalStage mechanismsOfAction{ rows{ mechanismOfAction actionType } } } diseases{ disease{ name } } } } } }'''
rows=[]
for g in genes:
    h=ot_q(Qs,{'q':g}).get('data',{}).get('search',{}).get('hits',[])
    hit=[x for x in h if x['name']==g]
    if not hit: rows.append(dict(gene=g)); continue
    t=ot_q(Qt,{'id':hit[0]['id']}).get('data',{}).get('target') or {}
    tr=[x['label'] for x in (t.get('tractability') or []) if x['value'] and x['modality'] in ('SM','AB')]
    dc=(t.get('drugAndClinicalCandidates') or {}).get('rows',[])
    drugs=[]
    for d in dc:
        dn=d['drug']['name'] if d.get('drug') else None
        moa=';'.join({m['actionType'] or '' for m in ((d.get('drug') or {}).get('mechanismsOfAction') or {}).get('rows',[])})
        drugs.append(f"{dn} [{d['maxClinicalStage']}; {moa}]")
    rows.append(dict(gene=g,ensembl=hit[0]['id'],ot_n_drugs=len(dc),ot_drugs='; '.join(sorted(set(drugs))[:15]),tractable_SM_AB='; '.join(tr[:4])))
    time.sleep(0.2)
T=pd.DataFrame(rows)
# DGIdb
Qd='query($n:[String!]){ genes(names:$n){ nodes{ name interactions{ drug{ name approved } interactionScore interactionTypes{ type directionality } } } } }'
dg={}
for i in range(0,len(genes),25):
    try: r=requests.post(DG,json={'query':Qd,'variables':{'n':genes[i:i+25]}},timeout=90).json()
    except Exception: continue
    for n in r.get('data',{}).get('genes',{}).get('nodes',[]):
        ap=[f"{x['drug']['name']}({','.join(t['type'] for t in x['interactionTypes']) or '?'})" for x in n['interactions'] if x['drug']['approved']]
        dg[n['name']]=dict(dgidb_n_approved=len(ap),dgidb_approved='; '.join(sorted(set(ap))[:15]))
T=T.merge(pd.DataFrame.from_dict(dg,orient='index').reset_index().rename(columns={'index':'gene'}),on='gene',how='left')
T=T.merge(mag,on='gene',how='left')
c1=cis.assign(gene=cis.gene.str.split(' ').str[0]).sort_values('p').drop_duplicates('gene')
T=T.merge(c1[['gene','category','b','se','p','fdr_in_category','b_mvp','p_mvp']].rename(columns={'b':'cisMR_logOR_perSD','p':'cisMR_p','fdr_in_category':'cisMR_fdr','category':'pqtl_source'}),on='gene',how='left')
if len(smr):
    s1=smr.sort_values('p_smr').drop_duplicates('gene')[['gene','tissue','smr_b','p_smr']].rename(columns={'tissue':'smr_tissue'})
    T=T.merge(s1,on='gene',how='left')
def direction(r):
    b=r.get('cisMR_logOR_perSD')
    if pd.notna(b): return 'protein up -> risk up (inhibit)' if b>0 else 'protein up -> risk down (activate/avoid inhibitors)'
    b=r.get('smr_b')
    if pd.notna(b) and r.get('p_smr',1)<0.05: return 'expression up -> risk up (inhibit)' if b>0 else 'expression up -> risk down (activate/avoid inhibitors)'
    return 'unknown'
T['genetic_direction']=T.apply(direction,axis=1)
T['has_approved_drug']=(T.ot_drugs.fillna('').str.contains('APPROVAL'))|(T.dgidb_n_approved.fillna(0)>0)
T['evidence_score']=(-np.log10(T.magma_p.fillna(1)))+(-np.log10(T.cisMR_p.fillna(1)))+np.where(T.p_mvp.fillna(1)<0.05,2,0)
T=T.sort_values(['has_approved_drug','evidence_score'],ascending=False)
T.to_csv(R+'drug_target_mapping.csv',index=False)
print(T[['gene','magma_fdr','cisMR_p','cisMR_logOR_perSD','p_mvp','genetic_direction','ot_n_drugs','dgidb_n_approved','has_approved_drug']].head(40).to_string())
