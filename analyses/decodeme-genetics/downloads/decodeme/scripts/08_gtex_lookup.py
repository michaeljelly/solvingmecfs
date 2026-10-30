# GTEx v8 portal API: significant single-tissue eQTLs and sQTLs for the lead variant of each DecodeME locus (plus top GWAS-1 SNP per locus)
import pandas as pd, requests, time
W='/tmp/work_decodeme/'
L=pd.read_csv(W+'results/lead_loci_all_gwas.csv')
leads=L[['locus','lead','gwas','log10p']].drop_duplicates('lead')
out=[]
for _,r in leads.iterrows():
    c,p,a,b=r.lead.split(':')
    vid=f'chr{c}_{p}_{a}_{b}_b38'
    for kind in ['singleTissueEqtl','singleTissueSqtl']:
        try:
            j=requests.get(f'https://gtexportal.org/api/v2/association/{kind}',params={'variantId':vid,'datasetId':'gtex_v8','itemsPerPage':500},timeout=60).json()
            for x in j.get('data',[]):
                out.append(dict(locus=r.locus,lead=r.lead,gwas=r.gwas,type=kind.replace('singleTissue',''),gene=x.get('geneSymbol'),tissue=x.get('tissueSiteDetailId'),p=x.get('pValue'),nes=x.get('nes')))
        except Exception as e: print(vid,kind,e)
        time.sleep(0.3)
R=pd.DataFrame(out); R.to_csv(W+'results/gtex_v8_lead_qtl_lookup.csv',index=False)
if len(R):
    s=R.groupby(['locus','lead','type','gene']).agg(n_tissues=('tissue','nunique'),min_p=('p','min'),nes_sign=('nes',lambda v: f"+{(v>0).sum()}/-{(v<0).sum()}")).reset_index()
    s=s.sort_values(['locus','type','n_tissues'],ascending=[True,True,False])
    s.to_csv(W+'results/gtex_v8_lead_qtl_summary.csv',index=False); print(s.to_string())
