# Assemble instrument tables (long format) from GWAS Catalog top associations, UKB-PPP and deCODE supplementary pQTL tables.
import pandas as pd, numpy as np, re, openpyxl
from scipy import stats
D='/tmp/work_gmr/dl/'; O='/tmp/work_gmr/'
a=pd.read_parquet(D+'gwascat_assoc.parquet'); s=pd.read_parquet(D+'gwascat_studies.parquet')
def nsize(t):
    return sum(int(x.replace(',','')) for x in re.findall(r'(\d[\d,]*) (?:[A-Za-z]|\()',str(t)) if x.replace(',','').isdigit()) or np.nan
def parse_cat(x, binary=False):
    x=x.copy()
    x['rs']=x['STRONGEST SNP-RISK ALLELE'].str.split('-').str[0].str.strip()
    x['ea']=x['STRONGEST SNP-RISK ALLELE'].str.split('-').str[-1].str.strip().str.upper()
    x=x[x.rs.str.match(r'^rs\d+$')&x.ea.str.match(r'^[ACGT]+$')]
    x['p']=pd.to_numeric(x['P-VALUE'],errors='coerce')
    b=pd.to_numeric(x['OR or BETA'],errors='coerce')
    ci=x['95% CI (TEXT)'].fillna('')
    lo_hi=ci.str.extract(r'\[\s*(-?[\d.]+(?:[eE]-?\d+)?)\s*-\s*(-?[\d.]+(?:[eE]-?\d+)?)\s*\]')
    lo=pd.to_numeric(lo_hi[0],errors='coerce'); hi=pd.to_numeric(lo_hi[1],errors='coerce')
    dec=ci.str.contains('decrease',case=False)
    isOR=binary | (~ci.str.contains('increase|decrease',case=False) & (b>0) & lo.notna())
    beta=np.where(isOR,np.log(b),np.where(dec,-b,b))
    se=np.where(isOR,(np.log(hi.abs().clip(lower=1e-9))-np.log(lo.abs().clip(lower=1e-9))).__abs__()/3.92,(hi-lo).abs()/3.92)
    # fallback: se from p when CI missing
    z=stats.norm.isf(x.p.clip(lower=1e-300)/2)
    se=np.where(np.isfinite(se)&(se>0),se,np.abs(beta)/z)
    x['beta']=beta; x['se']=se; x['eaf']=pd.to_numeric(x['RISK ALLELE FREQUENCY'],errors='coerce')
    x['n']=x['INITIAL SAMPLE SIZE'].map(nsize)
    x=x[np.isfinite(x.beta)&np.isfinite(x.se)&(x.se>0)]
    return x
rows=[]
def add(df,eid,cat,src,label):
    d=df[['rs','ea','eaf','beta','se','p','n']].copy(); d['oa']=np.nan; d['exposure']=eid; d['category']=cat; d['source']=src; d['label']=label; rows.append(d)
# metabolites
ch=a[a.PUBMEDID=='36635386']; ch=ch[ch['STUDY ACCESSION'].str.slice(4).astype(int).between(90199621,90201020)]
ch=ch[~ch['DISEASE/TRAIT'].str.contains(' ratio',case=False)]
pc=parse_cat(ch)
for acc,g in pc.groupby('STUDY ACCESSION'): add(g,'Chen:'+acc,'metabolite','Chen2023 (CLSA, n~8k, catalog lead SNPs)',g['DISEASE/TRAIT'].iloc[0].replace(' levels',''))
su=parse_cat(a[a.PUBMEDID=='36357675'])
for acc,g in su.groupby('STUDY ACCESSION'): add(g,'Sur:'+acc,'metabolite_rep','Surendran2022 (INTERVAL/EPIC-Norfolk, n~14k, catalog)',g['DISEASE/TRAIT'].iloc[0].replace(' levels',''))
orr=parse_cat(a[a.PUBMEDID=='32929287'])
for acc,g in orr.groupby('STUDY ACCESSION'): add(g,'Orru:'+acc,'immune_cell','Orru2020 (SardiNIA, n~3.7k, catalog)',g['DISEASE/TRAIT'].iloc[0])
rf={'GCST002783':('BMI (Locke 2015, non-UKB)','BMI'),'GCST006624':('SBP (Evangelou 2018)','SBP'),'GCST006630':('DBP (Evangelou 2018)','DBP'),
    'GCST011369':('Ferritin (Bell 2021, non-UKB)','Ferritin'),'GCST011367':('Serum iron (Bell 2021)','Iron'),'GCST011366':('Transferrin saturation (Bell 2021)','TSAT'),'GCST011368':('TIBC (Bell 2021)','TIBC'),
    'GCST90131901':('Insomnia (Watanabe 2022, log-odds)','Insomnia'),
    # added at coordinator request (FinnGen POSTVIRFAT lab pattern: low reticulocytes, low sTfR, higher TSAT; IgG subclass deficiency; adrenal insufficiency)
    'GCST90002405':('Reticulocyte count (Vuckovic 2020, UKB)','Retic_count'),'GCST90002406':('Reticulocyte % of red cells (Vuckovic 2020)','Retic_pct'),
    'GCST90002387':('Immature reticulocyte fraction (Vuckovic 2020)','IRF'),'GCST90002385':('High light scatter reticulocyte count (Vuckovic 2020)','HLSR_count'),
    'GCST90002310':('Haemoglobin (Chen MH 2020, EUR)','Hb'),'GCST90002304':('Haematocrit (Chen MH 2020, EUR)','Hct'),
    'GCST90434833':('Soluble transferrin receptor TFRC (Dhindsa 2023, UKB-PPP)','sTfR'),'GCST90860264':('Circulating EPO (Koprulu 2026)','EPO'),
    'GCST90429011':('Morning plasma cortisol (Chan 2024, CORNET)','Cortisol_am'),
    'GCST008568':('Serum IgA (Jonsson 2017, Iceland)','IgA'),'GCST008576':('Serum IgG (Jonsson 2017)','IgG'),'GCST008575':('Serum IgM (Jonsson 2017)','IgM')}
for acc,(lab,eid) in rf.items():
    g=parse_cat(a[a['STUDY ACCESSION']==acc]); add(g,eid,'risk_factor','GWAS Catalog '+acc,lab)
I=pd.concat(rows); I=I[I.p<5e-8]
# UKB-PPP cis pQTLs (discovery, ST9)
wb=openpyxl.load_workbook(D+'sun2023_supp.xlsx',read_only=True); ws=wb['ST9']
r=list(ws.iter_rows(min_row=5,values_only=True)); h=r[0]; u=pd.DataFrame(r[1:],columns=h)
u=u[u['cis/trans']=='cis'].copy()
u['ea']=u['Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)'].str.split(':').str[3]; u['oa']=u['Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)'].str.split(':').str[2]
u=u.rename(columns={'rsID':'rs','A1FREQ (discovery)':'eaf','BETA (discovery, wrt. A1)':'beta','SE (discovery)':'se'})
u['p']=10**-u['log10(p) (discovery)'].astype(float); u['n']=34557
u['gene']=u['Assay Target']; u=u.sort_values('p').drop_duplicates('UKBPPP ProteinID')
U=pd.DataFrame({'rs':u.rs,'ea':u.ea,'oa':u.oa,'eaf':u.eaf.astype(float),'beta':u.beta.astype(float),'se':u.se.astype(float),'p':u.p,'n':u.n,'exposure':'PPP:'+u['UKBPPP ProteinID'],'category':'protein_cis','source':'UKB-PPP Sun2023 ST9 discovery','label':u.gene})
# deCODE cis (sentinel, rank 1 per protein)
wb=openpyxl.load_workbook(D+'ferk2021_supp.xlsx',read_only=True); ws=wb['ST02']
r=list(ws.iter_rows(min_row=3,values_only=True)); h=[str(x).replace('\n',' ').strip() for x in r[0]]; f=pd.DataFrame(r[1:],columns=h)
f=f[f['cis/ trans']=='cis'].copy(); f['mlp']=pd.to_numeric(f['-Log10(P) (adj.)'],errors='coerce'); f['b']=pd.to_numeric(f['beta (adj.)'],errors='coerce')
f=f.sort_values('mlp',ascending=False).drop_duplicates('SeqId')
from scipy.stats import norm
z=norm.isf(np.clip(10**-f.mlp/2,1e-300,1)); z=np.where(f.mlp>300,np.sqrt(2*f.mlp*np.log(10)),z)
F=pd.DataFrame({'rs':f.variant,'ea':f.Amin,'oa':f.Amaj,'eaf':pd.to_numeric(f['MAF (%)'],errors='coerce')/100,'beta':f.b,'se':np.abs(f.b)/z,'p':10**-f.mlp,'n':35559,'exposure':'deCODE:'+f.SeqId.astype(str),'category':'protein_cis_rep','source':'deCODE Ferkingstad2021 ST2','label':f['gene  (prot.)'] if 'gene  (prot.)' in f else f.iloc[:,2]})
F=F[F.rs.astype(str).str.startswith('rs')&(F.p<5e-8)&(F.se>0)]
ALL=pd.concat([I,U,F],ignore_index=True)
ALL.to_parquet(O+'instruments_catalog_pqtl.parquet',index=False)
print(ALL.groupby('category').agg(n_exp=('exposure','nunique'),n_snp=('rs','size')))
