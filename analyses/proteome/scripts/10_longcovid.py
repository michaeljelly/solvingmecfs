"""Long COVID Olink datasets as an outside comparison: does the ME/CFS protein picture appear in long COVID?"""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
Z=W+'raw/zen/'
def save(ds,X,meta,feat):
    X=X.astype(float); X.columns=[str(c) for c in X.columns]; X.to_parquet(D+f'{ds}_X.parquet'); meta.loc[X.index].to_csv(D+f'{ds}_meta.csv'); feat.loc[X.columns].to_csv(D+f'{ds}_feat.csv')
    print(ds,X.shape,meta.loc[X.index].case.value_counts().to_dict())
# 1. UK Olink Explore 3072 (Zenodo 14772494): CA long COVID vs CO controls
parts=[]
for f in ['plate1_olink.csv','plate2_olink.csv','plate3_olink.csv']:
    p=pd.read_csv(Z+f,sep=None,engine='python',usecols=['SampleID','Sample_Type','OlinkID','Assay','NPX','QC_Warning','UniProt'])
    p=p[(p.Sample_Type=='SAMPLE')&(p.SampleID.str[:2].isin(['CA','CO']))&(p.QC_Warning.str.upper()=='PASS')]; parts.append(p)
p=pd.concat(parts); X=p.pivot_table(index='SampleID',columns='OlinkID',values='NPX',aggfunc='mean')
ann=p.drop_duplicates('OlinkID').set_index('OlinkID')
meta=pd.read_csv(Z+'meta_UKsamples.csv'); meta['sid']=meta.SampleID.str[:2]+meta.SampleID.str[2:].astype(int).astype(str); meta=meta.set_index('sid').reindex(X.index)
meta['case']=(X.index.str[:2]=='CA').astype(int); meta['male']=(meta.Sex=='m').astype(int); meta['age']=meta.Age; meta['bmi']=meta.BMI
save('lc_uk_explore',X,meta,pd.DataFrame({'gene':ann.Assay,'uniprot':ann.UniProt},index=ann.index))
# 2. IRIS (Zenodo 13363117): long COVID vs recovered, 3 months
d=pd.read_csv(Z+'IRIS_LongCOVID_olink_data.csv'); d=d[d.timing_month=='3 months']
X=d.pivot_table(index='IRIS_number',columns='Assay',values='NPX',aggfunc='mean'); X.index=X.index.astype(str)
m=d.drop_duplicates('IRIS_number').set_index('IRIS_number'); m.index=m.index.astype(str)
m['case']=m.long_covid.astype(int); m['male']=(m.Sex=='Male').astype(int)
save('lc_iris_3mo',X,m[['case','male']],pd.DataFrame({'gene':X.columns,'uniprot':''},index=X.columns))
# 3. Bansal/Olechnowicz (Zenodo 15237088): baseline, long COVID vs healthy
x=pd.read_excel(Z+'baseline57.xlsx',header=0)
x=x[x['QC_Warning'].astype(str).str.lower()=='pass'] if 'QC_Warning' in x else x
X=x.pivot_table(index='Study Participant',columns='Assay',values='NPX',aggfunc='mean'); X.index=X.index.astype(str)
m=x.drop_duplicates('Study Participant').set_index('Study Participant'); m.index=m.index.astype(str)
print(m['Long_COVID status'].value_counts().to_dict())
m=m[m['Long_COVID status'].isin(['Long COVID','Healthy','Recovered'])]
m['case']=(m['Long_COVID status']=='Long COVID').astype(int); m['male']=(m.gender_at_birth=='Male').astype(int); m['age']=m.age
save('lc_bansal',X.loc[m.index],m[['case','male','age','Long_COVID status']],pd.DataFrame({'gene':X.columns,'uniprot':''},index=X.columns))
# case-control
for ds,cov in [('lc_uk_explore',('male','age','bmi')),('lc_iris_3mo',('male',)),('lc_bansal',('male','age'))]:
    X,m,f=load(ds); G=gene_collapse(X,f); des=m[['case']+list(cov)]; ok=des.notna().all(1).values
    r=ols_many(G[ok],des[ok],'case'); r['n_case']=int(des[ok].case.sum()); r['n_ctrl']=int((des[ok].case==0).sum()); r.index.name='gene'; r.to_csv(R+f'cc_{ds}.csv')
    print(ds,'q<0.05',(r.q<0.05).sum(),'p<0.05',(r.p<0.05).sum(),'of',len(r),'| top:',', '.join(f"{g}({b:+.2f},p={p:.1g})" for g,b,p in r.sort_values('p').head(10)[['beta','p']].itertuples()))
