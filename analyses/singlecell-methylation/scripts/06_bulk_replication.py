# Replicate scRNA leads in independent bulk RNA cohorts: mitochondrial-transcript fraction, CD8 effector (KLRG1/CX3CR1), adaptive NK (KLRC2), ABC (FCRL5/ITGAX/TBX21), platelet (PPBP/PF4).
import pandas as pd, numpy as np, gzip, statsmodels.formula.api as smf
RAW='/tmp/work_geo/raw'; P='/mnt/project-files/mecfs-research/downloads/geo/data/processed'; R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
MT_SYM=['ND1','ND2','ND3','ND4','ND4L','ND5','ND6','COX1','COX2','COX3','ATP6','ATP8','CYTB']
MT_ENS=['ENSG00000198888','ENSG00000198763','ENSG00000198840','ENSG00000198886','ENSG00000212907','ENSG00000198786','ENSG00000198695','ENSG00000198804','ENSG00000198712','ENSG00000198938','ENSG00000198899','ENSG00000228253','ENSG00000198727']
MT_ID=[4535,4536,4537,4538,4539,4540,4541,4512,4513,4514,4508,4509,4519]
MARK={'KLRG1':('ENSG00000139187',10219),'CX3CR1':('ENSG00000168329',1524),'KLRC2':('ENSG00000205809',3822),'FCRL5':('ENSG00000143297',83416),'ITGAX':('ENSG00000140678',3687),
      'TBX21':('ENSG00000073861',30009),'PPBP':('ENSG00000163736',5473),'PF4':('ENSG00000163737',5196),'ZNF683':('ENSG00000176083',257101),'FCER1G':('ENSG00000158869',2207),'GZMB':('ENSG00000100453',3002)}
def smeta(g):
    meta={}
    for l in gzip.open(f'{RAW}/{g}_series_matrix.txt.gz','rt'):
        if l.startswith('!Sample_'):
            p=l.rstrip('\n').split('\t'); k=p[0][8:]
            while k in meta: k+='_'
            meta[k]=[x.strip('"') for x in p[1:]]
        if l.startswith('"ID_REF"') or l.startswith('!series_matrix_table_begin'): break
    return pd.DataFrame(meta)
def feats(cnt,mt,markers):
    lib=cnt.sum(0); f=pd.DataFrame({'mt_frac':cnt.loc[cnt.index.intersection(mt)].sum(0)/lib})
    cpm=np.log2(cnt/lib*1e6+1)
    for k,v in markers.items():
        if v in cpm.index: f[k]=cpm.loc[v]
    return f
rows=[]
def test(ds,F,covar=''):
    for c in F.columns:
        if c in('case','sex','subject','age','celltype'): continue
        d=F.dropna(subset=[c]).copy(); d['z']=(d[c]-d[c].mean())/d[c].std()
        if d.z.std()==0: continue
        m=smf.ols('z~case'+covar,d).fit()
        rows.append(dict(dataset=ds,feature=c,n_case=int(d.case.sum()),n_ctrl=int((d.case==0).sum()),mean_case=d.loc[d.case==1,c].mean(),mean_ctrl=d.loc[d.case==0,c].mean(),d=m.params['case'],p=m.pvalues['case']))
# GSE227375 Miami PBMC, baseline T0 (one per person)
cnt=pd.read_parquet(f'{P}/GSE227375_counts_symbol.parquet'); m=pd.read_csv(f'{P}/GSE227375_meta.csv',index_col=0)
m=m[m.time==0]; F=feats(cnt[m.index],MT_SYM,{k:k for k in MARK}); F['case']=m.case; F['sex']=m.sex
test('GSE227375_PBMC_T0',F,'+C(sex)')
# all timepoints person mean
F=feats(cnt,MT_SYM,{k:k for k in MARK}); mm=pd.read_csv(f'{P}/GSE227375_meta.csv',index_col=0); F=F.join(mm[['case','sex','subject']]); Fp=F.groupby('subject').mean(numeric_only=True); Fp['sex']=F.groupby('subject').sex.first()
test('GSE227375_PBMC_personmean',Fp,'+C(sex)')
# GSE268212 sorted CD8
cnt=pd.read_csv(f'{RAW}/GSE268212_mecfs_cd8_rna.rawCounts.csv.gz',index_col=0); cnt.index=cnt.index.str.split('.').str[0]; cnt=cnt.groupby(level=0).sum()
M=smeta('GSE268212'); M['col']=M['description']; M=M.set_index('col')
ct=M['characteristics_ch1'].str.split(': ').str[1]; ph=M['characteristics_ch1_'].str.split(': ').str[1]
for cellt in ['NVCD8','EMCD8']:
    cols=[c for c in ct.index[ct==cellt] if c in cnt.columns]
    F=feats(cnt[cols],MT_ENS,{k:v[0] for k,v in MARK.items()}); F['case']=(ph[cols]=='ME').astype(int).values
    test(f'GSE268212_{cellt}',F)
# GSE98139 adolescent whole blood (normalized counts)
cnt=pd.read_csv(f'{RAW}/GSE98139_All_normalized_counts.csv.gz',index_col=0,sep=None,engine='python'); cnt.index=cnt.index.str.split('.').str[0]
M=smeta('GSE98139'); M['id']=M.title.str.extract(r'(\d+)')[0]; M=M.set_index('id')
cond=M['characteristics_ch1'].str.split(': ').str[1]
cols=[c for c in cnt.columns if c in cond.index]
F=feats(cnt[cols],MT_ENS,{k:v[0] for k,v in MARK.items()}); F['case']=(cond[cols].str.contains('CFS')).astype(int).values
test('GSE98139_WB_adolescent',F)
# GSE128078 whole blood day1
cnt=pd.read_csv(f'{RAW}/GSE128078_ncbi_counts.tsv.gz',sep='\t',index_col=0)
M=smeta('GSE128078').set_index('geo_accession'); ds=M['characteristics_ch1_'].str.split(': ').str[1]; day=M['characteristics_ch1___'].str.split(': ').str[1]
cols=[c for c in cnt.columns if c in M.index and day[c]=='1']
F=feats(cnt[cols],MT_ID,{k:v[1] for k,v in MARK.items()}); F['case']=(ds[cols].str.contains('ME')).astype(int).values
test('GSE128078_WB_day1',F)
# GSE270045 long COVID whole blood
cnt=pd.read_csv(f'{RAW}/GSE270045_LC_counts.tsv.gz',sep='\t',index_col=0)
M=smeta('GSE270045'); sid=M.filter(like='characteristics').apply(lambda r:[x for x in r if str(x).startswith('sample id')][0].split(': ')[1],axis=1)
M['sid']=sid; M=M.set_index('sid'); cols=[c for c in cnt.columns if c in M.index]
F=feats(cnt[cols],['MT-'+x for x in ['ND1','ND2','ND3','ND4','ND4L','ND5','ND6','CO1','CO2','CO3','ATP6','ATP8','CYB']],{k:k for k in MARK})
F['case']=(~M.loc[cols,'title'].str.contains('Healthy')).astype(int).values
test('GSE270045_WB_longCOVID',F)
O=pd.DataFrame(rows); O.to_csv(f'{R}/bulk_replication_scRNA_leads.csv',index=False)
pd.set_option('display.width',200); print(O[O.feature.isin(['mt_frac','KLRG1','CX3CR1','KLRC2','FCRL5','ZNF683','PPBP'])].round(4).to_string())
