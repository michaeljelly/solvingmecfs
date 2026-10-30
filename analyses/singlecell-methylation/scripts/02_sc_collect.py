# Collect per-cell tables -> composition counts, per-sample x celltype mean scores; merge metadata.
import pandas as pd, numpy as np, glob, os
O='/tmp/work_scm/sc_out'; R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
meta=pd.read_csv('/mnt/project-files/mecfs-research/downloads/geo/data/processed/GSE214283_meta.csv')
meta['person']=meta.title.str.replace(r'-D[12]$','',regex=True); meta['day']=meta['collection day']
cells=[]
for f in sorted(glob.glob(f'{O}/*_cells.parquet')):
    c=pd.read_parquet(f); cells.append(c)
C=pd.concat(cells); C=C.reset_index().rename(columns={'index':'barcode'})
MAP={'Tcm/Naive helper T cells':'CD4_naive_CM','Tem/Effector helper T cells':'CD4_Tem','Tem/Effector helper T cells PD1+':'CD4_Tem_PD1',
 'Regulatory T cells':'Treg','Treg(diff)':'Treg','Type 1 helper T cells':'Th1','Type 17 helper T cells':'Th17','Follicular helper T cells':'Tfh',
 'Memory CD4+ cytotoxic T cells':'CD4_CTL','Tcm/Naive cytotoxic T cells':'CD8_naive_CM','Tem/Trm cytotoxic T cells':'CD8_Tem','Tem/Temra cytotoxic T cells':'CD8_Temra',
 'Trm cytotoxic T cells':'CD8_Tem','MAIT cells':'MAIT','gamma-delta T cells':'gdT','CRTAM+ gamma-delta T cells':'gdT','NKT cells':'NKT','Cycling T cells':'Cycling_T','Cycling NK cells':'Cycling_T',
 'CD8a/a':'CD8_other','CD8a/b(entry)':'CD8_other','Double-positive thymocytes':'T_other','Double-negative thymocytes':'T_other','T(agonist)':'T_other','Early lymphoid/T lymphoid':'T_other',
 'Naive B cells':'B_naive','Memory B cells':'B_memory','Age-associated B cells':'B_ABC','Transitional B cells':'B_naive','Plasmablasts':'Plasmablast','Plasma cells':'Plasmablast','Germinal center B cells':'B_memory','Cycling B cells':'Plasmablast','Follicular B cells':'B_naive','Proliferative germinal center B cells':'Plasmablast',
 'pDC':'pDC','pDC precursor':'pDC','DC2':'cDC2','DC1':'cDC1','DC':'cDC2','DC3':'cDC2','Transitional DC':'cDC2','Migratory DCs':'cDC2','Megakaryocytes/platelets':'Platelet','Early MK':'Platelet','Megakaryocyte precursor':'Platelet',
 'HSC/MPP':'HSPC','CMP':'HSPC','GMP':'HSPC','MEMP':'HSPC','Megakaryocyte-erythroid-mast cell progenitor':'HSPC','ILC':'ILC','ILC3':'ILC','ILC1':'ILC','ILC2':'ILC','ILC precursor':'ILC','Mast cells':'HSPC',
 'Erythrocytes':'Erythroid','Late erythroid':'Erythroid','Mid erythroid':'Erythroid','Early erythroid':'Erythroid','Neutrophils':'Granulocyte','Myelocytes':'Granulocyte','Promyelocytes':'Granulocyte','Granulocytes':'Granulocyte'}
f=C['fine'].astype(str)
C['type']=f.map(lambda x: x if x.startswith(('Mono_','NK_')) or x=='Doublet' else MAP.get(x,'Other'))
print(C.type.value_counts().to_string())
print('unmapped ->Other:',C.loc[C.type=='Other','fine'].value_counts().head(15).to_dict())
C.to_parquet('/tmp/work_scm/cells_all.parquet')
cnt=pd.crosstab(C.gsm,C.type); cnt.to_csv(f'{R}/sc_celltype_counts_per_sample.csv')
meta.to_csv(f'{R}/sc_sample_meta.csv',index=False)
