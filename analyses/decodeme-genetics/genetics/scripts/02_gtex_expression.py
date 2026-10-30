"""GTEx v8 median TPM for DecodeME candidate genes (all genes within 500 kb + Tier-1 eQTL genes).
Source file: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz (github mirror laurenfries/Adventitial_scRNAseq).
Reports top tissues and a tau tissue-specificity index."""
import pandas as pd, numpy as np
g = pd.read_csv('/mnt/project-files/mecfs-research/genetics/data/g.gz', sep='\t', skiprows=2)
g = g.drop(columns=['Name']).groupby('Description').max()
genes = ['RABGAP1L','DARS2','CENPL','ZBTB37','SERPINC1','KLHL20','RC3H1','ANKRD45','TNFSF4','PRDX6',
 'BTN2A2','BTN2A1','BTN3A1','BTN3A2','BTN1A1','TRIM38','ABT1','ZNF322',
 'FBXL4','POU3F2','MMS22L',
 'SUDS3','TAOK3','PEBP1','KSR2',
 'OLFM4','PCDH8',
 'CCPG1','UNC13C','RAB27A','PIGB','RSL24D1',
 'CA10',
 'ARFGEF2','CSE1L','PREX1','STAU1','ZNFX1','DDX27','KCNB1']
sel = g.reindex([x for x in genes if x in g.index])
missing = [x for x in genes if x not in g.index]
def tau(x):
    x = np.log2(x+1); m = x.max()
    return float(((1-x/m).sum()/(len(x)-1))) if m>0 else np.nan
brain = [c for c in g.columns if c.startswith('Brain')]
out=[]
for gene,row in sel.iterrows():
    top = row.sort_values(ascending=False).head(4)
    out.append(dict(gene=gene, max_tpm=round(top.iloc[0],1), tau=round(tau(row),2),
        brain_mean=round(row[brain].mean(),1), whole_blood=round(row['Whole Blood'],1),
        EBV_lymphocytes=round(row['Cells - EBV-transformed lymphocytes'],1),
        muscle=round(row['Muscle - Skeletal'],1), tibial_nerve=round(row['Nerve - Tibial'],1),
        colon_transverse=round(row['Colon - Transverse'],1), small_intestine=round(row['Small Intestine - Terminal Ileum'],1),
        top_tissues='; '.join(f"{k.replace('Brain - ','Br:')} {v:.0f}" for k,v in top.items())))
df = pd.DataFrame(out)
df.to_csv('/mnt/project-files/mecfs-research/genetics/candidate_gene_gtex_expression.csv', index=False)
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',90)
print(df.to_string(index=False)); print('missing:',missing)
