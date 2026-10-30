"""Load every per-sample ME/CFS protein dataset into a common format.
Output per dataset in data/: <ds>_X.parquet (samples x features, log2 scale), <ds>_meta.csv (index = sample),
<ds>_feat.csv (feature -> gene symbol, uniprot).
"""
import pandas as pd, numpy as np, openpyxl, gzip, re, warnings
warnings.filterwarnings('ignore')
W = '/tmp/work_proteome/'; D = W + 'data/'
LOC = '/mnt/project-files/mecfs-research/downloads/'

def rows(f, s=None):
    wb = openpyxl.load_workbook(f, read_only=True)
    ws = wb[s] if s else wb.worksheets[0]
    return [r for r in ws.iter_rows(values_only=True)]

def save(ds, X, meta, feat):
    X = X.astype(float); X.columns = [str(c) for c in X.columns]
    meta = meta.loc[X.index]; feat = feat.loc[X.columns]
    X.to_parquet(D + f'{ds}_X.parquet'); meta.to_csv(D + f'{ds}_meta.csv'); feat.to_csv(D + f'{ds}_feat.csv')
    print(ds, X.shape, meta.get('case', pd.Series(dtype=float)).value_counts().to_dict())

import os
RUN_ALL = os.environ.get('RUN_ALL')=='1'
# ---------- 1. Hoel/Fluge SomaScan 7k serum (PAD000026; supplement SuppData1 has clinical columns)
r = rows(W + 'raw/supp/PMC13006441/mmc2.xlsx', '2_log2Matrix')
df = pd.DataFrame(r[1:], columns=[str(c) for c in r[0]])
apt = [c for c in df.columns if c.startswith('seq.')]
meta = df[[c for c in df.columns if not c.startswith('seq.')]].copy()
meta.index = [f'hoel{i}' for i in range(len(df))]
meta['case'] = (meta.Group == 'ME').astype(int)
meta['male'] = (meta.Sex.astype(str) == '0').astype(int)  # Sex=1 is the majority (female) code
agemap = {a: i for i, a in enumerate(sorted(meta.Age.astype(str).unique()))}
meta['age'] = meta.Age.astype(str).map(agemap)
meta['bmi'] = pd.to_numeric(meta.BMI, errors='coerce')
X = df[apt].apply(pd.to_numeric, errors='coerce'); X.index = meta.index
r = rows(W + 'raw/supp/PMC13006441/mmc2.xlsx', '3_Aptamer_annotations')
ann = pd.DataFrame(r[1:], columns=r[0]).astype(str).set_index('AptName')
feat = pd.DataFrame({'gene': ann.EntrezGeneSymbol, 'uniprot': ann.UniProt, 'target': ann.Target}, index=ann.index)
save('hoel_serum', X, meta, feat)

# ---------- 2. Germain 2025 SomaScan 7k plasma, 2-day CPET, 4 timepoints
S = pd.read_excel(W + 'raw/supp/PMC12757482/mmc1.xlsx', sheet_name='Somalogic_dataset', engine='openpyxl') if False else None
r = rows(W + 'raw/supp/PMC12757482/mmc1.xlsx', 'Somalogic_dataset')
df = pd.DataFrame(r[1:], columns=[str(c) for c in r[0]])
demo = pd.read_csv(D + 'germain_Demographics.csv').set_index('cor_id')
cb = pd.read_csv(D + 'germain_Somalogic_codebook.csv').set_index('SeqId')
df.index = df.cor_id + '_' + df.timepoint
meta = df[['cor_id', 'timepoint']].join(demo, on='cor_id')
meta['case'] = (meta.phenotype == 'ME/CFS').astype(int); meta['male'] = (meta.sex == 'M').astype(int)
meta['bmi_bin'] = meta.bmi
X = df.drop(columns=['cor_id', 'timepoint']).apply(pd.to_numeric, errors='coerce')
print('germain value range', np.nanpercentile(X.values, [1, 50, 99]))
if np.nanmedian(X.values) > 100: X = np.log2(X)
feat = pd.DataFrame({'gene': cb.EntrezGeneSymbol.reindex(X.columns).values, 'uniprot': cb.UniProt.reindex(X.columns).values,
                     'target': cb.Target.reindex(X.columns).values}, index=X.columns)
save('germain_plasma', X, meta, feat)

# ---------- 3/4. Walitt NIH SomaScan 1.3k serum and CSF
def read_adat(path):
    lines = gzip.open(path, 'rt').read().split('\n')
    i0 = next(i for i, l in enumerate(lines) if l.startswith('^TABLE_BEGIN'))
    ann_lines = []; j = i0 + 1
    while lines[j].startswith('\t'):
        ann_lines.append(lines[j].split('\t')); j += 1
    k = next(i for i, p in enumerate(ann_lines[0]) if p == 'SeqId')
    ann = pd.DataFrame({p[k]: p[k + 1:] for p in ann_lines})
    header = lines[j].split('\t'); rr = [l.split('\t') for l in lines[j + 1:] if l.strip()]
    meta = pd.DataFrame([x[:k] for x in rr], columns=header[:k])
    vals = pd.DataFrame([x[k + 1:] for x in rr]).apply(pd.to_numeric, errors='coerce'); vals.columns = ann.SeqId.values
    return meta, vals, ann
for tissue, f in [('csf', 'GSE251790_adat.txt.gz'), ('serum', 'GSE254030_adat.txt.gz')]:
    m, v, a = read_adat(LOC + 'nih_walitt/raw/geo/' + f)
    keep = (m.SampleType == 'Sample').values; m = m[keep].reset_index(drop=True); v = v[keep].reset_index(drop=True)
    prot = a.Type.eq('Protein').values; v = np.log2(v.loc[:, prot].clip(lower=1)); a = a[prot].reset_index(drop=True)
    m.index = m.SampleId.values; v.index = m.index
    gc = [c for c in m.columns if m[c].isin(['Patient', 'Control']).any()][0]; m['Group'] = m[gc]
    sc = [c for c in m.columns if m[c].isin(['M', 'F', 'Male', 'Female']).any()][0]; m['Gender'] = m[sc].str[0]
    ac = [c for c in m.columns if c.lower() == 'age']; m['Age'] = m[ac[0]] if ac else np.nan
    m['case'] = (m.Group == 'Patient').astype(int); m['male'] = (m.Gender == 'M').astype(int)
    m['age'] = pd.to_numeric(m.Age, errors='coerce')
    feat = pd.DataFrame({'gene': a.EntrezGeneSymbol.values, 'uniprot': a.UniProt.values, 'target': a.Target.values}, index=a.SeqId.values)
    save(f'walitt_{tissue}', v, m[['SampleId', 'Group', 'Gender', 'case', 'male', 'age', 'PlateId', 'PlatePosition']], feat)

# ---------- 5. Rydland 2026 plasma EV (SEC + timsTOF, MaxQuant LFQ, 178 proteins; PXD072203)
r = rows(W + 'raw/supp/PMC13312568/mmc1.xlsx')
hdr = list(r[1]); body = pd.DataFrame(r[2:], columns=[f'c{i}' for i in range(len(hdr))])
grp = hdr[:98]
Xt = body.iloc[:, :98].apply(pd.to_numeric, errors='coerce')
genes = body.iloc[:, hdr.index('Gene names')].astype(str).str.split(';').str[0]
uni = body.iloc[:, hdr.index('Majority protein IDs')].astype(str).str.split(';').str[0]
X = Xt.T; X.index = [f'ev{i}' for i in range(98)]; X.columns = [f'p{i}' for i in range(len(body))]
meta = pd.DataFrame({'group': grp, 'case': [1 if g == 'ME' else 0 for g in grp]}, index=X.index)
feat = pd.DataFrame({'gene': genes.values, 'uniprot': uni.values}, index=X.columns)
keep = feat.gene.ne('None') & feat.gene.ne('nan'); X = X.loc[:, keep]
save('rydland_ev', X, meta, feat)

# ---------- 6. Giloteaux 2024 EV (TMT, 0h/15min/24h after CPET)
r = rows(W + 'raw/supp/PMC10764978/JEV2-13-12403-s002.xlsx')
uni = r[3][4:]; gsym = r[5][4:]; hdr = r[6]
body = pd.DataFrame(r[7:], columns=[str(h) for h in hdr]).dropna(subset=['ID'])
X = body.iloc[:, 4:].apply(pd.to_numeric, errors='coerce')
X.columns = [f'{g}|{i}' for i, g in enumerate(gsym)]
print('giloteaux value range', np.nanpercentile(X.values, [1, 50, 99]))
X = np.log2(X.clip(lower=1e-3))
X.index = body.ID.astype(int).astype(str) + '_' + body.time.astype(str)
meta = pd.DataFrame({'pid': body.ID.astype(int).astype(str).values, 'time': body.time.values, 'phenotype': body.phenotype.values,
                     'tmt': body['TMT _set'].values}, index=X.index)
meta['case'] = (meta.phenotype == 'ME/CFS').astype(int)
feat = pd.DataFrame({'gene': list(gsym), 'uniprot': list(uni)}, index=X.columns)
save('giloteaux_ev', X, meta, feat)

# ---------- 7. Bragee 2026 CSF (patients only; POTS labels)
r = rows(W + 'raw/supp/PMC13194927/41598_2026_46965_MOESM1_ESM.xlsx', 'VSN_ProteinQuant')
df = pd.DataFrame(r[1:], columns=[str(c) for c in r[0]]).set_index('protein')
X = df.apply(pd.to_numeric, errors='coerce').T
meta = pd.DataFrame(index=X.index)
meta['severity'] = [c.split('_')[1] for c in X.index]; meta['pots'] = [0 if c.endswith('NoPOTS') else 1 for c in X.index]
meta['case'] = 1
feat = pd.DataFrame({'gene': X.columns, 'uniprot': ''}, index=X.columns)
save('bragee_csf', X, meta, feat)

# ---------- 8. PXD065439 IgG-bound serum proteins (MaxQuant LFQ)
pg = pd.read_csv(W + 'raw/pride/pxd065439/txt/proteinGroups.txt', sep='\t', low_memory=False)
pg = pg[(pg['Reverse'] != '+') & (pg['Potential contaminant'] != '+') & (pg['Only identified by site'] != '+')]
lfq = [c for c in pg.columns if c.startswith('LFQ intensity ')]
if not lfq: lfq = [c for c in pg.columns if c.startswith('Intensity ') and c != 'Intensity']
X = np.log2(pg[lfq].replace(0, np.nan)).T
X.columns = [f'pg{i}' for i in range(len(pg))]
X.index = [c.replace('LFQ intensity ', '').replace('Intensity ', '') for c in lfq]
meta = pd.DataFrame(index=X.index)
meta['group'] = ['HC' if 'Healthy' in s else ('PCS' if 'PCSCFS' in s else 'ME') for s in X.index]
meta['batch'] = [s.split('_')[0] for s in X.index]
meta['case'] = (meta.group != 'HC').astype(int)
feat = pd.DataFrame({'gene': pg['Gene names'].astype(str).str.split(';').str[0].values,
                     'uniprot': pg['Majority protein IDs'].astype(str).str.split(';').str[0].values}, index=X.columns)
save('igg_serum', X, meta, feat)
