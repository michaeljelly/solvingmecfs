"""SomaScan 1.3k proteomics from GEO (GSE251790 CSF, GSE254030 serum): group tests and project-lead proteins (Th17 axis, kynurenine, catecholamine enzymes)."""
import pandas as pd, numpy as np, gzip, re, warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
B = '/mnt/project-files/mecfs-research/downloads/nih_walitt'

def read_adat(path):
    lines = gzip.open(path, 'rt').read().split('\n')
    i0 = next(i for i, l in enumerate(lines) if l.startswith('^TABLE_BEGIN'))
    ann_lines = []
    j = i0 + 1
    while lines[j].startswith('\t'):
        parts = lines[j].split('\t'); ann_lines.append(parts); j += 1
    k = next(i for i, p in enumerate(ann_lines[0]) if p == 'SeqId')
    ann = pd.DataFrame({p[k]: p[k + 1:] for p in ann_lines})
    header = lines[j].split('\t')
    rows = [l.split('\t') for l in lines[j + 1:] if l.strip()]
    meta = pd.DataFrame([r[:k] for r in rows], columns=header[:k])
    vals = pd.DataFrame([r[k + 1:] for r in rows]).apply(pd.to_numeric, errors='coerce')
    vals.columns = ann.SeqId.values
    return meta, vals, ann

out = []
for tissue, f in [('CSF', 'GSE251790_adat.txt.gz'), ('serum', 'GSE254030_adat.txt.gz')]:
    meta, vals, ann = read_adat(f'{B}/raw/geo/{f}')
    keep = (meta.SampleType == 'Sample').values
    meta, vals = meta[keep].reset_index(drop=True), vals[keep].reset_index(drop=True)
    grpcol = [c for c in meta.columns if meta[c].isin(['Patient', 'Control']).any()][0]
    g = meta[grpcol]
    prot = ann.Type.eq('Protein').values
    V = np.log2(vals.loc[:, prot].clip(lower=1)); A = ann[prot].reset_index(drop=True)
    print(tissue, 'samples', dict(g.value_counts()), 'proteins', V.shape[1])
    for i, c in enumerate(V.columns):
        a = V.loc[g == 'Control', c]; b = V.loc[g == 'Patient', c]
        out.append(dict(tissue=tissue, SeqId=c, Target=A.Target[i], gene=A.EntrezGeneSymbol[i], log2FC=b.mean() - a.mean(),
                        p=stats.mannwhitneyu(a, b).pvalue))
res = pd.DataFrame(out)
res['q'] = res.groupby('tissue').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
res.to_csv(f'{B}/results/I_somalogic_all_tests.csv', index=False)
for t, d in res.groupby('tissue'):
    print(t, 'p<0.05:', (d.p < 0.05).sum(), 'of', len(d), '(expected by chance ~%d); q<0.10:' % (0.05 * len(d)), (d.q < 0.10).sum())
    print(d.sort_values('p').head(8)[['Target', 'gene', 'log2FC', 'p', 'q']].round(4).to_string())
leads = r'^(IL17A|IL17F|IL17RA|IL17RC|IL17B|IL17D|IL22|IL23A|IL12B|CCL20|IL6|IL1B|TNF|IFNG|IDO1|KYNU|KMO|TDO2|DDC|DBH|TH|COMT|MAOA|MAOB|PNMT|FGF21|GDF15|CRP|RORC|IL21|IL10)\b'
lead = res[res.gene.fillna('').str.contains(leads)]
lead.to_csv(f'{B}/results/I_somalogic_lead_proteins.csv', index=False)
print('\nlead proteins:'); print(lead.sort_values(['tissue', 'p'])[['tissue', 'Target', 'gene', 'log2FC', 'p', 'q']].round(4).to_string())
