"""Exploratory locus-level tissue-specificity test (a crude stand-in for MAGMA gene-property analysis,
because DecodeME summary statistics could not be downloaded here).
For each of the 8 DecodeME lead variants take protein-coding genes within +/-250 kb (fallback: nearest gene).
Per gene: tissue specificity = log2(TPM+1) z-scored across 54 GTEx tissues. Per locus: max z per tissue.
Statistic: mean over loci. Null: 20,000 draws of 8 random genome positions (same rule, positions spaced >1 Mb).
Caveat: 8 loci, no LD modelling, genes chosen by position only. Hypothesis-generating."""
import pandas as pd, numpy as np
rng = np.random.default_rng(1)
loc = pd.read_csv('/mnt/project-files/mecfs-research/genetics/data/NCBI38.gene.loc', sep='\t', header=None,
                  names=['entrez','chr','start','end','strand','gene'])
loc = loc[~loc.gene.str.startswith(('LOC','LINC','MIR','SNOR'))]
loc = loc[loc.chr.astype(str).isin([str(i) for i in range(1,23)])]
loc['chr']=loc.chr.astype(int)
g = pd.read_csv('/mnt/project-files/mecfs-research/genetics/data/g.gz', sep='\t', skiprows=2).drop(columns=['Name']).groupby('Description').max()
g = g[g.max(axis=1)>=1]
L = np.log2(g+1); Z = L.sub(L.mean(1),axis=0).div(L.std(1)+1e-9,axis=0)
loc = loc[loc.gene.isin(Z.index)].reset_index(drop=True)
Zl = Z.loc[loc.gene].values
bychr = {c: d for c,d in loc.groupby('chr')}
def locus_vec(c,pos,w=250_000):
    d = bychr[c]; idx = d.index[(d.end>pos-w)&(d.start<pos+w)]
    if len(idx)==0:
        dist = np.minimum(abs(d.start-pos),abs(d.end-pos)); idx=[dist.idxmin()]
    return Zl[idx].max(0)
leads=[(1,173846152),(6,26239176),(6,97984426),(12,118202773),(13,53194927),(15,54866724),(17,52183006),(20,48914387)]
obs = np.mean([locus_vec(c,p) for c,p in leads],0)
# chromosome sampling weighted by gene span; exclude MHC for null? keep simple: sample uniformly within gene-bearing extent
ext = loc.groupby('chr').agg(lo=('start','min'),hi=('end','max')); w = (ext.hi-ext.lo)/ (ext.hi-ext.lo).sum()
N=20000; null=np.zeros((N,Z.shape[1]))
for i in range(N):
    cs = rng.choice(ext.index, size=8, p=w.values)
    vs=[locus_vec(c, int(rng.integers(ext.lo[c], ext.hi[c]))) for c in cs]
    null[i]=np.mean(vs,0)
p = ((null>=obs).sum(0)+1)/(N+1)
from statsmodels.stats.multitest import multipletests
q = multipletests(p,method='fdr_bh')[1]
res = pd.DataFrame(dict(tissue=Z.columns, obs=obs.round(2), null_mean=null.mean(0).round(2), p=p, q=q)).sort_values('p')
res.to_csv('/mnt/project-files/mecfs-research/genetics/locus_tissue_permutation.csv', index=False)
print(res.head(15).to_string(index=False))
