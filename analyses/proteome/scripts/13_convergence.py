"""Cross-omics pathway convergence: proteins (blood meta), blood RNA (4-cohort meta), DecodeME MAGMA genes, metabolites (metabolome-wide meta).
Per layer and pathway: competitive enrichment (unsigned: mean z^2 for protein/RNA/metabolite; mean MAGMA Z for genetics) vs random same-size sets.
Convergence: (1) pathway-score correlation between layers, null = gene-label permutation of one layer (keeps pathway overlap structure);
(2) number of pathways in the top 10% of >=2 / >=3 layers vs the same null; (3) Fisher-combined pathway ranking."""
import sys, re; sys.path.insert(0, '/tmp/work_proteome/scripts')
from common import *
rng = np.random.default_rng(2024)
NPERM = 300
# ---------------- layers
pm = pd.read_csv(R + 'meta_blood_genes.csv', index_col=0)
prot = (pm.re / pm.re_se).dropna(); prot = prot[~prot.index.duplicated()]
rn = pd.read_csv(R + 'rna_meta_genomewide.csv', index_col=0); rna = rn.z[rn.k >= 2].dropna()
mg = pd.read_csv('/mnt/project-files/mecfs-research/downloads/decodeme/results/magma_gwas_1_genes.csv')
magma = mg.drop_duplicates('SYMBOL').set_index('SYMBOL').ZSTAT.dropna()
met = pd.read_csv('/mnt/project-files/mecfs-research/deep2/metabolome_wide/meta_all.csv')
met = met[met.all_k >= 2].copy(); met['z'] = met.all_est / met.all_se
# ---------------- metabolite -> KEGG compound by name
K = W + 'raw/kegg/'
norm = lambda s: re.sub(r'[^a-z0-9]', '', str(s).lower())
cn = {}
for l in open(K + 'list_compound.txt'):
    cid, names = l.rstrip('\n').split('\t'); cid = cid.replace('cpd:', '')
    for n in names.split(';'):
        cn.setdefault(norm(n), cid)
def match(lbl):
    for v in [lbl, re.sub(r'\*$', '', lbl), re.sub(r'^(l|d|dl)-', '', lbl, flags=re.I), lbl.replace('ate', 'ic acid')]:
        if norm(v) in cn: return cn[norm(v)]
    return None
met['kegg'] = met.label.astype(str).map(match)
print('metabolites', len(met), 'mapped to KEGG', met.kegg.notna().sum())
metz = met.dropna(subset=['kegg']).groupby('kegg').z.mean()
# ---------------- pathways (KEGG, human)
pname = dict(l.rstrip('\n').split('\t') for l in open(K + 'list_pathway_hsa.txt'))
pname = {k.replace('path:', ''): v.replace(' - Homo sapiens (human)', '') for k, v in pname.items()}
e2s = {}
for l in open(K + 'list_hsa.txt'):
    p = l.rstrip('\n').split('\t')
    if len(p) >= 4: e2s[p[0]] = p[3].split(',')[0].split(';')[0].strip()
kg = {}
for l in open(K + 'link_pathway_hsa.txt'):
    g, p = l.rstrip('\n').split('\t'); p = p.replace('path:', '')
    if g in e2s: kg.setdefault(p, set()).add(e2s[g])
kc = {}
for l in open(K + 'link_pathway_compound.txt'):
    c, p = l.rstrip('\n').split('\t'); p = p.replace('path:map', 'hsa')
    if p in kg: kc.setdefault(p, set()).add(c.replace('cpd:', ''))
# ---------------- enrichment
def enrich(score, sets, stat, nperm=2000, min_n=5):
    v = score.values; idx = {g: i for i, g in enumerate(score.index)}; out = {}
    tr = v ** 2 if stat == 'sq' else v
    for s, genes in sets.items():
        ii = [idx[g] for g in genes if g in idx]
        if len(ii) < min_n: continue
        obs = tr[ii].mean()
        null = np.array([tr[rng.choice(len(tr), len(ii), replace=False)].mean() for _ in range(nperm)])
        out[s] = dict(n=len(ii), obs=obs, z=(obs - null.mean()) / null.std(), p=(np.sum(null >= obs) + 1) / (nperm + 1))
    return pd.DataFrame(out).T
layers = {'protein': (prot, kg, 'sq'), 'rna': (rna, kg, 'sq'), 'magma': (magma, kg, 'lin'), 'metabolite': (metz, kc, 'sq')}
E = {}
for L, (sc, sets, st) in layers.items():
    E[L] = enrich(sc, sets, st, nperm=1000); E[L]['q'] = bh(E[L].p.values.astype(float))
    print(L, 'features', len(sc), 'pathways tested', len(E[L]), 'p<0.05', int((E[L].p < 0.05).sum()), 'q<0.1', int((E[L].q < 0.1).sum()))
T = pd.concat({L: e[['n', 'z', 'p', 'q']] for L, e in E.items()}, axis=1)
T.columns = [f'{a}_{b}' for a, b in T.columns]; T['name'] = [pname.get(i, i) for i in T.index]
zc = [c for c in T.columns if c.endswith('_z')]
pc = [c for c in T.columns if c.endswith('_p')]
T['k_layers'] = T[pc].notna().sum(1)
T['n_layers_p05'] = (T[pc] < 0.05).sum(1)
T['fisher_p'] = [stats.combine_pvalues(r.dropna().astype(float).values, 'fisher')[1] if r.notna().sum() >= 2 else np.nan for _, r in T[pc].iterrows()]
T['fisher_q'] = bh(T.fisher_p.values)
T.sort_values('fisher_p').to_csv(R + 'convergence_kegg.csv')
print(T.sort_values('fisher_p').head(30)[['name'] + zc + ['n_layers_p05', 'fisher_p', 'fisher_q']].round(3).to_string())
# ---------------- convergence null: permute gene labels of one layer
def top_frac(zs, f=0.1):
    return zs >= zs.quantile(1 - f)
obs_rows = []; null_rows = []
Lnames = list(layers)
obsZ = T[zc].copy(); obsZ.columns = Lnames
def conv_stats(Z):
    r = {}
    for i, a in enumerate(Lnames):
        for b in Lnames[i + 1:]:
            ok = Z[a].notna() & Z[b].notna()
            r[f'rho_{a}_{b}'] = stats.spearmanr(Z.loc[ok, a], Z.loc[ok, b]).correlation
    top = pd.DataFrame({L: top_frac(Z[L].dropna()).reindex(Z.index) for L in Lnames})
    r['n_top10_2plus'] = int((top.fillna(False).sum(1) >= 2).sum()); r['n_top10_3plus'] = int((top.fillna(False).sum(1) >= 3).sum())
    return r
obs = conv_stats(obsZ); print('observed', {k: round(v, 3) for k, v in obs.items()})
nulls = []
for it in range(NPERM):
    Z = obsZ.copy()
    for L in Lnames:  # permute every layer's gene labels independently -> no convergence, pathway overlap kept
        sc, sets, st = layers[L]; sp = pd.Series(rng.permutation(sc.values), index=sc.index)
        e = enrich_fast = None
        v = sp.values ** 2 if st == 'sq' else sp.values; idx = {g: i for i, g in enumerate(sp.index)}
        mu, sd = v.mean(), v.std()
        zz = {}
        for s, genes in sets.items():
            ii = [idx[g] for g in genes if g in idx]
            if len(ii) < 5: continue
            zz[s] = (v[ii].mean() - mu) / (sd / np.sqrt(len(ii)))
        Z[L] = pd.Series(zz).reindex(Z.index)
    nulls.append(conv_stats(Z))
N = pd.DataFrame(nulls)
# observed stats recomputed with the same analytic z (to be comparable with the null)
Zo = obsZ.copy()
for L in Lnames:
    sc, sets, st = layers[L]; v = sc.values ** 2 if st == 'sq' else sc.values; idx = {g: i for i, g in enumerate(sc.index)}; mu, sd = v.mean(), v.std(); zz = {}
    for s, genes in sets.items():
        ii = [idx[g] for g in genes if g in idx]
        if len(ii) >= 5: zz[s] = (v[ii].mean() - mu) / (sd / np.sqrt(len(ii)))
    Zo[L] = pd.Series(zz).reindex(Zo.index)
obs2 = conv_stats(Zo)
res = pd.DataFrame({'observed': obs2, 'null_mean': N.mean(), 'null_95': N.quantile(0.95),
                    'p_perm': {k: (np.sum(N[k] >= obs2[k]) + 1) / (len(N) + 1) for k in obs2}})
res.to_csv(R + 'convergence_null_test.csv'); print(res.round(3).to_string())
