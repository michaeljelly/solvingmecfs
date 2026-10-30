"""CSF metabolomics re-analysis (Walitt 2024), plasma lipidomics, stool NMR: reproduce and test project leads.
Leads: indoles/IPA, kynurenine pathway, fatty acids/ketones/acylcarnitines, bile acids, catecholamine metabolites (VMA/HVA)."""
import pandas as pd, numpy as np, warnings, re
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
B = '/mnt/project-files/mecfs-research/downloads/nih_walitt'
T = f'{B}/results/tidy'; R = f'{B}/results'

m = pd.read_csv(f'{T}/csf_metabolomics_wide.csv')
meta = m[['group', 'sample_id', 'participant_no_if_female', 'female']]
X = m.iloc[:, 4:].astype(float)
print('CSF metabolomics:', X.shape, pd.crosstab(meta.group, meta.female).to_dict())
L = np.log2(X)
hv = (meta.group == 'HV').values; pi = ~hv

def test_block(L, mask_a, mask_b, label):
    rows = []
    for c in L.columns:
        a = L.loc[mask_a, c].values; b = L.loc[mask_b, c].values
        mw = stats.mannwhitneyu(a, b, method='exact' if (len(a) + len(b)) < 40 else 'asymptotic')
        sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
        rows.append(dict(metabolite=c, log2FC=b.mean() - a.mean(), d=(b.mean() - a.mean()) / sp if sp > 0 else np.nan, p=mw.pvalue))
    d = pd.DataFrame(rows); d['q'] = multipletests(d.p, method='fdr_bh')[1]; d['set'] = label
    return d

# 1. all participants, raw (as published, batch-normalised imputed)
res_all = test_block(L, hv, pi, 'all')
# 2. global shift check: per-sample median log2 abundance
samp_med = L.median(axis=1)
print('\nPer-sample median log2 level: HV %.3f, PI %.3f, MWU p=%.4g' % (samp_med[hv].median(), samp_med[pi].median(), stats.mannwhitneyu(samp_med[hv], samp_med[pi]).pvalue))
print('Median per-metabolite log2FC (PI-HV): %.3f; fraction of metabolites lower in PI: %.2f' % (res_all.log2FC.median(), (res_all.log2FC < 0).mean()))
# sign test on whether more metabolites go down than chance
print('Binomial test on direction:', stats.binomtest(int((res_all.log2FC < 0).sum()), len(res_all)).pvalue)
# 3. median-normalised (removes a global per-sample dilution/scale factor)
Ln = L.sub(samp_med, axis=0)
res_norm = test_block(Ln, hv, pi, 'all_mediannorm')
# 4. sex strata
res_f = test_block(L, hv & (meta.female == 1).values, pi & (meta.female == 1).values, 'female')
res_m = test_block(L, hv & (meta.female == 0).values, pi & (meta.female == 0).values, 'male')
allres = pd.concat([res_all, res_norm, res_f, res_m])
allres.to_csv(f'{R}/F_csf_metabolomics_all_tests.csv', index=False)
for s_, d in allres.groupby('set'):
    print(f'{s_}: n tests {len(d)}, p<0.05 {(d.p < 0.05).sum()}, q<0.05 {(d.q < 0.05).sum()}, q<0.10 {(d.q < 0.10).sum()}')

# paper's numbers: Fig 6I list and SD14A
sd14 = pd.read_excel(f'{B}/raw/nature_esm/suppdata/Supplementary Data 14.xlsx', header=None)
hdr = sd14.index[sd14[0] == 'CHEM_ID'][0]
end = sd14.index[(sd14.index > hdr) & sd14[0].isna()][0]
a = sd14.iloc[hdr + 1:end, :9]; a.columns = sd14.iloc[hdr, :9].tolist()
a = a.astype({'p-value': float, 'adjusted p-value': float, 'fold change': float})
a['BH_recomputed_from_listed_p'] = multipletests(a['p-value'], method='fdr_bh')[1]
print('\nSD14A rows:', len(a), '| listed adjusted p monotone increasing in row order:', bool(np.all(np.diff(a['adjusted p-value'].values) >= 0)),
      '| rows where listed adj p < listed p (impossible for BH):', int((a['adjusted p-value'] < a['p-value'] - 1e-12).sum()))
a.to_csv(f'{R}/F_SD14A_check.csv', index=False)
merged = a.merge(res_all, left_on='CHEMICAL_NAME', right_on='metabolite', how='left')
print('SD14A metabolites with our BH q<0.05 (all 445 tested):', int((merged.q < 0.05).sum()), 'of', len(merged))

# 5. lead panel
lead_patterns = {
 'tryptophan/indole': r'^tryptophan$|kynuren|anthranil|picolin|quinolin|5-hydroxyindole|indole|indoxyl|serotonin|tryptophan betaine|C-glycosyltryptophan|xanthuren',
 'catecholamine/tyrosine': r'^tyrosine$|^phenylalanine$|homovanill|methoxytyr|dopamine|methoxytyramine|vanillyl|4-hydroxyphenylpyruvate|3-\(4-hydroxyphenyl\)lactate|catechol sulfate|^DOPA',
 'ketones/FA oxidation': r'hydroxybutyrate \(BHBA\)|acetoacetate|carnitine|hydroxyoctanoate|hydroxydecanoate|^\S*noate \(|glycerol$',
 'bile acids/sterols': r'cholate|chenodeoxy|deoxychol|7-Hoca|cholesterol$',
 'microbial co-metabolites': r'p-cresol|phenol sulfate|hippurate|phenylacetate|butyrate/isobutyrate|trimethylamine|imidazole propionate|phenyllactate|3-phenylpropionate',
}
rows = []
for fam, pat in lead_patterns.items():
    for c in L.columns:
        if re.search(pat, c, re.I):
            r = {'family': fam, 'metabolite': c}
            for s_ in ['all', 'all_mediannorm', 'female', 'male']:
                q = allres[(allres.set == s_) & (allres.metabolite == c)].iloc[0]
                r[f'log2FC_{s_}'] = q.log2FC; r[f'p_{s_}'] = q.p; r[f'q_{s_}'] = q.q
            rows.append(r)
lead = pd.DataFrame(rows).drop_duplicates('metabolite')
lead.to_csv(f'{R}/F_csf_lead_metabolites.csv', index=False)
pd.set_option('display.width', 250)
print('\n=== lead metabolites in CSF (log2FC PI vs HV; q = BH over all 445) ===')
print(lead[['family', 'metabolite', 'log2FC_all', 'p_all', 'q_all', 'log2FC_all_mediannorm', 'p_all_mediannorm', 'log2FC_male', 'p_male']].round(3).to_string())

# 6. kynurenine/tryptophan ratio and indole sum
kt = L['kynurenine'] - L['tryptophan']
print('\nCSF log2(kyn/trp): HV %.2f PI %.2f MWU p=%.3f' % (kt[hv].mean(), kt[pi].mean(), stats.mannwhitneyu(kt[hv], kt[pi]).pvalue))
ka = L['kynurenate'] - L['kynurenine']
print('CSF log2(kynurenate/kyn): HV %.2f PI %.2f MWU p=%.3f' % (ka[hv].mean(), ka[pi].mean(), stats.mannwhitneyu(ka[hv], ka[pi]).pvalue))

# 7. what survives median normalisation
print('\nFDR<0.10 after median normalisation:'); print(res_norm[res_norm.q < 0.10].sort_values('p').round(4).to_string())
# sample order / batch check: sample ids vs group
m_ = meta.copy(); m_['num'] = m_.sample_id.str.extract(r'(\d+)').astype(int); m_['median_log2'] = samp_med.values
print(m_.sort_values('num')[['sample_id', 'group', 'female', 'median_log2']].to_string())
print('Spearman sample-number vs median level within PI:', stats.spearmanr(m_[pi].num, m_[pi].median_log2), ' within HV:', stats.spearmanr(m_[hv].num, m_[hv].median_log2))
# 8. plasma lipid classes (SD20) and species (Fig S19A)
lc = pd.read_csv(f'{T}/plasma_lipid_classes_long.csv')
rows = []
for f_, d in lc.groupby('feature'):
    a_ = d[d.group == 'HV'].value; b_ = d[d.group == 'PI-ME/CFS'].value
    rows.append(dict(lipid_class=f_, n_HV=len(a_), n_PI=len(b_), median_HV=a_.median(), median_PI=b_.median(), log2_ratio=np.log2(b_.mean() / a_.mean()), p_mwu=stats.mannwhitneyu(a_, b_).pvalue))
lcr = pd.DataFrame(rows); lcr['q'] = multipletests(lcr.p_mwu, method='fdr_bh')[1]
lcr.sort_values('p_mwu').to_csv(f'{R}/G_plasma_lipid_classes.csv', index=False)
print('\n=== plasma lipid classes ==='); print(lcr.sort_values('p_mwu').round(3).to_string())
ls = pd.read_csv(f'{T}/plasma_lipidomics_species_wide.csv')
spc = [c for c in ls.columns if c not in ('Group', 'Birth Sex', 'group')]
S = ls[spc].apply(pd.to_numeric, errors='coerce')
rows = []
for c in spc:
    a_ = S.loc[ls.group == 'HV', c].dropna(); b_ = S.loc[ls.group == 'PI-ME/CFS', c].dropna()
    if len(a_) < 5 or len(b_) < 5 or (a_.sum() + b_.sum()) == 0: continue
    rows.append(dict(species=c, lipid_class=re.match(r'[A-Za-z]+', c).group(0), log2_ratio=np.log2((b_.mean() + 1e-9) / (a_.mean() + 1e-9)), p=stats.mannwhitneyu(a_, b_).pvalue))
lsr = pd.DataFrame(rows); lsr['q'] = multipletests(lsr.p, method='fdr_bh')[1]
lsr.sort_values('p').to_csv(f'{R}/G_plasma_lipid_species.csv', index=False)
print('\nlipid species tested %d, p<0.05: %d, q<0.05: %d' % (len(lsr), (lsr.p < 0.05).sum(), (lsr.q < 0.05).sum()))
print(lsr[lsr.lipid_class == 'FFA'].sort_values('p').round(3).to_string())
print(lsr.groupby('lipid_class').agg(n=('p','size'), median_log2_ratio=('log2_ratio','median'), n_p05=('p', lambda x: (x<0.05).sum())).round(3).to_string())
print(lsr.sort_values('p').head(10).round(4).to_string())

# 9. stool NMR (SD21)
st = pd.read_csv(f'{T}/stool_nmr_metabolites_long.csv')
rows = []
for f_, d in st.groupby('feature'):
    a_ = d[d.group == 'HV'].value; b_ = d[d.group == 'PI-ME/CFS'].value
    rows.append(dict(metabolite=f_, n_HV=len(a_), n_PI=len(b_), median_HV=a_.median(), median_PI=b_.median(), p_mwu=stats.mannwhitneyu(a_, b_).pvalue))
str_ = pd.DataFrame(rows); str_['q'] = multipletests(str_.p_mwu, method='fdr_bh')[1]
str_.sort_values('p_mwu').to_csv(f'{R}/H_stool_nmr.csv', index=False)
print('\n=== stool NMR ==='); print(str_.sort_values('p_mwu').round(4).to_string())
