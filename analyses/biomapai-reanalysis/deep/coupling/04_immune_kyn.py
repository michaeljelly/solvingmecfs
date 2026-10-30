"""Task 4: IFN-g / granzyme T-cell activity (MAIT, gamma-delta, CD8, NK-like) vs tryptophan-kynurenine metabolites.
Between-person: person means (visits with both immune + plasma), covariates regressed out within group, Spearman,
controls vs patients Fisher z. Within-person (patients with >=2 paired visits): mixed model
z(metabolite) ~ immune_between + immune_within + covariates + timepoint + (1|person). BH-FDR per family."""
from cc import *
from scipy import stats
import statsmodels.formula.api as smf

X, ann = metab(0.5)
IM = immune()
imf = [c for c in IM.columns if (('MAIT' in c or 'CD8' in c or 'DN ' in c or 'CD3- CCR6-' in c) and
                                 any(k in c for k in ['INFg+', 'GzA+', 'Granzyme', 'TNFa+'])
                                 and 'IL-17' not in c and 'IL4' not in c and 'GMCSF' not in c)]
imf = [c for c in imf if 'total' in c or 'GzA+' in c or 'Granzyme' in c] + ['Total gd Cells (gdTCR CD3hi + CD3lo)', 'CD3+ gdTCR+ CD3hi', 'CD3+ gdTCR+ CD3lo', '% CD8+ MAIT d1', '%CD8+ MAIT d0']
imf = list(dict.fromkeys([c for c in imf if c in IM.columns]))
mets = [m for m in ['tryptophan', 'kynurenine', 'kynurenate', 'quinolinate', 'xanthurenate', 'picolinate', 'anthranilate',
                    '3-hydroxykynurenine', 'serotonin', 'indolepropionate', 'indolelactate', '3-indoxyl sulfate', 'indoleacetate'] if m in X.columns]
Mt = X[mets].copy(); Mt['KT_ratio'] = X.kynurenine - X.tryptophan
mets = list(Mt.columns)
print(len(imf), 'immune features;', mets)
common = IM.index.intersection(Mt.index).intersection(meta.index)
M = meta.loc[common]
D = pd.concat([IM.loc[common, imf], Mt.loc[common]], axis=1)
per = D.groupby(M.subject).mean()
pm = M.groupby('subject')[COV + ['y']].first()
Rc = resid(per[pm.y == 0], pm[pm.y == 0][COV]); Rp = resid(per[pm.y == 1], pm[pm.y == 1][COV])
nc, np_ = len(Rc), len(Rp)
rows = []
for i in imf:
    for m in mets:
        rc = stats.spearmanr(Rc[i], Rc[m])[0]; rp = stats.spearmanr(Rp[i], Rp[m])[0]
        z = (np.arctanh(rc) - np.arctanh(rp)) / np.sqrt(1.06 / (nc - 3) + 1.06 / (np_ - 3))
        rows.append(dict(immune=i, metabolite=m, r_ctrl=rc, r_pat=rp,
                         p_ctrl=2 * stats.norm.sf(abs(np.arctanh(rc)) * np.sqrt((nc - 3) / 1.06)),
                         p_pat=2 * stats.norm.sf(abs(np.arctanh(rp)) * np.sqrt((np_ - 3) / 1.06)),
                         z_diff=z, p_diff=2 * stats.norm.sf(abs(z))))
B = pd.DataFrame(rows)
for c in ['p_ctrl', 'p_pat', 'p_diff']:
    B[c.replace('p_', 'q_')] = bh(B[c])
B.to_csv(OUT + 'immune_kyn_between.csv', index=False)
pd.set_option('display.width', 230)
print('persons ctrl', nc, 'pat', np_, 'pairs', len(B))
print('q<.05 ctrl', (B.q_ctrl < .05).sum(), 'pat', (B.q_pat < .05).sum(), 'diff', (B.q_diff < .05).sum())
print(B.sort_values('p_pat').head(12).round(3).to_string())
print(B.sort_values('p_diff').head(8).round(3).to_string())
key = B[B.metabolite.isin(['kynurenine', 'KT_ratio', 'tryptophan'])]
print('kyn/KT/trp pairs: median r ctrl %.3f pat %.3f' % (key.r_ctrl.median(), key.r_pat.median()))

# within-person (patients)
DD = D.join(M[COV + ['y', 'subject', 'tp']])
P = DD[DD.y == 1].copy()
P = P[P.groupby('subject').subject.transform('size') >= 2]
print('within: patient visits', len(P), 'persons', P.subject.nunique())
sel = [i for i in imf if 'total INFg+' in i or 'GzA+' in i or 'Granzyme' in i or 'gd' in i]
rows = []
for i in sel:
    P['ib'] = P.groupby('subject')[i].transform('mean'); P['iw'] = P[i] - P['ib']
    sdv = P[i].std(); P['ib'] /= sdv; P['iw'] /= sdv
    for m in ['kynurenine', 'tryptophan', 'KT_ratio', 'kynurenate', 'quinolinate']:
        if m not in P: continue
        P['yv'] = (P[m] - P[m].mean()) / P[m].std()
        try:
            f = smf.mixedlm('yv ~ ib + iw + age + male + bmi_i + bmi_na + C(tp)', P, groups=P.subject).fit(reml=True)
            rows.append(dict(immune=i, metabolite=m, b_between=f.params['ib'], p_between=f.pvalues['ib'],
                             b_within=f.params['iw'], p_within=f.pvalues['iw']))
        except Exception as e:
            print('fail', i, m, e)
W = pd.DataFrame(rows); W['q_between'] = bh(W.p_between); W['q_within'] = bh(W.p_within)
W.to_csv(OUT + 'immune_kyn_within_patients.csv', index=False)
print('within tests', len(W), 'q_within<.05', (W.q_within < .05).sum(), 'q_between<.05', (W.q_between < .05).sum())
print(W.sort_values('p_within').head(8).round(4).to_string())
print(W.sort_values('p_between').head(8).round(4).to_string())
