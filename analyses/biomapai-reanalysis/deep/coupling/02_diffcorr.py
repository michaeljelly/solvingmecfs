"""Task 2: differential correlation (controls vs patients) of known microbe-derived plasma metabolites with their
producer taxa / gene functions. Unit = person (mean over visits that have both stool and plasma), covariates
(age, sex, BMI) regressed out within each group, Spearman r, Fisher z test (var 1.06/(n-3)), BH-FDR.
Pre-specified pairs are one FDR family; a data-driven screen (target metabolites x all 118 species CLR + gene sets
+ Shannon) is a second family."""
from cc import *
from scipy import stats

X, ann = metab(0.5)
spc, _ = species_clr()
full, tax = species_full()
klog, kraw = kegg_logcpm()

TARGETS = {
 'indolepropionate': ['gs_stickland', 'shannon'],
 'indolelactate': ['tx_Bifidobacterium', 'gs_stickland'],
 'indoleacetate': ['tx_Bacteroides', 'gs_iorA'],
 '3-indoxyl sulfate': ['gs_tnaA', 'tx_Ecoli', 'tx_Bacteroides'],
 '6-hydroxyindole sulfate': ['gs_tnaA', 'tx_Ecoli'],
 'p-cresol sulfate': ['gs_hpdB', 'tx_pcresol'],
 'p-cresol glucuronide*': ['gs_hpdB', 'tx_pcresol'],
 'phenol sulfate': ['gs_tpl'],
 'phenylacetylglutamine': ['gs_stickland', 'shannon'],
 '4-hydroxyphenylacetylglutamine': ['gs_stickland'],
 'hippurate': ['gs_hca_mhp', 'shannon'],
 'benzoate': ['gs_hca_mhp', 'shannon'],
 '3-phenylpropionate (hydrocinnamate)': ['gs_hca_mhp', 'shannon'],
 'cinnamoylglycine': ['gs_hca_mhp', 'shannon'],
 'deoxycholate': ['gs_bai', 'tx_Cscindens', 'gs_bsh'],
 'glycodeoxycholate 3-sulfate': ['gs_bai', 'tx_Cscindens', 'gs_bsh'],
 'taurodeoxycholate': ['gs_bai', 'tx_Cscindens', 'gs_bsh'],
 'glycolithocholate': ['gs_bai', 'tx_Cscindens'],
 'lithocholate sulfate (1)': ['gs_bai', 'tx_Cscindens'],
 'glycolithocholate sulfate*': ['gs_bai', 'tx_Cscindens'],
 'taurolithocholate 3-sulfate': ['gs_bai', 'tx_Cscindens'],
 'isoursodeoxycholate': ['tx_Caerofaciens', 'tx_Rgnavus', 'gs_hdhA'],
 'ursodeoxycholate': ['tx_Caerofaciens', 'tx_Rgnavus', 'gs_hdhA'],
 'glycoursodeoxycholate': ['tx_Caerofaciens', 'tx_Rgnavus'],
 'imidazole propionate': ['gs_urdA', 'gs_hutH', 'tx_Elenta'],
 'trimethylamine N-oxide': ['gs_cutC', 'gs_cntA', 'tx_TMA'],
 'butyrate/isobutyrate (4:0)': ['gs_butyrate', 'tx_butyrate'],
 'isovalerate (i5:0)': ['gs_bcfa', 'shannon'],
 'isobutyrylglycine': ['gs_bcfa', 'shannon'],
 'isovalerylglycine': ['gs_bcfa', 'shannon'],
 '4-ethylphenylsulfate': ['shannon'],
 '3-formylindole': ['gs_tnaA', 'shannon'],
}
GENES = {'gs_stickland': ['K04090', 'K00170', 'K00171', 'K20903', 'K20904'], 'gs_iorA': ['K04090'],
         'gs_tnaA': ['K01667'], 'gs_hpdB': ['K18427'], 'gs_tpl': ['K01668'],
         'gs_hca_mhp': ['K05708', 'K00529', 'K05712', 'K05713', 'K02554'],
         'gs_bai': ['K15868', 'K15870', 'K15871', 'K15873', 'K22605'], 'gs_bsh': ['K01442'], 'gs_hdhA': ['K00076'],
         'gs_urdA': ['K17363'], 'gs_hutH': ['K01745'], 'gs_cutC': ['K20038', 'K20037'], 'gs_cntA': ['K22443', 'K22444'],
         'gs_butyrate': ['K00929', 'K00634', 'K01034', 'K01035'], 'gs_bcfa': ['K00187', 'K03738']}
TAXA = {'tx_Bifidobacterium': r'^Bifidobacterium_', 'tx_Bacteroides': r'^(Bacteroides|Phocaeicola)_',
        'tx_Ecoli': r'^Escherichia_coli$', 'tx_pcresol': r'^(Blautia_hydrogenotrophica|Romboutsia_timonensis|Intestinibacter_bartlettii|Clostridioides_difficile)$',
        'tx_Cscindens': r'^Clostridium_scindens$', 'tx_Caerofaciens': r'^Collinsella_aerofaciens$',
        'tx_Rgnavus': r'^Ruminococcus_gnavus$', 'tx_Elenta': r'^Eggerthella_lenta$',
        'tx_TMA': r'^(Hungatella_hathewayi|Emergencia_timonensis|Escherichia_coli|Desulfovibrio_.*)$',
        'tx_butyrate': r'^(Faecalibacterium_prausnitzii|Eubacterium_rectale|Roseburia_.*|Anaerobutyricum_hallii|Anaerostipes_hadrus|Coprococcus_.*)$'}

P = pd.DataFrame(index=full.index)
for k, rx in TAXA.items():
    cols = [c for c in full.columns if pd.Series([c]).str.contains(rx).iloc[0]]
    P[k] = np.log10(full[cols].sum(axis=1) + 1e-5)
p = full.where(full > 0)
P['shannon'] = -(p * np.log(p)).sum(axis=1)
cpm = kraw.div(kraw.sum(axis=1), axis=0) * 1e6
for k, ks in GENES.items():
    ks = [x for x in ks if x in cpm.columns]
    P[k] = np.log2(cpm[ks].sum(axis=1) + 1)

common = X.index.intersection(P.index).intersection(spc.index).intersection(kraw.index).intersection(meta.index)
M = meta.loc[common]
tcols = [t for t in TARGETS if t in X.columns]
print('missing targets', [t for t in TARGETS if t not in X.columns])
D = pd.concat([X.loc[common, tcols], P.loc[common], spc.loc[common].add_prefix('sp_')], axis=1)
D['subject'] = M.subject
per = D.groupby('subject').mean()
pm = M.groupby('subject')[COV + ['y']].first()
pm['illness_duration'] = M.groupby('subject').illness_duration.first()
print('persons', pm.y.value_counts().to_dict())


def grp_resid(g):
    idx = pm.index[g]
    return resid(per.loc[idx], pm.loc[idx, COV])


Rc, Rp = grp_resid(pm.y == 0), grp_resid(pm.y == 1)
Rs = grp_resid(pm.illness_duration == 'Short')


def sp_r(A, a, b):
    return stats.spearmanr(A[a], A[b])[0]


def compare(pairs):
    rows = []
    nc, np_, ns = len(Rc), len(Rp), len(Rs)
    for m, q in pairs:
        rc, rp, rs = sp_r(Rc, m, q), sp_r(Rp, m, q), sp_r(Rs, m, q)
        z = (np.arctanh(rc) - np.arctanh(rp)) / np.sqrt(1.06 / (nc - 3) + 1.06 / (np_ - 3))
        zs = (np.arctanh(rc) - np.arctanh(rs)) / np.sqrt(1.06 / (nc - 3) + 1.06 / (ns - 3))
        pc = 2 * stats.norm.sf(abs(np.arctanh(rc)) / np.sqrt(1.06 / (nc - 3)))
        pp = 2 * stats.norm.sf(abs(np.arctanh(rp)) / np.sqrt(1.06 / (np_ - 3)))
        rows.append(dict(metabolite=m, predictor=q, r_ctrl=rc, p_ctrl=pc, r_pat=rp, p_pat=pp, r_short=rs,
                         z_diff=z, p_diff=2 * stats.norm.sf(abs(z)), p_diff_short=2 * stats.norm.sf(abs(zs)),
                         n_ctrl=nc, n_pat=np_, n_short=ns))
    O = pd.DataFrame(rows)
    for c in ['p_ctrl', 'p_pat', 'p_diff', 'p_diff_short']:
        O[c.replace('p_', 'q_')] = bh(O[c])
    return O


pre = compare([(m, q) for m in tcols for q in TARGETS[m]])
pre.to_csv(OUT + 'diffcorr_prespecified.csv', index=False)
pd.set_option('display.width', 250)
print(pre[['metabolite', 'predictor', 'r_ctrl', 'r_pat', 'r_short', 'z_diff', 'p_diff', 'q_diff', 'q_ctrl', 'q_pat']].round(3).to_string())

# data-driven screen
preds = [c for c in per.columns if c.startswith('sp_')] + list(GENES) + list(TAXA) + ['shannon']
scr = compare([(m, q) for m in tcols for q in preds])
scr.to_csv(OUT + 'diffcorr_screen.csv', index=False)
print('\nscreen pairs', len(scr), 'q_diff<0.05:', (scr.q_diff < .05).sum(), 'q_diff<0.1:', (scr.q_diff < .1).sum())
print('pairs significant (q<0.05) in ctrl', (scr.q_ctrl < .05).sum(), 'in pat', (scr.q_pat < .05).sum())
print(scr.sort_values('p_diff').head(15)[['metabolite', 'predictor', 'r_ctrl', 'r_pat', 'z_diff', 'p_diff', 'q_diff']].round(3).to_string())
# concordance of correlation patterns
print('overall correlation of r_ctrl vs r_pat across screen pairs: %.3f' % np.corrcoef(scr.r_ctrl, scr.r_pat)[0, 1])
sig = scr[(scr.q_ctrl < .05) | (scr.q_pat < .05)]
print('among pairs significant in either group (n=%d): same sign %.2f, mean |r| ctrl %.3f pat %.3f' % (
    len(sig), (np.sign(sig.r_ctrl) == np.sign(sig.r_pat)).mean(), sig.r_ctrl.abs().mean(), sig.r_pat.abs().mean()))
# permutation of group labels (persons) for the number of q_diff<0.1 and the mean |r_ctrl|-|r_pat| gap
rng = np.random.default_rng(7)
allidx = pm.index
obs_gap = (scr.r_ctrl.abs() - scr.r_pat.abs()).mean()
gaps = []
Rall = per.loc[allidx]
for it in range(200):
    lab = rng.permutation(pm.y.values)
    a = resid(Rall[lab == 0], pm.loc[allidx[lab == 0], COV]); b = resid(Rall[lab == 1], pm.loc[allidx[lab == 1], COV])
    ra = a[tcols].rank().T.values; rb = b[tcols].rank().T.values
    pa = a[preds].rank().values; pb = b[preds].rank().values
    ca = np.corrcoef(np.vstack([ra, pa.T]))[:len(tcols), len(tcols):]
    cb = np.corrcoef(np.vstack([rb, pb.T]))[:len(tcols), len(tcols):]
    gaps.append(np.nanmean(np.abs(ca)) - np.nanmean(np.abs(cb)))
gaps = np.array(gaps)
print('mean |r| ctrl - |r| pat over screen: %.4f ; label-permutation null mean %.4f sd %.4f ; p=%.3f' % (
    obs_gap, gaps.mean(), gaps.std(), (1 + (np.abs(gaps - gaps.mean()) >= abs(obs_gap - gaps.mean())).sum()) / 201))
pd.Series(dict(obs_gap=obs_gap, null_mean=gaps.mean(), null_sd=gaps.std())).to_csv(OUT + 'diffcorr_global_perm.csv')
