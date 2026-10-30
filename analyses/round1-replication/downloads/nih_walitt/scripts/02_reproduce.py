"""Reproduce key Walitt 2024 claims from person-level source data, with exact/permutation tests and BH-FDR.
Outputs results/*.csv and prints a summary."""
import pandas as pd, numpy as np, warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
B = '/mnt/project-files/mecfs-research/downloads/nih_walitt'
T = f'{B}/results/tidy'; R = f'{B}/results'
rng = np.random.default_rng(1)
NPERM = 20000

def hedges_g(a, b):
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    d = (np.mean(b) - np.mean(a)) / sp if sp > 0 else np.nan
    return d * (1 - 3 / (4 * (na + nb) - 9))

def perm_p(a, b, stat=lambda a, b: np.mean(b) - np.mean(a)):
    r = stats.permutation_test((a, b), stat, n_resamples=NPERM, alternative='two-sided', random_state=rng, vectorized=False)
    return r.pvalue

def compare(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    mw = stats.mannwhitneyu(a, b, alternative='two-sided', method='exact' if len(a) + len(b) < 60 else 'auto')
    pos = np.all(a > 0) and np.all(b > 0)
    return dict(n_HV=len(a), n_PI=len(b), median_HV=np.median(a), median_PI=np.median(b),
                log2_ratio_means=np.log2(np.mean(b) / np.mean(a)) if pos else np.nan, hedges_g=hedges_g(a, b),
                p_mwu=mw.pvalue, p_perm_meandiff=perm_p(a, b), p_welch=stats.ttest_ind(a, b, equal_var=False).pvalue)

# ---------- A. scalar measures ----------
paper_p = {'CSF_DOPA': 0.021, 'CSF_DOPAC': 0.025, 'CSF_DHPG': 0.006, 'Blood_B_naive_pct': 0.03, 'Blood_B_switched_memory_pct': 0.008,
           'CSF_CD8_PD1_pct': 0.033, 'Blood_CD8_CD226_pct': 0.055}
s = pd.read_csv(f'{T}/scalar_measures_long.csv')
res = []
for m, d in s.groupby('measure'):
    r = compare(d[d.group == 'HV'].value, d[d.group == 'PI-ME/CFS'].value); r['measure'] = m; r['paper_p'] = paper_p.get(m, np.nan); res.append(r)
A = pd.DataFrame(res).set_index('measure')
A['q_BH_mwu'] = multipletests(A.p_mwu, method='fdr_bh')[1]
A.sort_values('p_mwu').to_csv(f'{R}/A_scalar_group_tests.csv')
print('\n=== A. scalar measures (sorted by MWU p) ===')
print(A.sort_values('p_mwu')[['n_HV', 'n_PI', 'median_HV', 'median_PI', 'hedges_g', 'p_mwu', 'p_perm_meandiff', 'q_BH_mwu', 'paper_p']].round(4).to_string())

# ---------- B. catechol / behaviour correlations ----------
p = pd.read_csv(f'{T}/paired_measures_long.csv')
res = []
for (src, g), d in p.groupby(['source', 'group']):
    if len(d) < 5: continue
    x, y = d.x.values, d.y.values
    pr = stats.pearsonr(x, y); sr = stats.spearmanr(x, y)
    # permutation p for Pearson
    perm = np.array([stats.pearsonr(x, rng.permutation(y))[0] for _ in range(5000)])
    pp = (np.sum(np.abs(perm) >= abs(pr[0])) + 1) / (len(perm) + 1)
    # leave-one-out: worst-case Pearson p
    loo = max(stats.pearsonr(np.delete(x, i), np.delete(y, i))[1] for i in range(len(x)))
    res.append(dict(source=src, group=g, x=d.x_name.iloc[0], y=d.y_name.iloc[0], n=len(d), pearson_r=pr[0], p_pearson=pr[1],
                    p_perm=pp, spearman_rho=sr[0], p_spearman=sr[1], p_pearson_worst_leave1out=loo))
Bt = pd.DataFrame(res)
Bt['q_BH_spearman_all_tests'] = multipletests(Bt.p_spearman, method='fdr_bh')[1]
Bt.to_csv(f'{R}/B_catechol_behaviour_correlations.csv', index=False)
print('\n=== B. correlations (paper Fig 6D-G, S11, S12) ===')
print(Bt.round(3).to_string())

# ---------- C. tilt plasma catecholamines ----------
t = pd.read_csv(f'{T}/tilt_plasma_catecholamines_wide.csv')
mins = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
def last_time(r, pre):
    vals = [r[f'{pre}.{m}'] for m in mins]
    ok = [m for m, v in zip(mins, vals) if pd.notna(v)]
    return max(ok) if ok else np.nan
t['last_sample_min'] = t.apply(lambda r: last_time(r, 'NE'), axis=1)
ne_cols = [f'NE.{m}' for m in mins]; ep_cols = [f'EPI.{m}' for m in mins]
t['NE_peak_upright'] = t[ne_cols].max(axis=1); t['EPI_peak_upright'] = t[ep_cols].max(axis=1)
t['NE_upright_mean'] = t[ne_cols].mean(axis=1); t['EPI_upright_mean'] = t[ep_cols].mean(axis=1)
t['NE_delta_8min'] = t['NE.8'] - t['NE.SUPINE']; t['NE_ratio_8min'] = t['NE.8'] / t['NE.SUPINE']
t['EPI_delta_8min'] = t['EPI.8'] - t['EPI.SUPINE']
t['NE_delta_peak'] = t['NE_peak_upright'] - t['NE.SUPINE']
t['completed_40'] = t['NE.40'].notna() | t['EPI.40'].notna()
t['hyperadrenergic_NE_ge_3.55nM'] = t['NE_peak_upright'] >= 3.55  # 600 pg/mL, if units are nmol/L
t.to_csv(f'{R}/C_tilt_person_level_derived.csv', index=False)
res = []
for v in ['NE.SUPINE', 'EPI.SUPINE', 'NE.4', 'NE.8', 'EPI.4', 'EPI.8', 'NE_delta_8min', 'NE_ratio_8min', 'EPI_delta_8min', 'NE_peak_upright', 'EPI_peak_upright', 'NE_delta_peak', 'NE_upright_mean', 'EPI_upright_mean', 'last_sample_min']:
    a = t.loc[t.group == 'HV', v].dropna(); b = t.loc[t.group == 'PI-ME/CFS', v].dropna()
    r = compare(a, b); r['variable'] = v; res.append(r)
C = pd.DataFrame(res).set_index('variable'); C['q_BH_mwu'] = multipletests(C.p_mwu, method='fdr_bh')[1]
C.to_csv(f'{R}/C_tilt_group_tests.csv')
print('\n=== C. tilt plasma NE/EPI (Supp Fig 4 source) ===')
print(C[['n_HV', 'n_PI', 'median_HV', 'median_PI', 'hedges_g', 'p_mwu', 'p_perm_meandiff', 'q_BH_mwu']].round(4).to_string())
ct = pd.crosstab(t.group, t.completed_40); print(ct); print('Fisher completed-40 p =', stats.fisher_exact(ct.values)[1])
ct2 = pd.crosstab(t.group, t['hyperadrenergic_NE_ge_3.55nM']); print(ct2); print('Fisher NE>=600pg/mL p =', stats.fisher_exact(ct2.values)[1])
# does NE response relate to tilt tolerance (time to end of sampling) within patients?
for g in ['PI-ME/CFS', 'HV']:
    d = t[t.group == g].dropna(subset=['NE_delta_8min', 'last_sample_min'])
    print(g, 'Spearman NE delta 8 min vs last sample time: rho=%.2f p=%.3f n=%d' % (*stats.spearmanr(d.NE_delta_8min, d.last_sample_min), len(d)))
    d = t[t.group == g].dropna(subset=['EPI_delta_8min', 'last_sample_min'])
    print(g, 'Spearman EPI delta 8 min vs last sample time: rho=%.2f p=%.3f n=%d' % (*stats.spearmanr(d.EPI_delta_8min, d.last_sample_min), len(d)))

# ---------- D. flow cytometry ----------
f = pd.read_csv(f'{T}/flow_cytometry_long.csv')
res = []
for (ti, pop), d in f.groupby(['tissue', 'population']):
    a = d[d.group == 'HV'].value; b = d[d.group == 'PI-ME/CFS'].value
    if len(a) < 5 or len(b) < 5: continue
    r = compare(a, b); r.update(tissue=ti, population=pop); res.append(r)
D = pd.DataFrame(res)
for ti in D.tissue.unique():
    D.loc[D.tissue == ti, 'q_BH_mwu_within_tissue'] = multipletests(D.loc[D.tissue == ti, 'p_mwu'], method='fdr_bh')[1]
D.sort_values('p_mwu').to_csv(f'{R}/D_flow_cytometry_tests.csv', index=False)
print('\n=== D. flow cytometry, top 12 ===')
print(D.sort_values('p_mwu').head(12)[['tissue', 'population', 'n_HV', 'n_PI', 'median_HV', 'median_PI', 'hedges_g', 'p_mwu', 'p_perm_meandiff', 'q_BH_mwu_within_tissue']].round(4).to_string())

# ---------- E. EEfRT ----------
e = pd.read_csv(f'{T}/eefrt_trials.csv')
e = e[(e.Valid_is_1 == 1) & (e.Adjudicated_is_1 == 1) & (e.Trial >= 1)]
e_choice = e[e.Choice_Time != 5]  # SAS script: Choice_Time ^= 5 (computer-made choices excluded)
per = e_choice.groupby('mapID').agg(group=('Status_MECFS_is_1', 'first'), PHTC=('Trial_Difficulty_Hard_is_1', 'mean'), n_trials=('Trial', 'size'))
hard = e[e.Trial_Difficulty_Hard_is_1 == 1].groupby('mapID').Successful_Completion_Yes_is_1.mean().rename('hard_completion_rate')
easy = e[e.Trial_Difficulty_Hard_is_1 == 0].groupby('mapID').Successful_Completion_Yes_is_1.mean().rename('easy_completion_rate')
rate = e[e.Trial_Difficulty_Hard_is_1 == 1].groupby('mapID').Button_Press_Rate.mean().rename('hard_press_rate')
per = per.join(hard).join(easy).join(rate)
per.to_csv(f'{R}/E_eefrt_person_level.csv')
print('\n=== E. EEfRT person level ===')
for v in ['PHTC', 'hard_completion_rate', 'easy_completion_rate', 'hard_press_rate']:
    a = per.loc[per.group == 0, v].dropna(); b = per.loc[per.group == 1, v].dropna()
    r = compare(a, b); print(v, {k: round(x, 4) if isinstance(x, float) else x for k, x in r.items()})
print(per.groupby('group')[['PHTC', 'hard_completion_rate']].describe().round(3).T.to_string())
d1 = per[per.group == 1].dropna()
print('PI: Spearman PHTC vs hard-completion rate', stats.spearmanr(d1.PHTC, d1.hard_completion_rate))
# GEE logistic for hard choice, exchangeable within person
import statsmodels.api as sm, statsmodels.formula.api as smf
ee = e_choice.copy(); ee['ME'] = ee.Status_MECFS_is_1; ee['male'] = ee.Sex_Male_is_1; ee['EV'] = ee.Value_of_Reward * ee.Probability_of_Reward; ee['Trial'] = ee.Trial - 1
mod = smf.gee('Trial_Difficulty_Hard_is_1 ~ ME + Value_of_Reward + Probability_of_Reward + EV + Trial + male', 'mapID', ee, family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
print(mod.summary().tables[1])
print('OR(ME) = %.3f' % np.exp(mod.params['ME']))
# sensitivity: add each person's hard completion rate (ability) as covariate
ee = ee.merge(per[['hard_completion_rate']], left_on='mapID', right_index=True)
mod2 = smf.gee('Trial_Difficulty_Hard_is_1 ~ ME + Value_of_Reward + Probability_of_Reward + EV + Trial + male + hard_completion_rate', 'mapID', ee, family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
print('adjusted for hard completion rate: OR(ME)=%.3f p=%.4f' % (np.exp(mod2.params['ME']), mod2.pvalues['ME']))
