"""Task 3: (a) are microbe-derived plasma metabolites different in patients? Linear mixed models on all visits:
z(log2 metabolite) ~ group + age + sex + BMI + timepoint + (1|person); also Short vs control, Long vs control.
(b) within patients: symptom ~ metabolite_between(person mean) + metabolite_within(visit deviation) + covariates + timepoint + (1|person).
BH-FDR within each family."""
from cc import *
import statsmodels.formula.api as smf
from scipy import stats

X, ann = metab(0.5)
MICRO = ['indolepropionate', 'indolelactate', 'indoleacetate', '3-indoxyl sulfate', '6-hydroxyindole sulfate', '3-formylindole',
         'indoleacetylglutamine', 'p-cresol sulfate', 'p-cresol glucuronide*', 'phenol sulfate', '4-ethylphenylsulfate',
         'phenylacetylglutamine', '4-hydroxyphenylacetylglutamine', 'hippurate', '3-hydroxyhippurate', '4-hydroxyhippurate',
         'benzoate', '3-phenylpropionate (hydrocinnamate)', 'cinnamoylglycine', 'catechol sulfate',
         'deoxycholate', 'glycodeoxycholate 3-sulfate', 'taurodeoxycholate', 'glycolithocholate', 'lithocholate sulfate (1)',
         'glycolithocholate sulfate*', 'taurolithocholate 3-sulfate', 'isoursodeoxycholate', 'ursodeoxycholate', 'glycoursodeoxycholate',
         'hyocholate', 'imidazole propionate', 'trimethylamine N-oxide', 'butyrate/isobutyrate (4:0)', 'isovalerate (i5:0)',
         'isobutyrylglycine', 'isovalerylglycine', 'equol sulfate', 'enterolactone sulfate',
         'tryptophan', 'kynurenine', 'kynurenate', 'serotonin', 'N,N,N-trimethyl-5-aminovalerate']
MICRO = [m for m in MICRO if m in X.columns]
Z = X[MICRO].copy()
Z['KT_ratio'] = X['kynurenine'] - X['tryptophan']
# composite: microbial indoles (IPA, ILA, IAA, indoxyl sulfate)
Z = (Z - Z.mean()) / Z.std()
Z['indole_composite'] = Z[['indolepropionate', 'indolelactate', 'indoleacetate', '3-indoxyl sulfate']].mean(axis=1)
Z['secBA_composite'] = Z[[c for c in ['deoxycholate', 'glycodeoxycholate 3-sulfate', 'taurodeoxycholate', 'glycolithocholate', 'lithocholate sulfate (1)',
                                     'glycolithocholate sulfate*', 'taurolithocholate 3-sulfate'] if c in Z]].mean(axis=1)
feats = list(Z.columns)
safe = {f: 'm%d' % i for i, f in enumerate(feats)}
D = meta.join(Z.rename(columns=safe), how='inner')
D['grp'] = D.illness_duration.map({'Control': 'Control', 'Short': 'Short', 'Long': 'Long'})
print('visits', len(D), D.y.value_counts().to_dict(), 'persons', D.groupby('y').subject.nunique().to_dict())

rows = []
for f in feats:
    v = safe[f]
    try:
        m = smf.mixedlm(f'{v} ~ y + age + male + bmi_i + bmi_na + C(tp)', D, groups=D.subject).fit(reml=True)
        ci = m.conf_int().loc['y']
        m2 = smf.mixedlm(f'{v} ~ C(grp, Treatment("Control")) + age + male + bmi_i + bmi_na + C(tp)', D, groups=D.subject).fit(reml=True)
        ks = 'C(grp, Treatment("Control"))[T.Short]'; kl = 'C(grp, Treatment("Control"))[T.Long]'
        rows.append(dict(metabolite=f, beta_SD=m.params['y'], lo=ci[0], hi=ci[1], p=m.pvalues['y'],
                         beta_short=m2.params[ks], p_short=m2.pvalues[ks], beta_long=m2.params[kl], p_long=m2.pvalues[kl]))
    except Exception as e:
        print('fail', f, e)
L = pd.DataFrame(rows)
for c in ['p', 'p_short', 'p_long']:
    L['q' + c[1:]] = bh(L[c])
L.to_csv(OUT + 'levels_mixed.csv', index=False)
pd.set_option('display.width', 220)
print(L.sort_values('p').round(3).head(20).to_string())

# (b) symptoms within patients
S = score.copy()
sym = list(S.columns)
ssafe = {s: 's%d' % i for i, s in enumerate(sym)}
P = D[D.y == 1].join(S.rename(columns=ssafe), how='inner')
P[list(ssafe.values())] = (P[list(ssafe.values())] - P[list(ssafe.values())].mean()) / P[list(ssafe.values())].std()
print('patient visits with scores', len(P), 'persons', P.subject.nunique(), 'with >=2 visits', (P.groupby('subject').size() >= 2).sum())
rows = []
for f in feats:
    v = safe[f]
    P['mb'] = P.groupby('subject')[v].transform('mean'); P['mw'] = P[v] - P['mb']
    for s in sym:
        sv = ssafe[s]
        try:
            m = smf.mixedlm(f'{sv} ~ mb + mw + age + male + bmi_i + bmi_na + C(tp)', P, groups=P.subject).fit(reml=True)
            rows.append(dict(metabolite=f, symptom=s, b_between=m.params['mb'], p_between=m.pvalues['mb'],
                             b_within=m.params['mw'], lo_within=m.conf_int().loc['mw', 0], hi_within=m.conf_int().loc['mw', 1],
                             p_within=m.pvalues['mw']))
        except Exception as e:
            print('fail', f, s, e)
Q = pd.DataFrame(rows)
Q['q_between'] = bh(Q.p_between); Q['q_within'] = bh(Q.p_within)
Q.to_csv(OUT + 'symptoms_between_within.csv', index=False)
print('tests', len(Q), 'q_between<.05', (Q.q_between < .05).sum(), 'q_within<.05', (Q.q_within < .05).sum(),
      'q_within<.1', (Q.q_within < .1).sum())
print(Q.sort_values('p_within').head(15).round(4).to_string())
print(Q.sort_values('p_between').head(12).round(4).to_string())
print('\nPrior hint check (indolelactate / indoxyl sulfate vs OI and severity):')
print(Q[Q.metabolite.isin(['indolelactate', '3-indoxyl sulfate', 'indole_composite']) &
        Q.symptom.isin(['NLT.Orthostatic Intolerance', 'severity'])].round(4).to_string())
