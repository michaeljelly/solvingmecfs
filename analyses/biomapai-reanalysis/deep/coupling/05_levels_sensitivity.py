"""Sensitivity analyses for the case-control differences in microbe-derived metabolites (task 3):
diet adjustment, excluding antibiotic/probiotic visits, tp1-only OLS, authors' processed metabolome table,
and gut 'mediation' (adding Shannon diversity + gut-predicted metabolite score from species) in visits with stool."""
from cc import *
import statsmodels.formula.api as smf
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold

X, ann = metab(0.5)
KEY = ['indolepropionate', 'indolelactate', 'hippurate', 'catechol sulfate', 'cinnamoylglycine', '3-phenylpropionate (hydrocinnamate)',
       'p-cresol glucuronide*', 'glycodeoxycholate 3-sulfate', 'taurodeoxycholate', 'deoxycholate', 'kynurenate']
Z = X[KEY]; Z = (Z - Z.mean()) / Z.std()
safe = {f: 'm%d' % i for i, f in enumerate(KEY)}
D = meta.join(Z.rename(columns=safe), how='inner')
for c in ['diet_meat', 'diet_sugar', 'diet_veg', 'diet_grains', 'diet_fruit']:
    D[c] = D[c].fillna(D[c].median())
full, _ = species_full()
p = full.where(full > 0)
shannon = -(p * np.log(p)).sum(axis=1)
D['shannon'] = shannon.reindex(D.index)
spc, _ = species_clr()
auth = load_auth = pd.read_csv(I + 'metabolome.csv', index_col=0).T
base = 'age + male + bmi_i + bmi_na + C(tp)'
diet = ' + diet_meat + diet_sugar + diet_veg + diet_grains + diet_fruit'
rows = []


def fit(f, data, form, term='y'):
    m = smf.mixedlm(form, data, groups=data.subject).fit(reml=True)
    return m.params[term], m.pvalues[term]


for f in KEY:
    v = safe[f]; r = dict(metabolite=f)
    r['base_b'], r['base_p'] = fit(f, D, f'{v} ~ y + {base}')
    r['diet_b'], r['diet_p'] = fit(f, D, f'{v} ~ y + {base}{diet}')
    Dn = D[(D.antibiotics == 'No') & (D.probiotics == 'No')]
    r['noabx_prob_b'], r['noabx_prob_p'] = fit(f, Dn, f'{v} ~ y + {base}{diet}')
    T1 = D[D.tp == 1]
    o = smf.ols(f'{v} ~ y + age + male + bmi_i + bmi_na{diet}', T1).fit()
    r['tp1_b'], r['tp1_p'] = o.params['y'], o.pvalues['y']
    if f in auth.columns:
        Da = meta.join(((auth[[f]] - auth[f].mean()) / auth[f].std()).rename(columns={f: 'a'}), how='inner')
        r['authors_table_b'], r['authors_table_p'] = fit(f, Da, f'a ~ y + {base}')
    # gut mediation within visits that have stool
    Ds = D.dropna(subset=['shannon'])
    Ds = Ds[Ds.index.isin(spc.index)]
    r['stool_n'] = len(Ds)
    r['stool_b'], r['stool_p'] = fit(f, Ds, f'{v} ~ y + {base}{diet}')
    r['stool_shannon_b'], r['stool_shannon_p'] = fit(f, Ds, f'{v} ~ y + shannon + {base}{diet}')
    # out-of-fold gut-predicted metabolite (ridge on species CLR, trained across both groups, person-grouped CV)
    Xs = spc.loc[Ds.index].values; yv = Ds[v].values; oof = np.zeros(len(yv))
    for tr, te in GroupKFold(5).split(Xs, yv, Ds.subject):
        mu, sd = Xs[tr].mean(0), Xs[tr].std(0) + 1e-9
        mdl = RidgeCV(alphas=np.logspace(-1, 5, 13)).fit((Xs[tr] - mu) / sd, yv[tr])
        oof[te] = mdl.predict((Xs[te] - mu) / sd)
    Ds = Ds.assign(gutpred=oof)
    r['gutpred_r'] = np.corrcoef(oof, yv)[0, 1]
    r['stool_gutpred_b'], r['stool_gutpred_p'] = fit(f, Ds, f'{v} ~ y + gutpred + shannon + {base}{diet}')
    r['attenuation_pct'] = 100 * (1 - r['stool_gutpred_b'] / r['stool_b'])
    rows.append(r)
O = pd.DataFrame(rows)
O.to_csv(OUT + 'levels_sensitivity.csv', index=False)
pd.set_option('display.width', 260)
print(O.round(3).to_string())
