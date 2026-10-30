"""Check of the within-person kynurenate ~ CD8 IFN-g/GzA signal: controls too? raw immune percentages? tp adjusted?"""
from cc import *
import statsmodels.formula.api as smf
X, ann = metab(0.5)
out = []
for raw in [False, True]:
    IM = immune(raw=raw)
    for grp in [1, 0]:
        common = IM.index.intersection(X.index).intersection(meta.index)
        M = meta.loc[common]; M = M[M.y == grp]
        M = M[M.groupby('subject').subject.transform('size') >= 2]
        for i in ['% CD8+ RO+ CCR6- total INFg+ d1', '% CD8+ nonMAIT total GzA+ d1', '% CD8+ total INFg+ d1 P8', '% CD8+ total GzA+ d1', '% CD8+ MAIT total INFg+ d1']:
            for m in ['kynurenate', 'kynurenine', 'KT_ratio']:
                P = M.copy()
                v = IM.loc[P.index, i].astype(float)
                if raw: v = np.log(v + 0.1)
                P['iv'] = (v - v.mean()) / v.std()
                y = (X.kynurenine - X.tryptophan) if m == 'KT_ratio' else X[m]
                P['yv'] = ((y - y.mean()) / y.std()).loc[P.index]
                P['ib'] = P.groupby('subject').iv.transform('mean'); P['iw'] = P.iv - P.ib
                f = smf.mixedlm('yv ~ ib + iw + age + male + bmi_i + bmi_na + C(tp)', P, groups=P.subject).fit()
                out.append(dict(raw=raw, group='patients' if grp else 'controls', n_vis=len(P), n_pers=P.subject.nunique(), immune=i, met=m,
                                b_within=f.params['iw'], p_within=f.pvalues['iw'], b_between=f.params['ib'], p_between=f.pvalues['ib']))
O = pd.DataFrame(out); O.to_csv(OUT + 'immune_kyn_within_check.csv', index=False)
pd.set_option('display.width', 220); print(O.round(4).to_string())
