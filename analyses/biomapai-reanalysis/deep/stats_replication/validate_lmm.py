import numpy as np, pandas as pd, statsmodels.formula.api as smf, time
from data import *
from lmmfast import fit_many
m=load_meta()
for L in ['metabolome','species']:
    X=load_layer(L); idx=X.index.intersection(m.index); M=m.loc[idx]; X=X.loc[idx]
    X=(X-X.mean())/X.std()
    cols=list(X.columns[:6])
    C=np.c_[np.ones(len(M)),M.case,M.age_c,M.male,M.bmi_c,M.bmi_missing,M.tp2,M.tp3]
    t=time.time(); r=fit_many(X.values,C,M.person.values,L=np.eye(8)[[1]]); print(L,X.shape,'fast time',round(time.time()-t,2))
    for j,c in enumerate(cols):
        d=M.assign(yv=X[c].values)
        f=smf.mixedlm('yv~case+age_c+male+bmi_c+bmi_missing+tp2+tp3',d,groups=d.person).fit(reml=True)
        print(c[:30].ljust(30), 'sm b=%.4f se=%.4f p=%.4g | fast b=%.4f se=%.4f p=%.4g icc=%.2f'%(f.params['case'],f.bse['case'],f.pvalues['case'],r['est'][j,0],r['se'][j,0],r['p'][j,0],r['icc'][j]))
    # symptom model check
    P=M[(M.case==1)&M.severity.notna()]; XP=X.loc[P.index]
    for c in cols[:3]:
        bt=XP[c].groupby(P.person).transform('mean'); wi=XP[c]-bt
        d=P.assign(b=bt.values,w=wi.values)
        f=smf.mixedlm('severity~b+w+age_c+male',d,groups=d.person).fit(reml=True)
        Fm=np.stack([bt.values,wi.values],1)[None]
        r2=fit_many(P.severity.values,np.c_[np.ones(len(P)),P.age_c,P.male],P.person.values,F=Fm)
        print(' sym',c[:20],'sm',f.params[['b','w']].round(4).tolist(),f.bse[['b','w']].round(4).tolist(),'fast',r2['est'][0].round(4).tolist(),r2['se'][0].round(4).tolist())
