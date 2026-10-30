"""Part 4: candidate axes - continuous or grouped? correlated with each other? related to symptoms / case status?"""
from axes import *
from clust_tools import *
import diptest, statsmodels.formula.api as smf
D,info=build()
T=D[D.tp==1].copy(); ref=T.index[T.y==1]
T=axis_scores(T,ref); AX=['axis_lipolysis','axis_kyn_inflam','axis_dysbiosis']
P=T[T.y==1]; Ctl=T[T.y==0]
SYM=['NLT.Orthostatic Intolerance','PSQI.sleeping_problem','severity','RAND36.Energy/Fatigue','GSRS.gastrointestinal']
V=COMP+AX+SYM
def partial_spearman(G,a,b):
    g=G[[a,b,'age','male','bmi_i']].dropna(); R=g.rank()
    Dm=np.column_stack([np.ones(len(g)),R[['age','male','bmi_i']].values])
    ra=R[a]-Dm@np.linalg.lstsq(Dm,R[a],rcond=None)[0]; rb=R[b]-Dm@np.linalg.lstsq(Dm,R[b],rcond=None)[0]
    r,p=stats.pearsonr(ra,rb); return r,p,len(g)
rows=[]
for grp,G in [('patients',P),("controls",Ctl)]:
    for i,a in enumerate(V):
        for b in V[i+1:]:
            if a in SYM and b in SYM: continue
            if grp=='controls' and (a in SYM or b in SYM): continue
            r,p,n=partial_spearman(G,a,b); rows.append((grp,a,b,n,r,p))
R=pd.DataFrame(rows,columns=['group','a','b','n','partial_rho','p']); R['q']=R.groupby('group').p.transform(bh)
R.to_csv(OUT+'results/axes_correlations.csv',index=False)
pd.set_option('display.width',220)
print('--- patients, q<0.1'); print(R[(R.group=='patients')&(R.q<0.1)].sort_values('p').round(4).to_string())
print('--- axis-axis'); print(R[R.a.isin(AX)&R.b.isin(AX)].round(3).to_string())
# unimodality
u=[]
for c in COMP+AX:
    x=P[c].dropna().values; x=(x-x.mean())/x.std()
    d,pd_=diptest.diptest(x)
    b1=GaussianMixture(1).fit(x[:,None]).bic(x[:,None]); b2=GaussianMixture(2,n_init=3,random_state=0).fit(x[:,None]).bic(x[:,None])
    u.append((c,len(x),d,pd_,b1-b2,stats.skew(x)))
U=pd.DataFrame(u,columns=['var','n','dip','dip_p','BIC1_minus_BIC2','skew']); U.to_csv(OUT+'results/axes_unimodality.csv',index=False); print(U.round(3).to_string())
# clustering in 3-axis space
A=P[AX].dropna(); X=((A-A.mean())/A.std()).values; rng=np.random.default_rng(5); cr=[]
for k in (2,3,4):
    s1,_,lab=stability(X,k,'kmeans',B=30); sil=silhouette_score(X,lab)
    ns=[silhouette_score(N,fit(N,k,'kmeans')) for N in [gaussian_null(X,rng) for _ in range(30)]]
    cr.append((k,len(X),s1,sil,np.mean(ns),np.quantile(ns,.95)))
print('3-axis kmeans (k,n,stability,sil,null_sil,null_sil95)',np.round(cr,3).tolist())
# overlap of 'abnormal' patients (beyond control 90th pct in disease direction)
thr={a:Ctl[a].quantile(0.9) for a in AX}; ab=pd.DataFrame({a:P[a]>thr[a] for a in AX}).loc[A.index]
print('fraction of patients above control 90th pct:',ab.mean().round(3).to_dict(),' controls by construction 0.10')
for i,a in enumerate(AX):
    for b in AX[i+1:]:
        t=pd.crosstab(ab[a],ab[b]); o,pf=stats.fisher_exact(t) if t.shape==(2,2) else (np.nan,np.nan)
        print(f'overlap {a} & {b}: both={int((ab[a]&ab[b]).sum())} expected={ab[a].mean()*ab[b].mean()*len(ab):.1f} OR={o:.2f} p={pf:.3f}')
print('patients abnormal on 0/1/2/3 axes:',ab.sum(1).value_counts().sort_index().to_dict(),'n',len(ab))
# case/control and duration per axis, adjusted
cc=[]
for a in AX+COMP:
    g=T[[a,'illness_duration','age','male','bmi_i']].dropna().rename(columns={a:'v'}); g['v']=(g.v-g.v[g.illness_duration=='Control'].mean())/g.v[g.illness_duration=='Control'].std()
    m=smf.ols("v ~ C(illness_duration,Treatment('Control'))+age+male+bmi_i",g).fit()
    for lv in ['Short','Long']:
        k=f"C(illness_duration, Treatment('Control'))[T.{lv}]"; ci=m.conf_int().loc[k]
        cc.append((a,lv,m.params[k],ci[0],ci[1],m.pvalues[k]))
CC=pd.DataFrame(cc,columns=['var','group_vs_control','diff_SD','lo','hi','p']); CC['q']=bh(CC.p); CC.to_csv(OUT+'results/axes_case_control.csv',index=False)
print(CC.round(3).to_string())
# axes vs symptoms (adjusted OLS, z units)
sy=[]
for a in AX:
    for s in SYM:
        g=P[[a,s,'age','male','bmi_i']].dropna(); g=(g-g.mean())/g.std(); g.columns=['x','y','age','male','bmi']
        m=smf.ols('y~x+age+male+bmi',g).fit(); ci=m.conf_int().loc['x']; sy.append((a,s,len(g),m.params.x,ci[0],ci[1],m.pvalues.x))
SY=pd.DataFrame(sy,columns=['axis','symptom','n','std_beta','lo','hi','p']); SY['q']=bh(SY.p); SY.to_csv(OUT+'results/axes_vs_symptoms.csv',index=False); print(SY.round(3).to_string())
T.to_csv(OUT+'results/axes_tp1_values.csv')
print('R_gnavus detected (patients, controls):',(P.R_gnavus>-5+1e-9).mean().round(3),(Ctl.R_gnavus>-5+1e-9).mean().round(3))
