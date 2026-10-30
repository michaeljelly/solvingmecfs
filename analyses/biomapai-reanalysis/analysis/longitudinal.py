from common import *
from scipy import stats
from statsmodels.stats.multitest import multipletests
s=score.copy()
s['RAND36.Role Limitations due to Physical Health']/=4; s['RAND36.Role Limitations due to Emotional Problems']/=3
s['severity']=s.mean(axis=1)
s=s.join(meta[['subject','timepoints','illness_duration','y','age','gender']],how='inner')
s['tp']=s.timepoints.str[2:].astype(int)
p=s[s.y==1].sort_values(['subject','tp'])
# trajectories
first=p.groupby('subject').first(); last=p.groupby('subject').last(); n=p.groupby('subject').size()
d=(last['severity']-first['severity'])[n>=2]
dpf=(last['RAND36.Physical Functioning']-first['RAND36.Physical Functioning'])[n>=2]
dur=first.loc[d.index,'illness_duration']
print('patients with >=2 visits:',len(d))
sd0=p[p.tp==1].severity.std()
print('baseline severity SD in patients:',round(sd0,3))
for grp in ['Short','Long']:
    x=d[dur==grp]
    print(grp,'n',len(x),'mean change',round(x.mean(),3),'improved >0.5SD',(x<-0.5*sd0).sum(),'worsened >0.5SD',(x>0.5*sd0).sum(),'stable',((x.abs()<=0.5*sd0)).sum())
ctrl=s[s.y==0].severity
print('control severity mean',round(ctrl.mean(),3),'patient baseline mean',round(p[p.tp==1].severity.mean(),3))
# how many patients reach control range (below control 95th pct)
thr=ctrl.quantile(0.95)
print('patients whose last visit is within control 95th pct:',(last.loc[d.index,'severity']<=thr).sum())
traj=pd.DataFrame({'n_visits':n[d.index],'duration':dur,'sev_first':first.loc[d.index,'severity'],'sev_last':last.loc[d.index,'severity'],'delta_sev':d,'delta_PF':dpf})
traj.to_csv('patient_trajectories.csv')
# test within-patient: short vs long mean change
print('Short vs Long delta MWU p=',stats.mannwhitneyu(d[dur=='Short'],d[dur=='Long']).pvalue)
# symptom-by-symptom mean change
cols=[c for c in score.columns]
chg=(last[cols]-first[cols])[n>=2]
print((chg.mean()/p[p.tp==1][cols].std()).round(2).sort_values())

# ---- within-person association: demeaned features vs demeaned severity
out=[]
for L in ['immune','specie','kegg','metabolome','quest']:
    X=load(L)
    if L=='specie': X=np.log10(X+1e-5)
    if L=='kegg': X=np.log10(X+1)
    idx=X.index.intersection(p.index)
    pp=p.loc[idx]; keep=pp.groupby('subject').subject.transform('size')>=2
    pp=pp[keep]; Xp=X.loc[pp.index]
    # drop rare/low-var
    Xp=Xp.loc[:,Xp.std()>0]
    if L=='specie': Xp=Xp.loc[:,(X.loc[pp.index]>-4.9).mean()>0.2]
    Xd=Xp-Xp.groupby(pp.subject).transform('mean')
    for target in ['severity','RAND36.Physical Functioning','RAND36.Energy/Fatigue','NLT.Orthostatic Intolerance','PSQI.sleeping_problem','DANA.Cognitive Efficiency','RAND36.Pain','GSRS.gastrointestinal']:
        yd=pp[target]-pp.groupby('subject')[target].transform('mean')
        for f in Xd.columns:
            r,pv=stats.spearmanr(Xd[f],yd)
            out.append((L,target,f,len(pp),pp.subject.nunique(),r,pv))
W=pd.DataFrame(out,columns=['layer','target','feature','n_samples','n_subjects','rho_within','p'])
W=W.dropna()
W['q']=np.nan
for (L,t),g in W.groupby(['layer','target']):
    W.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
W.sort_values('p').to_csv('within_person_associations.csv',index=False)
print(W.groupby(['layer','target']).apply(lambda g:(g.q<0.1).sum()).unstack())
print(W[W.q<0.1].sort_values('p').head(60).to_string())
