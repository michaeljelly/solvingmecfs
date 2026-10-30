from common import *
from scipy import stats
X=np.log10(load('specie')+1e-5)
bp=[c for c in X.columns if any(k in c for k in ['Faecalibacterium','Roseburia','Coprococcus','Agathobacter','Eubacterium_rectale','Anaerostipes','Butyricicoccus'])]
print(len(bp),bp)
s=score.copy(); s['RAND36.Role Limitations due to Physical Health']/=4; s['RAND36.Role Limitations due to Emotional Problems']/=3; s['severity']=s.mean(axis=1)
M=meta.join(s,how='inner'); M=M.loc[M.index.intersection(X.index)]; M['tp']=M.timepoints.str[2:].astype(int)
raw=load('specie').loc[M.index]
M['butyrate_producers_pct']=raw[bp].sum(axis=1)*100 if raw[bp].sum(axis=1).max()<=1.01 else raw[bp].sum(axis=1)
M['R_gnavus']=X.loc[M.index,'Ruminococcus_gnavus']
print(raw.sum(axis=1).describe())
T=M[M.tp==1]
print(T.groupby('illness_duration').butyrate_producers_pct.median())
for g in ['Short','Long']: print(g,stats.mannwhitneyu(T[T.illness_duration==g].butyrate_producers_pct,T[T.illness_duration=='Control'].butyrate_producers_pct).pvalue)
P=T[T.y==1]
for c in ['butyrate_producers_pct','R_gnavus']:
    for t in ['severity','GSRS.gastrointestinal','RAND36.Energy/Fatigue','PSQI.sleeping_problem','NLT.Orthostatic Intolerance','RAND36.Emotional wellbeing']:
        r,p=stats.spearmanr(P[c],P[t]); 
        PP=M[(M.y==1)]; PP=PP[PP.groupby('subject').subject.transform('size')>=2]
        r2,p2=stats.spearmanr(PP[c]-PP.groupby('subject')[c].transform('mean'),PP[t]-PP.groupby('subject')[t].transform('mean'))
        print(f"{c:24s} {t:30s} between r={r:+.2f} p={p:.3f} | within r={r2:+.2f} p={p2:.3f}")
# probiotics / antibiotics use and severity
print(M[M.y==1].groupby('probiotics').severity.agg(['mean','size']))
print(M[M.y==1].groupby('antivirals').severity.agg(['mean','size']))
