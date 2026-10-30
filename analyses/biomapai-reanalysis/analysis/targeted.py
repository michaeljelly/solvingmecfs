from common import *
from scipy import stats
import statsmodels.formula.api as smf
mm=pd.read_csv('/mnt/project-files/biomap/BioMapAI/codes/AI/input/feature_meta/metabolome_feature_meta.csv',index_col=0)
X=load('metabolome')
s=score.copy(); s['RAND36.Role Limitations due to Physical Health']/=4; s['RAND36.Role Limitations due to Emotional Problems']/=3
s['severity']=s.mean(axis=1)
M=meta.join(s,how='inner').join(X,how='inner',rsuffix='_m')
M['tp']=M.timepoints.str[2:].astype(int)
def z(v): return (v-v.mean())/v.std()
sp=mm.SUB_PATHWAY.reindex(X.columns)
ffa=[c for c in X.columns if sp.get(c) in ('Long Chain Polyunsaturated Fatty Acid (n3 and n6)','Long Chain Saturated Fatty Acid','Long Chain Monounsaturated Fatty Acid')]
hfa=[c for c in X.columns if sp.get(c)=='Fatty Acid, Monohydroxy']
acyl=[c for c in X.columns if str(sp.get(c)).startswith('Fatty Acid Metabolism (Acyl Carnitine')]
print('n FFA',len(ffa),'hydroxyFA',len(hfa),'acylcarn',len(acyl))
print([c for c in X.columns if 'tryptophan' in c.lower() or 'kynuren' in c.lower() or 'indol' in c.lower() or 'serotonin' in c.lower()])
M['FFA']=X.loc[M.index,ffa].apply(z).mean(axis=1)
M['hydroxyFA']=X.loc[M.index,hfa].apply(z).mean(axis=1)
M['acylcarnitine']=X.loc[M.index,acyl].apply(z).mean(axis=1)
M['ketone_BHB']=X.loc[M.index,'3-hydroxybutyrate (BHBA)']
M['KT']=X.loc[M.index,'kynurenine']-X.loc[M.index,'tryptophan']  # values already scaled/logged? difference of scaled values
ind=[c for c in ['indolelactate','3-indoxyl sulfate','6-hydroxyindole sulfate','indolepropionate','indoleacetate','indole-3-carboxylic acid'] if c in X.columns]
print('indoles used',ind)
M['microbial_indoles']=X.loc[M.index,ind].apply(z).mean(axis=1)
comp=['FFA','hydroxyFA','acylcarnitine','ketone_BHB','KT','kynurenine','tryptophan','microbial_indoles']
targets=['NLT.Orthostatic Intolerance','severity','RAND36.Energy/Fatigue','PSQI.sleeping_problem','DANA.Cognitive Efficiency']
out=[]
for grp,G in [('patients_tp1',M[(M.y==1)&(M.tp==1)]),('controls_tp1',M[(M.y==0)&(M.tp==1)])]:
    for c in comp:
        for t in targets:
            r,p=stats.spearmanr(G[c],G[t]); out.append((grp,c,t,len(G),round(r,3),p))
P=M[M.y==1].copy(); P=P[P.groupby('subject').subject.transform('size')>=2]
for c in comp:
    cd=P[c]-P.groupby('subject')[c].transform('mean')
    for t in targets:
        td=P[t]-P.groupby('subject')[t].transform('mean')
        r,p=stats.spearmanr(cd,td); out.append(('patients_within_person',c,t,len(P),round(r,3),p))
O=pd.DataFrame(out,columns=['group','composite','target','n','rho','p'])
O.to_csv('targeted_composites.csv',index=False)
pd.set_option('display.width',200)
print(O.pivot_table(index=['composite','target'],columns='group',values='rho').round(2).to_string())
print(O[O.p<0.01].to_string())
# case-control levels at tp1
T=M[M.tp==1]
for c in comp:
    a=T[T.illness_duration=='Short'][c];l=T[T.illness_duration=='Long'][c];h=T[T.illness_duration=='Control'][c]
    print(f"{c:18s} ctrl {h.mean():+.2f} short {a.mean():+.2f} (p={stats.mannwhitneyu(a,h).pvalue:.3g}) long {l.mean():+.2f} (p={stats.mannwhitneyu(l,h).pvalue:.3g})")
# mixed model: OI ~ FFA adjusting age, gender, bmi in patients all tps
mm_=smf.mixedlm('Q("NLT.Orthostatic Intolerance") ~ FFA + age + C(gender)',P.dropna(subset=['FFA']),groups=P.dropna(subset=['FFA']).subject).fit()
print(mm_.summary().tables[1])
mm_=smf.mixedlm('severity ~ KT + microbial_indoles + age + C(gender)',P,groups=P.subject).fit()
print(mm_.summary().tables[1])
