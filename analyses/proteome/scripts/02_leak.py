"""Cell-leak / pre-analytical artefact: define a score in Hoel, test whether the same score differs in other SomaScan cohorts."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
import openpyxl
X,m,f=load('hoel_serum')
# location annotation from Hoel SuppData3
r=list(openpyxl.load_workbook(W+'raw/supp/PMC13006441/mmc4.xlsx',read_only=True)['2_ME vs. HC'].iter_rows(values_only=True))
loc=pd.DataFrame(r[1:],columns=r[0]).set_index('AptName')
print(loc.Loc.value_counts())
seed=['ALDOA','LDHB','VIM','H2AZ1','H1-2','H1-10','H2BC12','ATP5IF1']
Z=(X-X.mean())/X.std()
seedapt=[a for a in X.columns if str(f.loc[a,'gene']) in seed]
s0=Z[seedapt].mean(1)
# correlate all intracellular aptamers with seed score within controls+cases (adjusted for group)
intra=[a for a in X.columns if loc.Loc.get(a)=='Intracellular']
res=Z[intra].apply(lambda c: np.corrcoef(c, s0)[0,1])
top=res.sort_values(ascending=False).index[:100].tolist()
print('top leak aptamers',[f.loc[a,'gene'] for a in top[:30]])
pd.Series(top).to_csv(R+'leak_aptamers_top100.csv',index=False)
leak=Z[top].mean(1)
print('Hoel leak score d (case-control, raw diff/SD):', round((leak[m.case==1].mean()-leak[m.case==0].mean())/leak.std(),2))
m['leak']=leak; m[['leak']].to_csv(D+'hoel_serum_leak.csv')
# PCA
Zf=Z.fillna(0).values; U,S,Vt=np.linalg.svd(Zf-Zf.mean(0),full_matrices=False); pc1=U[:,0]*S[0]
print('Hoel PC1 var', round(S[0]**2/(S**2).sum(),3),'r(PC1,leak)',round(np.corrcoef(pc1,leak)[0,1],2),'r(PC1,case)',round(np.corrcoef(pc1,m.case)[0,1],2))
out=[]
for ds in ['germain_plasma','walitt_serum','walitt_csf']:
    Xg,mg,fg=load(ds)
    if ds=='germain_plasma':
        sel=mg.timepoint=='D1PRE'; Xg=Xg[sel]; mg=mg[sel]
    common=[a for a in top if a in Xg.columns]
    if len(common)<10:
        # map by gene symbol (1.3k uses different seq suffixes)
        gtop=set(str(f.loc[a,'gene']) for a in top)
        common=[a for a in Xg.columns if str(fg.loc[a,'gene']) in gtop]
    Zg=(Xg-Xg.mean())/Xg.std(); L=Zg[common].mean(1)
    d=(L[mg.case==1].mean()-L[mg.case==0].mean())/L.std(); p=stats.mannwhitneyu(L[mg.case==1],L[mg.case==0]).pvalue
    print(ds,'n aptamers',len(common),'leak d',round(d,2),'p',round(p,4))
    out.append(dict(dataset=ds,n_aptamers=len(common),leak_d=d,p=p))
    mg['leak']=L; mg[['leak']].to_csv(D+f'{ds}_leak.csv')
pd.DataFrame(out).to_csv(R+'leak_score_other_cohorts.csv',index=False)
