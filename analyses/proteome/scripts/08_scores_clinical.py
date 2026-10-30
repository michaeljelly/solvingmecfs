"""Per-person tissue scores: confounding by activity and BMI, links to severity, illness duration and exertion response."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
import statsmodels.formula.api as smf
rd=lambda f: set(pd.read_csv(W+'raw/hpa/'+f,sep='\t').iloc[:,0].astype(str))
loc=pd.read_csv(D+'hoel_gene_location.csv',index_col=0).iloc[:,0]
liver=rd('T_liver.tsv')|rd('C_hepatocytes.tsv'); muscle=rd('T_skeletal_muscle.tsv'); plate=rd('C_platelets.tsv'); entero=rd('C_enterocytes.tsv')|rd('T_intestine.tsv')
SETS={'liver_intra':{g for g in liver if loc.get(g)=='Intracellular'},'liver_secr':{g for g in liver if loc.get(g)=='Secreted'},
      'muscle_intra':{g for g in muscle if loc.get(g)=='Intracellular'},'platelet_intra':{g for g in plate if loc.get(g)=='Intracellular'},
      'enterocyte':entero}
def scores(G):
    return pd.DataFrame({k:G[[g for g in s if g in G.columns]].mean(1) for k,s in SETS.items()})
out=[]
def fit(df,formula,label,term):
    try:
        r=smf.ols(formula,df).fit(); out.append(dict(set=formula.split('~')[0].strip(),analysis=label,term=term,beta=r.params.get(term,np.nan),p=r.pvalues.get(term,np.nan),n=int(r.nobs)))
    except Exception as e: print(label,e)
# ---- Hoel
X,m,f=load('hoel_serum'); m=m.join(pd.read_csv(D+'hoel_serum_leak.csv',index_col=0).leak)
G=gene_collapse(X,f); S=scores(G); S=(S-S.mean())/S.std(); df=m.join(S)
df['steps_k']=df.Daily_Step_Count/1000
for s in SETS:
    fit(df,f'{s}~case+male+age+bmi+leak','hoel: case adj sex/age/BMI/leak','case')
    fit(df[df.case==1],f'{s}~steps_k+male+age+bmi+leak','hoel: steps within ME','steps_k')
    fit(df[df.case==1],f'{s}~SF36PF+male+age+bmi+leak','hoel: SF36PF within ME','SF36PF')
    fit(df,f'{s}~bmi+case+male+age+leak','hoel: BMI','bmi')
print('hoel steps by group', df.groupby('case').steps_k.describe().round(1).to_string())
# ---- Germain
X,m,f=load('germain_plasma'); G=gene_collapse(X,f); S=scores(G)
df=m.join(S); df['bas']=pd.to_numeric(df.bas_score,errors='coerce'); df['recl']=pd.to_numeric(df.q_reclined,errors='coerce'); df['PCS']=pd.to_numeric(df.PCS,errors='coerce')
df['long']=(df.mecfs_duration.astype(str).str.contains('>=')).astype(float); df.loc[df.case==0,'long']=np.nan
base=df[df.timepoint=='D1PRE'].copy()
for s in SETS: base[s]=(base[s]-base[s].mean())/base[s].std()
for s in SETS:
    fit(base,f'{s}~case+male+age+C(bmi_bin)+C(site)','germain: case','case')
    fit(base,f'{s}~case+male+age+C(site)','germain: case, no BMI','case')
    fit(base[base.case==1],f'{s}~bas+male+age+C(bmi_bin)+C(site)','germain: Bell activity within ME','bas')
    fit(base[base.case==1],f'{s}~recl+male+age+C(bmi_bin)+C(site)','germain: % reclined within ME','recl')
    fit(base[base.case==1],f'{s}~PCS+male+age+C(bmi_bin)+C(site)','germain: SF36 PCS within ME','PCS')
    fit(base[base.case==1],f'{s}~long+male+age+C(bmi_bin)+C(site)','germain: illness >=5y vs <5y','long')
    fit(base,f'{s}~bmi_bin+case+male+age+C(site)','germain: BMI bin (linear)','bmi_bin')
# exertion: change D1POST-D1PRE and D2PRE-D1PRE, case vs control
W2=df.pivot_table(index='cor_id',columns='timepoint',values=list(SETS))
dm=df.drop_duplicates('cor_id').set_index('cor_id')
for s in SETS:
    for a,b,lab in [('D1POST','D1PRE','acute (D1 post - pre)'),('D2PRE','D1PRE','24h recovery (D2 pre - D1 pre)'),('D2POST','D2PRE','2nd test acute'),('D2POST','D1POST','post D2 vs post D1')]:
        ch=(W2[(s,a)]-W2[(s,b)]).rename('ch'); d2=dm.join(ch).dropna(subset=['ch']); d2['ch']=d2.ch/d2.ch.std()
        d2=d2.rename(columns={'ch':s+'_ch'}) if False else d2
        for lbl,sub,term,fm in [('case-vs-control difference in change',d2,'case','ch~case+male+age+C(site)'),('mean change HC',d2[d2.case==0],'Intercept','ch~1'),('mean change ME',d2[d2.case==1],'Intercept','ch~1')]:
            r=smf.ols(fm,sub).fit(); out.append(dict(set=s,analysis=f'germain exertion {lab}: {lbl}',term=term,beta=r.params[term],p=r.pvalues[term],n=int(r.nobs)))
O=pd.DataFrame(out)
# attach set names in order
O.to_csv(R+'tissue_scores_clinical.csv',index=False)

print(O.pivot_table(index='analysis',columns='set',values='beta').round(2).to_string())
print(O.pivot_table(index='analysis',columns='set',values='p').round(3).to_string())
