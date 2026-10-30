import pandas as pd, numpy as np, warnings
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
fm=pd.read_pickle(OUT+'fm.pkl'); gm=pd.read_pickle(OUT+'gm.pkl'); X=pd.read_pickle(OUT+'X_scaled.pkl'); meta=pd.read_pickle(OUT+'meta.pkl')
X=X.loc[:,X.columns!='nan']
Lall=np.log2(X)
print('== Paper-style: Wilcoxon per timepoint and sex, all 1157 features, BH')
rows=[]
for sex in ['F','M']:
    for t in ['D1-PRE','D1-POST','D2-PRE','D2-POST']:
        s=(gm.Timepoint==t)&(gm.Sex==sex); A=Lall[s&(gm.Phenotype=='CFS')]; B=Lall[s&(gm.Phenotype=='Control')]
        p=np.array([stats.mannwhitneyu(A[f],B[f]).pvalue if Lall.loc[s,f].nunique()>1 else 1 for f in Lall.columns])
        q=multipletests(p,method='fdr_bh')[1]
        rows.append((sex,t,len(A),len(B),(q<0.05).sum(),(q<0.15).sum(),(p<0.05).sum()))
        if sex=='F' and t=='D1-PRE':
            top=pd.DataFrame({'p':p,'q':q,'l2fc':A.median().values-B.median().values},index=Lall.columns).sort_values('p').head(12); print(top.round(4).to_string())
print(pd.DataFrame(rows,columns=['sex','tp','nME','nCtrl','q<.05','q<.15','p<.05']).to_string(index=False))
# paper's key named metabolites
C=pd.read_pickle(OUT+'contrasts.pkl')
Lf=Lall
def wide(t): s=gm.Timepoint==t; return Lf[s].set_index(gm.ID[s].values)
W={t:wide(t) for t in ['D1-PRE','D1-POST','D2-PRE','D2-POST']}
print('\n== Named metabolites from paper: within-person change, females, Wilcoxon ME vs ctrl, and mean log2 change per group')
for f in ['alpha-ketoglutarate','lactate','pyruvate','1,5-anhydroglucitol (1,5-AG)','tridecenedioate (C13:1-DC)*','indoleacetoylcarnitine','glutamate','glutamine','X - 15245','X - 16397','X - 07765']:
    if f not in Lf.columns: print(f,'absent'); continue
    out=[f]
    for nm,(a,b) in {'R1':('D1-POST','D1-PRE'),'REC':('D2-PRE','D1-POST'),'R2':('D2-POST','D2-PRE')}.items():
        d=(W[a][f]-W[b][f]).dropna(); m=meta.loc[d.index]
        for sex in ['F','M']:
            me=d[(m.case==1)&(m.Sex==sex)]; ct=d[(m.case==0)&(m.Sex==sex)]
            out.append(f'{nm}-{sex}: ctrl {ct.mean():+.2f} ME {me.mean():+.2f} p={stats.mannwhitneyu(me,ct).pvalue:.3f}')
    print(' | '.join(out))
# global hemoconcentration: median log2 change across all robust metabolites per person
L=pd.read_pickle(OUT+'L.pkl')
print('\n== Global shift (median log2 change across 782 well-detected metabolites) per person, group means')
for k in ['R1 (D1 post-pre)','R2 (D2 post-pre)','REC (D2pre-D1pre)','D1post->D2pre']:
    g=C[k].median(axis=1); m=meta.loc[g.index]
    print(f'{k:20s} ctrl {g[m.case==0].mean():+.3f}  ME {g[m.case==1].mean():+.3f}  p={stats.mannwhitneyu(g[m.case==1],g[m.case==0]).pvalue:.3f}')
# lactate/pyruvate ratio
for nm in ['LEVEL D1-PRE','LEVEL D1-POST','LEVEL D2-PRE','LEVEL D2-POST','R1 (D1 post-pre)','R2 (D2 post-pre)']:
    Y=C[nm]; r=(Y['lactate']-Y['pyruvate']).dropna(); m=meta.loc[r.index]
    print(f'lactate/pyruvate log2 {nm:18s} ctrl {r[m.case==0].mean():+.3f} ME {r[m.case==1].mean():+.3f} p={stats.mannwhitneyu(r[m.case==1],r[m.case==0]).pvalue:.3f}')
