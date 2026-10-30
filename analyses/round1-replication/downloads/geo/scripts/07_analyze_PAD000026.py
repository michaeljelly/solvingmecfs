import pandas as pd, numpy as np, sys
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
R='/tmp/work_geo/raw/'; O='/tmp/work_geo/results/'
L=open(R+'PAD000026.adat').read().split('\n')
# annotation rows 42.. : first nonempty col holds row name
ann={}
for i in range(41,60):
    f=L[i].split('\t'); ann[f[33]]=f
seq=ann['SeqId']; start=seq.index('SeqId')+1
sym=ann.get('EntrezGeneSymbol'); tgt=ann.get('Target') or ann.get('TargetFullName')
A=pd.DataFrame({'seq':seq[start:],'sym':sym[start:] if sym else None})
print(list(ann.keys()))
A['col']='seq.'+A.seq.str.replace('-','.'); A=A.set_index('col')
d=pd.read_csv(R+'PAD000026_Log2Matrix.tsv',sep='\t',decimal=',')
meta=d.iloc[:,:9].copy(); X=d.iloc[:,9:].apply(pd.to_numeric,errors='coerce')
meta['case']=(meta.Group=='ME').astype(int); meta['male']=(meta.Sex.astype(str)=='0').astype(int) if meta.Sex.nunique()>1 else 0
age_map={a:i for i,a in enumerate(sorted(meta.Age.unique()))}; meta['age']=meta.Age.map(age_map)
print(meta.Sex.value_counts().to_dict(), meta.Group.value_counts().to_dict())
meta.index=[f's{i}' for i in range(len(meta))]; X.index=meta.index
# collapse aptamers to gene symbol: keep all aptamers but label with symbol
X=X.T; X['sym']=A.reindex(X.index)['sym'].values
lm=X.drop(columns='sym'); lm.index=[f"{s}|{i}" for s,i in zip(X.sym,X.index)]
# gene-level: mean of aptamers per symbol (after z within aptamer)
z=zscore_rows(lm); z['sym']=[i.split('|')[0] for i in z.index]; gl=z.groupby('sym').mean()
g,s,sc,miss=run_case_control('PAD000026_SomaScan',gl,meta,covars=('male','age'))
g.to_csv(O+'PAD000026_proteins_case_control.csv'); s.to_csv(O+'PAD000026_sets_case_control.csv')
print(s.round(3)[['beta','d','p','q','n_genes']].to_string()); print(g[g.p<0.1].sort_values('p').round(3)[['beta','d','p','q']].to_string())
print('missing',miss)
# SF36 association within cases
cm=meta[meta.case==1].copy(); rows=[]
for f in list(sc.columns):
    y=sc.loc[cm.index,f]; ok=cm.SF36PF.notna(); r=stats.spearmanr(y[ok],cm.SF36PF[ok]); rows.append(dict(feature=f,rho_SF36PF=r.correlation,p=r.pvalue,n=ok.sum()))
r=pd.DataFrame(rows); r['q']=bh(r.p); print(r.round(3).to_string()); r.to_csv(O+'PAD000026_sets_vs_SF36PF_cases.csv',index=False)
# genome-wide aptamer-level test for context
gw=test_matrix(zscore_rows(lm),meta,('male','age')); gw['sym']=[i.split('|')[0] for i in gw.index]
gw.to_csv(O+'PAD000026_all_aptamers_case_control.csv')
print('aptamers',len(gw),'q<0.05:',(gw.q<0.05).sum(),' q<0.1:',(gw.q<0.1).sum())
print(gw.sort_values('p').head(30).round(4)[['sym','d','p','q']].to_string())
tg=gw[gw.sym.isin(TARGET_GENES)].sort_values('p'); print(tg.head(15).round(4)[['sym','d','p','q']].to_string())
cp=competitive(gl,meta,('male','age')); cp.to_csv(O+'PAD000026_sets_competitive.csv'); print(cp.round(3).to_string())
# is intracellular-leak axis the driver? correlate FAO score with a generic intracellular-protein score (mean of top-abundance cytosolic markers)
# sensitivity: adjust for leading principal components (pre-analytical / cell-leak axis)
Z=zscore_rows(lm).fillna(0); U,Sv,Vt=np.linalg.svd(Z.values-Z.values.mean(1,keepdims=True),full_matrices=False)
pcs=pd.DataFrame(Vt[:3].T,index=Z.columns,columns=['PC1','PC2','PC3']); meta2=meta.join(pcs)
print('var expl',(Sv[:5]**2/np.sum(Sv**2)).round(3)); print('PC vs case r:',[round(np.corrcoef(meta2[p],meta2.case)[0,1],2) for p in ['PC1','PC2','PC3']])
leak=[i for i in lm.index if i.split('|')[0] in ('ALDOA','LDHB','VIM','H2AZ1','H1-2','H1-10','H2BC12','ATP5IF1')]
meta2['leak']=zscore_rows(lm.loc[leak]).mean(0); print('leak score d case:',round(ols_test(meta2.leak,meta2,'case',('male','age'))['d'],2))
g2,s2,_,_=run_case_control('PAD000026_SomaScan_leakadj',gl,meta2,covars=('male','age','leak'),extra_label='_leak_adjusted')
g2.to_csv(O+'PAD000026_proteins_case_control_leakadj.csv'); s2.to_csv(O+'PAD000026_sets_case_control_leakadj.csv')
print(s2.round(3)[['beta','d','p','q']].to_string()); print(g2[g2.p<0.05].sort_values('p').round(3)[['beta','d','p','q']].to_string())
