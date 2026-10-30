from common import *
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
rows=[]
for L in LAYERS:
    X=load(L); X=X.loc[X.index.intersection(meta.index)]
    y=meta.loc[X.index,'y'].values; g=meta.loc[X.index,'subject'].values
    if L in('specie','kegg'): X=np.log10(X+1e-6) if L=='specie' else np.log10(X+1)
    X=X.fillna(0).values
    for name,cv in [('sample-split (as paper)',StratifiedKFold(5,shuffle=True,random_state=0)),('person-split',StratifiedGroupKFold(5,shuffle=True,random_state=0))]:
        for mname,mk in [('RF',lambda:RandomForestClassifier(500,min_samples_leaf=3,n_jobs=-1,random_state=0)),('LR',lambda:make_pipeline(StandardScaler(),LogisticRegression(C=0.05,penalty='l2',max_iter=5000)))]:
            aucs=[]
            for rep in range(3):
                cv.random_state=rep
                p=np.zeros(len(y))
                for tr,te in cv.split(X,y,g):
                    m=mk().fit(X[tr],y[tr]); p[te]=m.predict_proba(X[te])[:,1]
                aucs.append(roc_auc_score(y,p))
            rows.append((L,X.shape[0],X.shape[1],name,mname,np.mean(aucs)))
            print(rows[-1],flush=True)
r=pd.DataFrame(rows,columns=['layer','n_samples','n_features','split','model','AUC'])
r.to_csv('leakage_results.csv',index=False)
print(r.pivot_table(index=['layer','model'],columns='split',values='AUC').round(3))
