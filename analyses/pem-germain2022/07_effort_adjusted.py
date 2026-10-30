# Patients reach lower workloads; lactate rise is a within-sample proxy for exercise intensity reached.
# Re-test case vs control response differences with each person's own lactate rise as covariate.
import pandas as pd, numpy as np, warnings
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
OUT='/mnt/project-files/mecfs-research/pem/'
C=pd.read_pickle(OUT+'contrasts.pkl'); meta=pd.read_pickle(OUT+'meta.pkl'); fm=pd.read_pickle(OUT+'fm.pkl')
sub=pd.Series(fm['SUB PATHWAY'].values,index=fm.BIOCHEMICAL.astype(str))
out=[]
for k in ['R1 (D1 post-pre)','R2 (D2 post-pre)']:
    Y=C[k].dropna(); m=meta.loc[Y.index].copy(); m['lac']=Y['lactate']
    rows=[]
    for f in Y.columns:
        if f=='lactate' or Y[f].std()<1e-6: continue
        r0=sm.OLS(Y[f],sm.add_constant(m[['case','male','AgeBin','BMIBin']].astype(float))).fit(cov_type='HC3')
        r1=sm.OLS(Y[f],sm.add_constant(m[['case','male','AgeBin','BMIBin','lac']].astype(float))).fit(cov_type='HC3')
        rows.append((f,r0.params['case'],r0.pvalues['case'],r1.params['case'],r1.pvalues['case'],r1.params['lac'],r1.pvalues['lac']))
    R=pd.DataFrame(rows,columns=['f','case_unadj','p_unadj','case_adj','p_adj','lac_beta','p_lac']).set_index('f')
    for c in ['p_unadj','p_adj','p_lac']: R[c.replace('p_','q_')]=multipletests(R[c],method='fdr_bh')[1]
    R['sub']=sub.reindex(R.index).values; R['contrast']=k; out.append(R)
    print(f'\n{k}: metabolites whose exercise change tracks lactate rise (q<0.05): {(R.q_lac<0.05).sum()} of {len(R)}; case effect q<0.1 unadj {(R.q_unadj<0.1).sum()} adj {(R.q_adj<0.1).sum()}; p<0.05 unadj {(R.p_unadj<0.05).sum()} adj {(R.p_adj<0.05).sum()}')
    print(R.sort_values('p_adj').head(10)[['case_unadj','p_unadj','case_adj','p_adj','q_adj','lac_beta','sub']].round(4).to_string())
    print(R.loc[[x for x in ['glucose','pyruvate','cysteinylglycine disulfide*','cys-gly, oxidized','citraconate/glutaconate','1-carboxyethylleucine','succinate','creatinine','N-acetylarginine','alanine','hypoxanthine','inosine','cortisol'] if x in R.index]][['case_unadj','p_unadj','case_adj','p_adj','lac_beta','p_lac']].round(4).to_string())
pd.concat(out).to_csv(OUT+'effort_adjusted_response.csv')
