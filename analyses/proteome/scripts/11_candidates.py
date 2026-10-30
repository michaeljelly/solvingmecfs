"""Named candidate proteins for the requested angles: effect in every dataset + random-effects meta across ME/CFS blood cohorts."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
C=['ACE','ACE2','AGT','REN','EDN1','ADM','NPPA','NPPB','NPPC','CORIN','EPO','HAMP','TFRC','FTL','TF','HPX','HP','VWF','PF4','PPBP','SELP','THBS1','FGA','FGB','FGG','F13A1','SERPINE1','PLAT','ICAM1','VCAM1','SELE','ANGPT1','ANGPT2','TEK','VEGFA','KDR','FLT1','ESM1','MB','CKM','MYL3','NEFL','GFAP','UCHL1','IDO1','KYNU','KMO','TGFB1','LTBP2','GDF15','FGF21','CXCL10','CXCL9','LGALS3BP','SIGLEC1','ISG15','TNFSF13B','CXCL13','CD27','IL6','CRP','C3','C4A','C5','CFB','CFD','CFH','C1QA','C7','ADH4','GSTA1','KRT18','LEP','FABP4','ADIPOQ','ITGAV','NOTCH3','MYOC','LUM','COMP','FAP','MMP2','MMP3','CD93','SPINK9']
DS={'hoel':'hoel_serum_leakadj','germain':'germain_D1PRE','walitt':'walitt_serum','ukb':'ukb_olink','walitt_csf':'walitt_csf','rydland_ev':'rydland_ev','igg':'igg_serum','lc_uk':'lc_uk_explore','bragee_POTS':'bragee_csf_POTS'}
T={k:pd.read_csv(R+f'cc_{v}.csv',index_col=0) for k,v in DS.items()}
rows=[]
for g in C:
    r={'gene':g}
    for k,t in T.items():
        if g in t.index: r[f'd_{k}']=t.loc[g,'beta']; r[f'p_{k}']=t.loc[g,'p']
    d=[r.get(f'd_{k}') for k in ['hoel','germain','walitt','ukb']]; se=[T[k].loc[g,'se'] if g in T[k].index else np.nan for k in ['hoel','germain','walitt','ukb']]
    mm=meta_re([x if x is not None else np.nan for x in d],se); r.update({'meta_k':mm.get('k'),'meta_d':mm.get('re'),'meta_p':mm.get('re_p'),'I2':mm.get('I2')})
    rows.append(r)
X=pd.DataFrame(rows).set_index('gene'); X['meta_q_within_list']=bh(X.meta_p.values); X.to_csv(R+'candidate_proteins.csv')
pd.set_option('display.width',260)
print(X[[c for c in X.columns if c.startswith('d_')]+['meta_d','meta_p','meta_q_within_list']].round(2).to_string())
