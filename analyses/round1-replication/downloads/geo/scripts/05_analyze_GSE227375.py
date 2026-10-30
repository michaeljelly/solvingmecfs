import pandas as pd, numpy as np, sys
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
O='/tmp/work_geo/results/'
m=pd.read_csv('/tmp/work_geo/GSE227375_meta.csv',index_col=0); m['male']=(m.sex=='male').astype(int); m['plat']=(m.platform=='GPL21290').astype(int)
c=pd.read_parquet('/tmp/work_geo/GSE227375_counts.parquet')[m.index]
lm=filter_expr(logcpm(c),c,min_cpm=0.5,min_frac=0.25)
G=[];S=[];P=[];W=[]
for t,lab in [(0,'T0_baseline'),(1,'T1_maxexertion'),(2,'T2_recovery')]:
    mm=m[m.time==t]; g,s,sc,miss=run_case_control('GSE227375',lm,mm,covars=('male','plat'),extra_label='_'+lab); G.append(g);S.append(s)
for t0,t1 in [(0,1),(0,2),(1,2)]:
    r,w=run_paired('GSE227375',lm,m,'subject','time',t0,t1,covars=('male',),label=f'T{t0}->T{t1}'); P.append(r);W.append(w)
G=pd.concat(G);S=pd.concat(S);P=pd.concat(P);W=pd.concat(W)
G.to_csv(O+'GSE227375_genes_case_control.csv'); S.to_csv(O+'GSE227375_sets_case_control.csv'); P.to_csv(O+'GSE227375_change_case_vs_ctrl.csv'); W.to_csv(O+'GSE227375_change_within_group.csv',index=False)
print('missing',miss)
print(S.round(3)[['contrast','beta','d','p','q','n_case','n_ctrl']].to_string())
print('\nCHANGE sets p<0.1'); print(P[P.is_set&(P.p<0.1)].round(3)[['contrast','beta','d','p','q','n_case','n_ctrl']].to_string())
print('\nGENES q<0.1'); print(G[G.q<0.1].round(3)[['contrast','beta','d','p','q']].to_string())
print('\nCHANGE genes q<0.1'); print(P[(~P.is_set)&(P.q<0.1)].round(3)[['contrast','beta','p','q']].to_string())
