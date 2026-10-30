import pandas as pd, numpy as np, sys
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
O='/tmp/work_geo/results/'; import os; os.makedirs(O,exist_ok=True)
meta=pd.read_csv('/tmp/work_geo/GSE214283_meta.csv').set_index('geo_accession')
meta['case']=(meta['disease state']=='case').astype(int); meta['subject']=meta['title'].str[:8]; meta['day']=meta['collection day']
meta['male']=(meta['Sex']=='M').astype(int); meta['age']=meta['age bin'].astype(float); meta['batch']=meta['chromium batch']
G=[];S=[];P=[];W=[]
for ct in ['ALL','CD4T_other','CD8T','NK','B','cMono','ncMono']:
    c=pd.read_parquet(f'/tmp/work_geo/pb/pb_{ct}.parquet'); c=c[[x for x in meta.index if x in c.columns]]
    c=c.loc[:,c.sum(0)>2e5]  # drop very low depth pseudobulks
    lm=filter_expr(logcpm(c),c,min_cpm=0.5,min_frac=0.25)
    d1=meta[(meta.day=='D1')&meta.index.isin(lm.columns)]; d2=meta[(meta.day=='D2')&meta.index.isin(lm.columns)]
    for lab,mm in [('D1_baseline',d1),('D2_24h_postCPET',d2)]:
        g,s,sc,miss=run_case_control('GSE214283_'+ct,lm,mm,covars=('male','age'),extra_label='_'+lab)
        G.append(g.assign(celltype=ct));S.append(s.assign(celltype=ct))
    r,w=run_paired('GSE214283_'+ct,lm,meta[meta.index.isin(lm.columns)],'subject','day','D1','D2',covars=('male','age'))
    P.append(r.assign(celltype=ct)); W.append(w.assign(celltype=ct))
    print(ct,lm.shape,'missing targets:',len(miss),flush=True)
G=pd.concat(G);S=pd.concat(S);P=pd.concat(P);W=pd.concat(W)
G.to_csv(O+'GSE214283_genes_case_control.csv'); S.to_csv(O+'GSE214283_sets_case_control.csv'); P.to_csv(O+'GSE214283_change_case_vs_ctrl.csv'); W.to_csv(O+'GSE214283_change_within_group.csv',index=False)
# cell proportions
cp=pd.read_csv('/tmp/work_geo/pb/celltype_counts.csv').set_index('gsm')
types=['CD4T_other','CD8T','NK','B','cMono','ncMono','pDC']
prop=cp[types].div(cp['n_cells'],axis=0); lp=np.log(prop+1e-3).T
rows=[]
for lab,mm in [('D1',meta[meta.day=='D1']),('D2',meta[meta.day=='D2'])]:
    t=test_matrix(lp,mm,('male','age')); t['contrast']='case_vs_control_'+lab; rows.append(t)
dl,sm_=paired_change(lp,meta.rename_axis('sample'),'subject','day','D1','D2'); t=test_matrix(dl,sm_,('male','age')); t['contrast']='change_case_minus_ctrl'; rows.append(t)
pr=pd.concat(rows); pr.to_csv(O+'GSE214283_celltype_props.csv'); print(pr.round(3).to_string())
print('\nSETS D1 (q<0.25 or p<0.05):'); print(S[(S.p<0.05)].round(3)[['celltype','contrast','beta','d','p','q','n_case','n_ctrl']].to_string())
print('\nSETS change case-minus-ctrl p<0.05:'); print(P[P.is_set&(P.p<0.05)].round(3)[['celltype','beta','d','p','q','n_case','n_ctrl']].to_string())
print('\nGENES q<0.1:'); print(G[G.q<0.1].round(3)[['celltype','contrast','beta','p','q']].to_string())
print('\nCHANGE genes q<0.1:'); print(P[(~P.is_set)&(P.q<0.1)].round(3)[['celltype','beta','p','q']].to_string())
