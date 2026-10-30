"""Germain 2025 SomaScan: protein-level exertion response (change scores) patients vs controls."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
X,m,f=load('germain_plasma'); G=gene_collapse(X,f)
G['cor_id']=m.cor_id.values; G['tp']=m.timepoint.values
dm=m.drop_duplicates('cor_id').set_index('cor_id')
res={}
for a,b,lab in [('D1POST','D1PRE','acute_D1'),('D2PRE','D1PRE','recovery_24h'),('D2POST','D2PRE','acute_D2'),('D2POST','D1POST','post2_vs_post1')]:
    A=G[G.tp==a].set_index('cor_id').drop(columns='tp'); B=G[G.tp==b].set_index('cor_id').drop(columns='tp')
    ids=A.index.intersection(B.index); C=A.loc[ids]-B.loc[ids]
    des=dm.loc[ids,['case','male','age']].join(pd.get_dummies(dm.loc[ids,'site'].astype(str),prefix='site',drop_first=True).astype(float))
    r=ols_many(C,des,'case'); r.to_csv(R+f'germain_exertion_{lab}.csv')
    # within-group mean change (paired t)
    for grp in [0,1]:
        Cg=C[des.case==grp]; t=Cg.mean()/(Cg.std()/np.sqrt(len(Cg))); r[f'mean_change_{grp}']=Cg.mean(); r[f'p_change_{grp}']=2*stats.t.sf(t.abs(),len(Cg)-1)
    r.to_csv(R+f'germain_exertion_{lab}.csv')
    print(lab,'n',len(ids),'case-vs-ctrl diff in change: p<0.05',(r.p<0.05).sum(),'q<0.05',(r.q<0.05).sum(),'| HC changed q<0.05',(bh(r.p_change_0)<0.05).sum(),'ME changed q<0.05',(bh(r.p_change_1)<0.05).sum())
    print('   top:',', '.join(f"{g}({b:+.2f},p={p:.1g})" for g,b,p in r.sort_values('p').head(10)[['beta','p']].itertuples()))
