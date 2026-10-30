"""Per-dataset case-control tests at gene level (SD units), covariate adjusted."""
import sys; sys.path.insert(0,'/tmp/work_proteome/scripts')
from common import *
allres={}
def run(name, ds, sel=None, covars=(), test='case', leak=False, extra=None):
    X,m,f=load(ds)
    if sel is not None: keep=sel(m); X=X[keep.values]; m=m[keep.values]
    if leak:
        L=pd.read_csv(D+f'{ds}_leak.csv',index_col=0).leak; m=m.join(L.rename('leak'),how='left'); covars=tuple(covars)+('leak',)
    if extra is not None: m=extra(m)
    G=gene_collapse(X,f)
    des=m[[test]+list(covars)].copy()
    for c in des.columns:
        if des[c].dtype==object: 
            dm=pd.get_dummies(des[c],prefix=c,drop_first=True).astype(float); des=des.drop(columns=c).join(dm)
    ok=des.notna().all(1); des=des[ok]; G=G[ok.values]
    r=ols_many(G,des,test)
    n1=int((des[test]==1).sum()); n0=int((des[test]==0).sum())
    r['n_case']=n1; r['n_ctrl']=n0; r['dataset']=name
    r.index.name='gene'; r.to_csv(R+f'cc_{name}.csv')
    print(f'{name}: n={n1}/{n0} genes={len(r)} p<0.05={int((r.p<0.05).sum())} q<0.05={int((r.q<0.05).sum())} q<0.1={int((r.q<0.1).sum())}')
    print('   top:', ', '.join(f"{g}({b:+.2f},q={q:.2g})" for g,b,q in r.sort_values('p').head(12)[['beta','q']].itertuples()))
    return r
run('hoel_serum_raw','hoel_serum',covars=('male','age','bmi'))
run('hoel_serum_leakadj','hoel_serum',covars=('male','age','bmi'),leak=True)
run('germain_D1PRE','germain_plasma',sel=lambda m:m.timepoint=='D1PRE',covars=('male','age','bmi_bin','site'))
run('germain_D1PRE_nobmi','germain_plasma',sel=lambda m:m.timepoint=='D1PRE',covars=('male','age','site'))
run('walitt_serum','walitt_serum',covars=('male','age'))
run('walitt_csf','walitt_csf',covars=('male','age'))
run('rydland_ev','rydland_ev')
run('giloteaux_ev_0h','giloteaux_ev',sel=lambda m:m.time=='0h',covars=('tmt',))
run('igg_serum','igg_serum',covars=('batch',))
run('igg_serum_MEonly','igg_serum',sel=lambda m:m.group!='PCS',covars=('batch',))
run('bragee_csf_POTS','bragee_csf',covars=('severity',),test='pots')
