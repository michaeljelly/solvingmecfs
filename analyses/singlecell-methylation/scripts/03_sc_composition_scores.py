# Composition (CLR, compositional) + marker-gated subsets + pre-specified per-cell-type scores. Baseline (D1) case vs control; 24h change (D2-D1) case vs control.
import pandas as pd, numpy as np, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
R='/mnt/project-files/mecfs-research/deep2/singlecell_methyl/results'
CE=pd.read_parquet('/tmp/work_scm/cells_all.parquet')
meta=pd.read_csv(f'{R}/sc_sample_meta.csv').set_index('geo_accession')
meta['case']=(meta['disease state']=='case').astype(int); meta['age']=meta['age bin'].astype(float)
CE=CE[~CE.type.isin(['Doublet','Other','Erythroid'])]
# ---- marker-gated subsets (fractions within parent) ----
T=CE.type
isT=T.str.startswith(('CD4','CD8','Treg','MAIT','Th','gdT','NKT','Cycling_T','T_other','Tfh'))
sub={}
def frac(name,parent_mask,child_mask):
    sub[name]=(parent_mask,child_mask&parent_mask)
g=lambda x:CE['g_'+x].values>0
cd8=T.isin(['CD8_Tem','CD8_Temra']); cd4m=T.isin(['CD4_Tem','CD4_Tem_PD1']); nk=T.str.startswith('NK_'); B=T.isin(['B_naive','B_memory','B_ABC']); mono=T.str.startswith('Mono_')
frac('gdT_of_T',isT,g('TRDC')&~g('SLC4A10'))
frac('Th17like_CCR6+KLRB1+_of_CD4mem',cd4m,g('CCR6')&g('KLRB1'))
frac('RORC+_of_CD4mem',cd4m,g('RORC'))
frac('Th1like_TBX21+_of_CD4mem',cd4m,g('TBX21'))
frac('Treg_FOXP3+_of_CD4',T.str.startswith('CD4')|(T=='Treg'),g('FOXP3'))
frac('NKbright_FCGR3A-GZMK+_of_NK',nk,~g('FCGR3A')&g('GZMK'))
frac('NKadaptive_KLRC2+FCER1G-_of_NK',nk,g('KLRC2')&~g('FCER1G'))
frac('NKG2C_KLRC2+_of_NK',nk,g('KLRC2'))
frac('CX3CR1+_of_CD8mem',cd8,g('CX3CR1'))
frac('KLRG1+CX3CR1+_of_CD8mem',cd8,g('KLRG1')&g('CX3CR1'))
frac('ZNF683+_of_CD8mem',cd8,g('ZNF683'))
frac('B3GAT1(CD57)+_of_CD8mem',cd8,g('B3GAT1'))
frac('Exhausted_PDCD1+TOX+_of_CD8mem',cd8,g('PDCD1')&g('TOX'))
frac('PDCD1+_of_CD8mem',cd8,g('PDCD1'))
frac('HAVCR2+_of_CD8mem',cd8,g('HAVCR2'))
frac('ABC_TBX21+orITGAX+_of_B',B,g('TBX21')|g('ITGAX'))
frac('FCRL5+_of_B',B,g('FCRL5'))
frac('CR2(CD21)-_of_memB',T.isin(['B_memory','B_ABC']),~g('CR2'))
frac('ABC_celltypist_of_B',B,T.values=='B_ABC')
frac('Plasmablast_of_B',B|(T=='Plasmablast'),T.values=='Plasmablast')
frac('PlateletRNA+_of_Mono',mono,g('PPBP')|g('PF4'))
frac('PlateletRNA+_of_CD4naive',T=='CD4_naive_CM',g('PPBP')|g('PF4'))
frac('Mono_intermediate_of_Mono',mono,T.values=='Mono_intermediate')
frac('Mono_nonclassical_of_Mono',mono,T.values=='Mono_nonclassical')
frac('Ki67+_of_T',isT,g('MKI67'))
frac('CD38+_of_CD8mem',cd8,g('CD38'))
frac('IFNG+_of_CD8mem',cd8,g('IFNG'))
rows=[]
for k,(pm,cm) in sub.items():
    df=pd.DataFrame({'gsm':CE.gsm.values,'p':np.asarray(pm),'c':np.asarray(cm)}).groupby('gsm')[['p','c']].sum()
    df['logit']=np.log((df.c+0.5)/(df.p-df.c+0.5)); df['frac']=df.c/df.p
    for gsm,r in df.iterrows(): rows.append(dict(gsm=gsm,feature='sub:'+k,value=r.logit,raw=r.frac,n_parent=r.p))
# ---- CLR composition of major types ----
cnt=pd.crosstab(CE.gsm,CE.type); keep=[c for c in cnt.columns if (cnt[c]>=3).mean()>0.5]
cc=cnt[keep]+0.5; clr=np.log(cc).sub(np.log(cc).mean(1),axis=0); prop=cnt[keep].div(cnt.sum(1),axis=0)
for c in keep:
    for gsm in cnt.index: rows.append(dict(gsm=gsm,feature='clr:'+c,value=clr.loc[gsm,c],raw=prop.loc[gsm,c],n_parent=cnt.loc[gsm].sum()))
# ---- pre-specified mean scores per cell type ----
scols=[c for c in CE.columns if c.startswith('s_')]+['pct_mt','pct_ribo','ngenes']
S=CE.groupby(['gsm','type'])[scols].mean(); N=CE.groupby(['gsm','type']).size()
PRE={'pDC':['s_ISG'],'cDC2':['s_ISG','s_MHCII'],'Mono_classical':['s_ISG','s_hypoxia','s_glycolysis','s_oxphos','s_FAO','s_inflam','s_platelet','s_MHCII','s_IEG_stress','s_heatshock','s_adrenergic','s_kynurenine','pct_mt','pct_ribo'],
     'Mono_intermediate':['s_ISG','s_hypoxia','s_inflam','s_platelet','s_MHCII','s_kynurenine'],'Mono_nonclassical':['s_ISG','s_hypoxia','s_inflam','s_platelet'],
     'CD8_Temra':['s_cyto_effector','s_exhaustion','s_ISG','pct_mt','s_oxphos','s_glycolysis'],'CD8_Tem':['s_cyto_effector','s_exhaustion','s_ISG','pct_mt'],
     'CD4_Tem':['s_Th17','s_exhaustion','s_ISG','pct_mt'],'MAIT':['s_ISG','s_cyto_effector','pct_mt'],'NK_CD56dim':['s_ISG','s_cyto_effector','s_adrenergic','pct_mt','s_oxphos'],
     'B_memory':['s_ABC','s_ISG','pct_mt'],'B_naive':['s_ISG','pct_mt'],'CD4_naive_CM':['s_ISG','pct_mt','pct_ribo','s_IEG_stress','s_heatshock','s_oxphos']}
for ct,feats in PRE.items():
    for f in feats:
        for gsm in meta.index:
            if (gsm,ct) in S.index and N[(gsm,ct)]>=5: rows.append(dict(gsm=gsm,feature=f'score:{ct}:{f}',value=S.loc[(gsm,ct),f],raw=S.loc[(gsm,ct),f],n_parent=N[(gsm,ct)]))
L=pd.DataFrame(rows).merge(meta[['person','day','case','Sex','age','chromium batch']],left_on='gsm',right_index=True)
L.to_csv(f'{R}/sc_features_long.csv',index=False)
# ---- tests ----
out=[]
for feat,d in L.groupby('feature'):
    d1=d[d.day=='D1'].copy()
    d1['z']=(d1.value-d1.value.mean())/d1.value.std()
    if d1.z.std()>0 and len(d1)>20:
        m=smf.ols('z ~ case + C(Sex) + age',d1).fit()
        out.append(dict(feature=feat,contrast='baseline_D1',n_case=int(d1.case.sum()),n_ctrl=int((1-d1.case).sum()),mean_case=d1.loc[d1.case==1,'raw'].mean(),mean_ctrl=d1.loc[d1.case==0,'raw'].mean(),d=m.params['case'],p=m.pvalues['case']))
    w=d.pivot_table(index='person',columns='day',values='value'); w=w.dropna()
    if 'D2' in w:
        dd=d[d.day=='D1'].set_index('person')[['case','Sex','age']].join((w.D2-w.D1).rename('chg'),how='inner')
        sd=w.D1.std(); dd['chg']=dd.chg/sd
        if dd.chg.std()>0:
            m=smf.ols('chg ~ case + C(Sex) + age',dd).fit()
            out.append(dict(feature=feat,contrast='change_D2minusD1',n_case=int(dd.case.sum()),n_ctrl=int((1-dd.case).sum()),mean_case=dd.loc[dd.case==1,'chg'].mean(),mean_ctrl=dd.loc[dd.case==0,'chg'].mean(),d=m.params['case'],p=m.pvalues['case']))
O=pd.DataFrame(out); O['family']=O.feature.str.split(':').str[0]
O['q']=O.groupby(['family','contrast']).p.transform(lambda p:multipletests(p,method='fdr_bh')[1])
O.to_csv(f'{R}/sc_composition_scores_tests.csv',index=False)
pd.set_option('display.width',250)
print(O.sort_values('p').head(45).round(4).to_string())
