import pandas as pd, numpy as np, sys
from scipy import stats
sys.path.insert(0,'/tmp/work_geo/scripts'); from common import *
O='/tmp/work_geo/results/'
def load(f,**flt):
    d=pd.read_csv(O+f,index_col=0)
    for k,v in flt.items(): d=d[d[k]==v]
    return d
src={
 'GSE214283 PBMC scRNA-pseudobulk D1 (30/28)':load('GSE214283_genes_case_control.csv',celltype='ALL',contrast='case_vs_control_D1_baseline'),
 'GSE227375 PBMC bulk T0 (23/25)':load('GSE227375_genes_case_control.csv',contrast='case_vs_control_T0_baseline'),
 'GSE98139 adolescent WB (28/19)':load('small_genes_case_control.csv',dataset='GSE98139_adolescent_WB'),
 'GSE128078 WB day1 (13/10)':load('small_genes_case_control.csv',dataset='GSE128078_WB'),
 'GSE293840 plasma cfRNA (93/75)':load('GSE293840_genes_case_control.csv'),
 'GSE293840 cfRNA global-adj':load('GSE293840_genes_case_control_globadj.csv'),
 'PAD000026 serum SomaScan (50/29)':load('PAD000026_proteins_case_control.csv'),
 'PAD000026 leak-adj':load('PAD000026_proteins_case_control_leakadj.csv'),
}
tw=pd.read_csv(O+'small_paired.csv',index_col=0); tw=tw[tw.dataset=='GSE16059_twins_all_44pairs']
src['GSE16059 twins affected-cotwin (44 pairs)']=tw
ssrc={k.replace('genes','sets'):None for k in src}
sets={
 'GSE214283 PBMC scRNA-pseudobulk D1 (30/28)':load('GSE214283_sets_case_control.csv',celltype='ALL',contrast='case_vs_control_D1_baseline'),
 'GSE227375 PBMC bulk T0 (23/25)':load('GSE227375_sets_case_control.csv',contrast='case_vs_control_T0_baseline'),
 'GSE98139 adolescent WB (28/19)':load('small_sets_case_control.csv',dataset='GSE98139_adolescent_WB'),
 'GSE128078 WB day1 (13/10)':load('small_sets_case_control.csv',dataset='GSE128078_WB'),
 'GSE293840 plasma cfRNA (93/75)':load('GSE293840_sets_case_control.csv'),
 'GSE293840 cfRNA global-adj':load('GSE293840_sets_case_control_globadj.csv'),
 'PAD000026 serum SomaScan (50/29)':load('PAD000026_sets_case_control.csv'),
 'PAD000026 leak-adj':load('PAD000026_sets_case_control_leakadj.csv'),
 'GSE16059 twins affected-cotwin (44 pairs)':tw[tw.is_set],
}
def fmt(r): 
    if r is None or pd.isna(r.get('d',np.nan)): return 'na'
    star='**' if r['q']<0.05 else ('*' if r['p']<0.05 else '')
    return f"{r['d']:+.2f}{star}"
def table(S,feats):
    rows=[]
    for f in feats:
        row={'feature':f}
        for k,d in S.items(): row[k]=fmt(d.loc[f]) if f in d.index else 'nd'
        rows.append(row)
    return pd.DataFrame(rows).set_index('feature')
T1=table(sets,list(GENESETS.keys())); T1.to_csv(O+'SUMMARY_sets_baseline_by_dataset.csv')
keyg=['IL17A','IL17F','RORC','IL23R','CCR6','KLRB1','IL22','PDCD1','TOX','HAVCR2','IDO1','IDO2','TDO2','KMO','KYNU','HAAO','AHR','CYP1A1','CYP1B1','AHRR','CPT1A','CPT2','ACADVL','HMGCS2','OXCT1','PDK4','PPARA','ADRB2','ADRB1','ADRA2A','COMT','NR4A1','IFI27','IFI44L','ISG15','SIGLEC1','MS4A1','CD79A','FCRL1','NCR1','KLRF1','GNLY']+GENESETS['DecodeME_named']
T2=table(src,keyg); T2.to_csv(O+'SUMMARY_genes_baseline_by_dataset.csv')
# Stouffer meta across independent blood transcriptome cohorts (not cfRNA/proteome): 214283, 227375, 98139, 128078, twins
cohorts=['GSE214283 PBMC scRNA-pseudobulk D1 (30/28)','GSE227375 PBMC bulk T0 (23/25)','GSE98139 adolescent WB (28/19)','GSE128078 WB day1 (13/10)','GSE16059 twins affected-cotwin (44 pairs)']
N={cohorts[0]:58,cohorts[1]:48,cohorts[2]:47,cohorts[3]:23,cohorts[4]:44}
def meta(S,feats):
    rows=[]
    for f in feats:
        zs=[];ws=[];dirs=[]
        for k in cohorts:
            d=S[k]
            if f in d.index and pd.notna(d.loc[f,'p']):
                z=np.sign(d.loc[f,'d'])*stats.norm.isf(d.loc[f,'p']/2); zs.append(z); ws.append(np.sqrt(N[k])); dirs.append(np.sign(d.loc[f,'d']))
        if len(zs)>=3:
            Z=np.dot(ws,zs)/np.sqrt(np.sum(np.square(ws))); rows.append(dict(feature=f,k=len(zs),n_up=int(np.sum(np.array(dirs)>0)),n_down=int(np.sum(np.array(dirs)<0)),stouffer_z=Z,p=2*stats.norm.sf(abs(Z))))
    r=pd.DataFrame(rows).set_index('feature'); r['q']=bh(r.p); return r
M1=meta(sets,list(GENESETS.keys())); M2=meta(src,TARGET_GENES)
M1.to_csv(O+'META_sets_blood_transcriptome.csv'); M2.to_csv(O+'META_genes_blood_transcriptome.csv')
pd.set_option('display.width',250)
short={k:k.split(' ')[0]+('_adj' if 'adj' in k else '') for k in list(src)}
print(T1.rename(columns=short).to_string()); print(); print(T2.rename(columns=short).to_string()); print(); print(M1.round(4).to_string()); print(); print(M2.sort_values('p').head(20).round(4).to_string())
