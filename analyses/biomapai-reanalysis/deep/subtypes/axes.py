"""Candidate biological axes (per visit): lipolysis/orthostatic, inflammatory kynurenine, gut dysbiosis."""
from sub_common import *
QC='/mnt/project-files/biomap/BioMapAI/data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv'
def z(v): return (v-v.mean())/v.std()
def build():
    q=pd.read_csv(QC,index_col=0); samp=[c for c in q.columns if '_tp' in c]
    sp=q.SUB_PATHWAY; nm=q.CHEMICAL_NAME
    L=np.log2(q[samp].astype(float)).T; L.columns=nm.values
    sub=pd.Series(sp.values,index=nm.values)
    ffa=[c for c in L.columns if sub.get(c) in ('Long Chain Polyunsaturated Fatty Acid (n3 and n6)','Long Chain Saturated Fatty Acid','Long Chain Monounsaturated Fatty Acid')]
    acyl=[c for c in L.columns if str(sub.get(c)).startswith('Fatty Acid Metabolism (Acyl Carnitine')]
    bhb=[c for c in L.columns if c.startswith('3-hydroxybutyrate')][0]
    A=pd.DataFrame(index=L.index)
    A['FFA']=L[ffa].apply(z).mean(1); A['acylcarnitine']=L[acyl].apply(z).mean(1); A['BHB']=L[bhb]
    A['KT_ratio']=L['kynurenine']-L['tryptophan']  # log2(kyn/trp)
    im=load('immune')
    ifng=['% CD4+ RO+ total INFg+ d1 P8','% CD4+ RO+ total INFg+ d6','% CD8+ total INFg+ d6','% CD8+ nonMAIT total INFg+ d1']
    ifng=[c for c in ifng if c in im.columns]
    I=pd.DataFrame(index=im.index)
    I['IFNg_T']=int_transform(im[ifng]).mean(1)
    I['MAIT_CD8']=im['%CD8+ MAIT d0']; I['gd_T']=im['Total gd Cells (gdTCR CD3hi + CD3lo)']
    s=load('specie')
    bp=[c for c in s.columns if any(k in c for k in ['Faecalibacterium','Roseburia','Coprococcus','Agathobacter','Eubacterium_rectale','Anaerostipes','Butyricicoccus','Anaerobutyricum'])]
    G=pd.DataFrame(index=s.index); G['butyrate_producers']=np.log10(s[bp].sum(1)+1e-4); G['R_gnavus']=np.log10(s['Ruminococcus_gnavus']+1e-5)
    D=meta.join(score,how='left').join(A,how='left').join(I,how='left').join(G,how='left')
    return D,dict(ffa=len(ffa),acyl=len(acyl),bhb=bhb,ifng=ifng,bp=bp)
COMP=['FFA','acylcarnitine','BHB','KT_ratio','IFNg_T','MAIT_CD8','gd_T','butyrate_producers','R_gnavus']
def axis_scores(D,ref):
    """axis = mean of z-scored components, z computed on reference rows (tp1 patients)"""
    Zc={c:(D[c]-D.loc[ref,c].mean())/D.loc[ref,c].std() for c in COMP}
    # rank-normalise heavy tailed ones before averaging
    D=D.copy()
    D['axis_lipolysis']=(Zc['FFA']+Zc['acylcarnitine']+Zc['BHB'])/3
    D['axis_kyn_inflam']=(Zc['KT_ratio']+Zc['IFNg_T'])/2
    D['axis_dysbiosis']=(Zc['R_gnavus']-Zc['butyrate_producers'])/2
    return D
