import pandas as pd, numpy as np
F='/mnt/project-files/biomap/BioMapAI/data/indenpendent_cohorts/jci.insight.2022.xlsx'
OUT='/mnt/project-files/mecfs-research/pem/'
def load(sheet='ScaledImpDataZeroDrug&Tobacco'):
    d=pd.read_excel(F,sheet_name=sheet,header=None)
    fm=d.iloc[8:,:15].copy(); fm.columns=d.iloc[7,:15].values
    fm=fm.reset_index(drop=True)
    gm=d.iloc[:7,16:].T.copy(); gm.columns=['ID','Phenotype','Timepoint','Sex','AgeBin','BMIBin','Bell']
    gm=gm.reset_index(drop=True)
    X=d.iloc[8:,16:].astype(float).T.reset_index(drop=True)
    X.columns=fm.BIOCHEMICAL.astype(str).values
    return fm,gm,X
if __name__=='__main__':
    fm,gm,X=load()
    fm.to_pickle(OUT+'fm.pkl'); gm.to_pickle(OUT+'gm.pkl'); X.to_pickle(OUT+'X_scaled.pkl')
    fo,go,Xo=load('Original Dataset'); Xo.to_pickle(OUT+'X_orig.pkl'); go.to_pickle(OUT+'gm_orig.pkl')
    print(X.shape, Xo.shape, (go.ID.values==gm.ID.values).all(), (go.Timepoint.values==gm.Timepoint.values).all())
    print(gm.groupby(['Phenotype','Sex']).ID.nunique())
    print(gm.groupby(['Phenotype','Timepoint']).size().unstack())
    c=gm.groupby('ID').size(); print(c.value_counts())
    print(gm.drop_duplicates('ID').groupby('Phenotype').Bell.describe())
    print(gm.drop_duplicates('ID').groupby(['Phenotype']).AgeBin.value_counts().unstack())
    print('dup names',X.columns.duplicated().sum(), 'missing orig frac', Xo.isna().mean().describe())
    print(fm['SUPER PATHWAY'].value_counts())
    # Bell constant within person?
    print('Bell varies within ID:', (gm.groupby('ID').Bell.nunique()>1).sum())
