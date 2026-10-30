# Convert GEO series-matrix betas to parquet + parse sample metadata
import pandas as pd, numpy as np, gzip, re, sys
W='/tmp/work_scm/meth'
def parse(gse):
    f=f'{W}/{gse}_series_matrix.txt.gz'; meta={}; start=None
    with gzip.open(f,'rt') as fh:
        for i,l in enumerate(fh):
            if l.startswith('!Sample_'):
                p=l.rstrip('\n').split('\t'); k=p[0][8:]; v=[x.strip('"') for x in p[1:]]
                while k in meta: k=k+'_'
                meta[k]=v
            if l.startswith('"ID_REF"'): start=i; break
    M=pd.DataFrame(meta)
    B=pd.read_csv(f,sep='\t',skiprows=start,index_col=0,comment='!',dtype={0:str},low_memory=False)
    B=B.apply(pd.to_numeric,errors='coerce').astype(np.float32)
    B.index=B.index.astype(str)
    return M,B
for gse in sys.argv[1:]:
    M,B=parse(gse)
    M.to_csv(f'{W}/{gse}_meta_raw.csv',index=False)
    B.to_parquet(f'{W}/{gse}_beta.parquet')
    print(gse,B.shape,B.isna().mean().mean(),flush=True)
