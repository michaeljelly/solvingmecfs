# Lift QC'd DecodeME variants (GRCh38) to GRCh37 and map to 1000G EUR (MAGMA ref) rsIDs by chr:pos + alleles
import pandas as pd, numpy as np
from pyliftover import LiftOver
W='/tmp/work_decodeme/'
lo=LiftOver(W+'ref/hg38ToHg19.over.chain.gz')
ids=pd.read_csv(W+'ss/gwas_qced.var.gz')['ID']
sp=ids.str.split(':',expand=True); sp.columns=['chr','pos','a0','a1']
sp['ID']=ids; sp['chr']=sp.chr.astype(int); sp['pos']=sp.pos.astype(int)
out=[]
bim_it=pd.read_csv(W+'magma/g1000_eur.bim',sep='\t',header=None,names=['chr','rs','cm','pos','b1','b2'],usecols=[0,1,3,4,5],chunksize=3_000_000)
bim=pd.concat([b for b in bim_it]); print('bim',len(bim),flush=True)
for c in range(1,23):
    s=sp[sp.chr==c].copy()
    p19=np.full(len(s),-1,dtype=np.int64); cn='chr%d'%c
    for i,p in enumerate(s.pos.values):
        r=lo.convert_coordinate(cn,int(p)-1)
        if r and len(r)==1 and r[0][0]==cn: p19[i]=r[0][1]+1
    s['pos19']=p19; s=s[s.pos19>0]
    b=bim[bim.chr==c]
    m=s.merge(b,left_on='pos19',right_on='pos',suffixes=('','_b'))
    ok=((m.a0==m.b1)&(m.a1==m.b2))|((m.a0==m.b2)&(m.a1==m.b1))
    m=m[ok].drop_duplicates('ID').drop_duplicates('rs')
    print(c,len(s),len(m),flush=True)
    out.append(m[['ID','rs','chr','pos19']])
o=pd.concat(out); o.to_csv(W+'pq/id_rs_map.tsv.gz',sep='\t',index=False)
print('mapped',len(o))
