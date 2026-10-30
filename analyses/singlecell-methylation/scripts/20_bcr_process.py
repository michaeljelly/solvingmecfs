# PRJEB83659 (UK, BCR heavy chain amplicons, MiSeq). Stream first N R1 reads per run, assign IGHV/IGHJ by k-mers (IMGT), isotype from C primer, CDR3 by motif, SHM proxy.
import pandas as pd, numpy as np, requests, gzip, io, re, os, sys, collections
W='/tmp/work_scm/bcr'; N=10000; K=11
def rc(s): return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def readfa(f):
    out={}; name=None
    for l in open(f):
        if l.startswith('>'): name=l.split('|')[1]; out[name]=''
        else: out[name]+=l.strip().upper().replace('.','')
    return out
V=readfa(f'{W}/IGHV.fasta'); J=readfa(f'{W}/IGHJ.fasta')
V={k:v for k,v in V.items() if len(v)>250}
def kidx(ref):
    d=collections.defaultdict(set)
    for g,s in ref.items():
        gene=g.split('*')[0]
        for i in range(len(s)-K+1): d[s[i:i+K]].add(gene)
    return d
VI=kidx(V); JI=kidx(J)
Vk={g.split('*')[0]:set() for g in V}
for g,s in V.items(): Vk[g.split('*')[0]].update(s[i:i+K] for i in range(len(s)-K+1))
CODON={a+b+c:aa for (a,b,c),aa in zip([(x,y,z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
def tr(s): return ''.join(CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
PRIM={'GGGCGGATGCACTCC':'IgM','GGCCCTTGGTGGAGG':'IgG'}
def qtrim(seq,qual,thr=20,win=10):
    q=np.frombuffer(qual.encode(),dtype=np.uint8)-33
    if len(q)<win: return seq
    m=np.convolve(q,np.ones(win)/win,'valid'); bad=np.where(m<thr)[0]
    return seq[:bad[0]+win//2] if len(bad) else seq
def process(read,qual):
    iso=PRIM.get(read[:15],'other')
    s=qtrim(read,qual)
    if len(s)<150: return None
    r=rc(s)
    vc=collections.Counter(); jc=collections.Counter()
    for i in range(len(r)-K+1):
        km=r[i:i+K]
        for g in VI.get(km,()): vc[g]+=1
        for g in JI.get(km,()): jc[g]+=1
    if not vc or vc.most_common(1)[0][1]<10: return None
    top=vc.most_common(3); v=top[0][0]; vn=top[0][1]
    ties=sorted(g for g,c in top if c==vn)
    j=jc.most_common(1)[0][0] if jc else None
    # V region = read part before J match; SHM proxy = fraction of read V kmers absent from assigned germline
    m=re.search('TGGGG[CG][CA][AG]GG',r)  # WGQG / WGRG
    cdr3=None
    if m:
        wpos=m.start()
        for fr in range(3):
            if (wpos-fr)%3==0:
                aa=tr(r[fr:wpos]); mm=list(re.finditer('Y[YFHC]C',aa))
                if mm:
                    c=aa[mm[-1].end():]
                    if 3<=len(c)<=35 and '*' not in c and 'X' not in c: cdr3=c
                break
        vreg=r[:max(0,wpos-40)]
    else: vreg=r[:max(0,len(r)-80)]
    kms=[vreg[i:i+K] for i in range(0,len(vreg)-K+1)]
    shm=np.mean([km not in Vk[v] for km in kms]) if len(kms)>=40 else np.nan
    return dict(iso=iso,v=v,vtie='/'.join(ties),j=j,cdr3=cdr3,shm=shm,vlen=len(vreg))
runs=pd.read_csv(f'{W}/runs_ftp.tsv',sep='\t')
runs=runs.sort_values('run_accession').groupby('library_name').head(1)   # one run per library (A and B = replicate PCRs)
os.makedirs(f'{W}/parsed',exist_ok=True)
import subprocess
from concurrent.futures import ThreadPoolExecutor
def fetch(r):
    url='https://'+r.fastq_ftp.split(';')[0]
    txt=subprocess.run(f"curl -s --retry 3 -m 600 '{url}' | gunzip 2>/dev/null | head -{N*4}",shell=True,capture_output=True,text=True).stdout.split('\n')
    return r, txt
todo=[r for _,r in runs.iterrows() if not os.path.exists(f'{W}/parsed/{r.run_accession}.parquet')]
with ThreadPoolExecutor(4) as ex:
    for r,txt in ex.map(fetch,todo):
        rows=[]
        for k in range(0,len(txt)-3,4):
            x=process(txt[k+1].strip(),txt[k+3].strip())
            if x: rows.append(x)
        d=pd.DataFrame(rows); d['run']=r.run_accession; d['sample']=r.sample_alias; d['lib']=r.library_name
        d.to_parquet(f'{W}/parsed/{r.run_accession}.parquet'); print(r.run_accession,r.library_name,len(txt)//4,len(d),flush=True)
print('DONE')
