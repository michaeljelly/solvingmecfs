# PRJNA940635 (UK ME/CFS Biobank; gDNA TCR-locus capture, sorted CD4 and CD8). Stream first 400k R1 reads/run, keep reads with a TRBJ core 20-mer (grep -F),
# assign TRBV/TRBJ by 11-mers (IMGT), extract productive CDR3 (C...F junction, VDJdb format).
import pandas as pd, numpy as np, subprocess, re, os, collections
W='/tmp/work_scm/tcr'; K=11; NREADS=150000
def rc(s): return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def readfa(f):
    out={}; name=None
    for l in open(f):
        if l.startswith('>'): name=l.split('|')[1]; out[name]=''
        else: out[name]+=l.strip().upper().replace('.','')
    return out
V=readfa(f'{W}/TRBV.fasta'); J=readfa(f'{W}/TRBJ.fasta')
def kidx(ref):
    d=collections.defaultdict(set)
    for g,s in ref.items():
        for i in range(len(s)-K+1): d[s[i:i+K]].add(g.split('*')[0])
    return d
VI=kidx(V); JI=kidx(J)
jpats=[l.strip() for l in open(f'{W}/trbj_pats.txt') if l.strip()]
fw=set(); 
for k,s in J.items():
    core=s[-30:-10] if len(s)>=40 else s[-25:-5]; fw.add(core)
CODON={a+b+c:aa for (a,b,c),aa in zip([(x,y,z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
def tr(s): return ''.join(CODON.get(s[i:i+3],'X') for i in range(0,len(s)-2,3))
def parse(r):
    if not any(p in r for p in fw): r=rc(r)
    if not any(p in r for p in fw): return None
    vc=collections.Counter(); jc=collections.Counter()
    for i in range(len(r)-K+1):
        km=r[i:i+K]
        for g in VI.get(km,()): vc[g]+=1
        for g in JI.get(km,()): jc[g]+=1
    if not vc or vc.most_common(1)[0][1]<4 or not jc: return None
    v=vc.most_common(1)[0][0]; j=jc.most_common(1)[0][0]
    m=None
    for mm in re.finditer('TT[TC]GG[ACGT][ACGT][ACGT]GG',r): m=mm   # F G x G
    if not m: return dict(v=v,j=j,cdr3=None,prod=False)
    fpos=m.start()
    if fpos<30: return None
    aa=tr(r[fpos%3:fpos+3]); cpos=[x.start() for x in re.finditer('C',aa)]
    if not cpos: return dict(v=v,j=j,cdr3=None,prod=False)
    for c in reversed(cpos):
        cdr3=aa[c:]
        if 8<=len(cdr3)<=25: break
    prod=('*' not in cdr3) and ('X' not in cdr3) and cdr3.endswith('F') and 8<=len(cdr3)<=25 and cdr3.startswith('C')
    return dict(v=v,j=j,cdr3=cdr3,prod=prod)
runs=pd.read_csv(f'{W}/runs.tsv',sep='\t'); os.makedirs(f'{W}/parsed',exist_ok=True)
open(f'{W}/pats_all.txt','w').write('\n'.join(jpats)+'\n')
from concurrent.futures import ThreadPoolExecutor
def fetch(r):
    url='https://'+r.fastq_ftp.split(';')[0]
    cmd=f"curl -s --retry 3 -m 900 '{url}' | gunzip 2>/dev/null | head -{NREADS*4} | awk 'NR%4==2' | grep -F -f {W}/pats_all.txt"
    return r, subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout.split()
todo=[r for _,r in runs.iterrows() if not os.path.exists(f'{W}/parsed/{r.run_accession}.parquet')]
with ThreadPoolExecutor(4) as ex:
    for r,seqs in ex.map(fetch,todo):
        rows=[x for x in (parse(s) for s in seqs) if x]
        d=pd.DataFrame(rows) if rows else pd.DataFrame(columns=['v','j','cdr3','prod']); d['run']=r.run_accession; d['sample']=r.sample_accession
        if len(seqs)>0: d.to_parquet(f'{W}/parsed/{r.run_accession}.parquet')
        print(r.run_accession,len(seqs),len(d),int(d['prod'].sum()) if len(d) else 0,flush=True)
print('DONE')
