# Stream a full GWAS summary file and keep (a) rows with p<5e-6 and (b) optionally all MHC rows (chr6 25-35 Mb, either build).
# usage: 03_fetch_exposure.py <id> <url> [mhc]
import sys, subprocess, gzip, io, math, os
eid,url=sys.argv[1],sys.argv[2]; mhc=len(sys.argv)>3
out='/tmp/work_gmr/exp/'
cmd=f"curl -sS --retry 5 '{url}'"
if url.endswith('.zip'): cmd=f"curl -sS --retry 5 -o /tmp/work_gmr/exp/{eid}.zip '{url}' && unzip -p /tmp/work_gmr/exp/{eid}.zip | grep -v '^#' ; rm -f /tmp/work_gmr/exp/{eid}.zip"
elif url.endswith('.gz') or url.endswith('.bgz'): cmd+=" | gzip -dc"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,text=True,bufsize=1<<20)
raw=p.stdout.readline().rstrip('\n')
sep='\t' if '\t' in raw else None
hdr=raw.split('\t') if sep else raw.split()
low=[h.lower() for h in hdr]
def pick(*names):
    for n in names:
        if n.lower() in low: return low.index(n.lower())
    return None
ix=dict(rs=pick('hm_rsid','rsid','snp','rsids','variant_id','SNPID','MarkerName'),chr=pick('hm_chrom','chromosome','chr','#chrom','CHR'),pos=pick('hm_pos','base_pair_location','pos','bp','POS'),
        ea=pick('hm_effect_allele','effect_allele','allele1','alt','ea','A1','Effect_allele'),oa=pick('hm_other_allele','other_allele','allele0','ref','oa','A2','A0','Other_allele'),
        eaf=pick('hm_effect_allele_frequency','effect_allele_frequency','a1freq','af_alt','eaf','Freq1','A1FREQ'),beta=pick('hm_beta','beta','BETA','b','BETA_SLEEPDURATION'),se=pick('standard_error','se','sebeta','SE','SE_SLEEPDURATION'),
        p=pick('p_value','p','pval','P_BOLT_LMM','P_BOLT_LMM_INF','P','P_SLEEPDURATION'),n=pick('n','N','sample_size'))
sys.stderr.write(f'{eid} header {hdr[:16]} -> {ix}\n')
HM3=set(open('/tmp/work_gmr/ref/hm3_rs.txt').read().split()); fh=gzip.open(out+eid+'.hm3.tsv.gz','wt',compresslevel=1)
fo=open(out+eid+'.tsv','w'); fm=open(out+eid+'.mhc.tsv','w') if mhc else None
cols=['rs','chr','pos','ea','oa','eaf','beta','se','p','n']
fo.write('\t'.join(cols)+'\n'); fh.write('\t'.join(cols)+'\n')

if fm: fm.write('\t'.join(cols)+'\n')
nr=0; nk=0
for line in p.stdout:
    f=line.rstrip('\n').split('\t') if sep else line.split()
    nr+=1
    try: pv=float(f[ix['p']])
    except: continue
    ch=f[ix['chr']].replace('chr','') if ix['chr'] is not None else ''
    try: ps=int(float(f[ix['pos']])) if ix['pos'] is not None else 0
    except: ps=0
    if ix['rs'] is not None and f[ix['rs']] in HM3:
        r_=[f[ix[c]] if ix[c] is not None else 'NA' for c in cols]; r_[1]=ch; fh.write('\t'.join(r_)+'\n')
    inm= fm is not None and ch=='6' and 25_000_000<=ps<=35_000_000
    if pv<5e-6 or inm:
        row=[f[ix[c]] if ix[c] is not None else 'NA' for c in cols]
        row[1]=ch
        s='\t'.join(row)+'\n'
        if pv<5e-6: fo.write(s); nk+=1
        if inm: fm.write(s)
fo.close(); fh.close(); 
if fm: fm.close()
sys.stderr.write(f'{eid} rows {nr} kept {nk}\n')
