# Remote tabix queries of GWAS Catalog harmonised files (bgzip + .tbi on EBI FTP over https)
import pysam, subprocess, pandas as pd, io, zlib
F='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics'
def acc_url(acc):
    n=int(acc[4:]); lo=(n-1)//1000*1000+1
    rng=(f'GCST{lo:08d}-GCST{lo+999:08d}' if len(acc)>10 else f'GCST{lo:06d}-GCST{lo+999:06d}')
    base=f'{F}/{rng}/{acc}/harmonised/'
    html=subprocess.run(['curl','-sS','-m','60',base],capture_output=True,text=True).stdout
    import re
    fs=re.findall(r'href="([^"]+\.h\.tsv\.gz)"',html)
    return base+fs[0] if fs else None
def header(url):
    raw=subprocess.run(['curl','-sS','-m','60','-r','0-65535',url],capture_output=True).stdout
    d=zlib.decompressobj(31); txt=d.decompress(raw).decode(errors='ignore')
    return txt.split('\n')[0].split('\t')
def region(url,chrom,start,end):
    h=header(url); tb=pysam.TabixFile(url)
    rows=[r.split('\t') for r in tb.fetch(str(chrom),int(start),int(end))]
    return pd.DataFrame(rows,columns=h[:len(rows[0])] if rows else h)
def std(df):
    low={c.lower():c for c in df.columns}
    def g(*ns):
        for n in ns:
            if n in low: return df[low[n]]
        return pd.Series([None]*len(df))
    o=pd.DataFrame({'rs':g('hm_rsid','rsid','variant_id'),'chr':g('hm_chrom','chromosome'),'pos':g('hm_pos','base_pair_location'),'ea':g('hm_effect_allele','effect_allele'),'oa':g('hm_other_allele','other_allele'),
                    'eaf':g('hm_effect_allele_frequency','effect_allele_frequency'),'beta':g('hm_beta','beta'),'se':g('standard_error'),'p':g('p_value')})
    for c in ['pos','eaf','beta','se','p']: o[c]=pd.to_numeric(o[c],errors='coerce')
    o['ea']=o.ea.str.upper(); o['oa']=o.oa.str.upper()
    return o
