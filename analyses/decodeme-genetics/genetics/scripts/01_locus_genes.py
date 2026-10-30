"""List genes within +/-500 kb of each DecodeME GWAS-1 lead variant (GRCh38).
Lead variants and ORs from DecodeME preprint Table (medRxiv 10.1101/2025.08.06.25333109).
Gene coordinates: MAGMA NCBI38.gene.loc (via github comorment/magma)."""
import pandas as pd
loc = pd.read_csv('/mnt/project-files/mecfs-research/genetics/data/NCBI38.gene.loc', sep='\t', header=None,
                  names=['entrez','chr','start','end','strand','gene'])
leads = pd.DataFrame([
 ('1q25.1','1',173846152,'T:C',0.927,0.325,7.59,'RABGAP1L'),
 ('6p22.2','6',26239176,'A:G',1.086,0.261,8.39,'BTN2A2'),
 ('6q16.1','6',97984426,'C:CA',0.934,0.546,7.31,'FBXL4'),
 ('12q24.23','12',118202773,'C:T',1.100,0.139,6.78,'SUDS3'),
 ('13q14.3','13',53194927,'GT:G',1.077,0.287,6.94,'OLFM4'),
 ('15q21.3','15',54866724,'A:G',1.082,0.312,8.12,'CCPG1'),
 ('17q22','17',52183006,'C:T',1.084,0.330,8.67,'CA10'),
 ('20q13.13','20',48914387,'T:TA',1.095,0.634,11.02,'ARFGEF2/CSE1L'),
], columns=['band','chr','pos','alleles','OR','A1freq','mlog10p','decodeme_gene'])
rows=[]
for _,l in leads.iterrows():
    g = loc[(loc.chr.astype(str)==l.chr)&(loc.end>l.pos-500_000)&(loc.start<l.pos+500_000)].copy()
    g['dist_kb'] = g.apply(lambda r: 0 if r.start<=l.pos<=r.end else min(abs(r.start-l.pos),abs(r.end-l.pos))/1000, axis=1)
    g = g[~g.gene.str.startswith(('LOC','LINC','MIR','SNOR'))].sort_values('dist_kb')
    for _,r in g.iterrows():
        rows.append(dict(band=l.band, lead_pos=l.pos, decodeme_gene=l.decodeme_gene, gene=r.gene, dist_kb=round(r.dist_kb,1)))
    print(f"{l.band} (lead {l.chr}:{l.pos}, OR {l.OR}, {l.decodeme_gene}): " + ', '.join(f"{r.gene}({r.dist_kb:.0f})" for _,r in g.head(12).iterrows()))
pd.DataFrame(rows).to_csv('/mnt/project-files/mecfs-research/genetics/locus_genes_500kb.csv', index=False)
leads.to_csv('/mnt/project-files/mecfs-research/genetics/decodeme_lead_variants.csv', index=False)
# per-allele risk: absolute risk shift at 1% prevalence for heterozygote / homozygote of risk allele
import numpy as np
print("\nAbsolute effect of one risk allele at baseline prevalence ~0.9% (DecodeME-like UK estimate range 0.2-1%):")
for _,l in leads.iterrows():
    orr = l.OR if l.OR>1 else 1/l.OR
    p0=0.009; odds=p0/(1-p0)
    p1=odds*orr/(1+odds*orr); p2=odds*orr**2/(1+odds*orr**2)
    print(f"  {l.decodeme_gene:14s} OR/allele {orr:.3f}: 0.90% -> {p1*100:.2f}% (1 copy) -> {p2*100:.2f}% (2 copies)")
