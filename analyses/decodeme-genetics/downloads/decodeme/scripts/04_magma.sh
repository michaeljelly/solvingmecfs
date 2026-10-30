#!/bin/bash
# MAGMA v1.10 pipeline for a DecodeME GWAS. Usage: 04_magma.sh <gwas_name>
set -e
G=$1; W=/tmp/work_decodeme; M=$W/magma/magma; O=$W/magma_out; mkdir -p $O
python3 - <<P
import pandas as pd, numpy as np
W='$W/'
m=pd.read_csv(W+'pq/id_rs_map.tsv.gz',sep='\t')
import os
if not os.path.exists(W+'magma_out/snploc.txt'):
    m[['rs','chr','pos19']].to_csv(W+'magma_out/snploc.txt',sep='\t',index=False,header=False)
d=pd.read_parquet(W+'pq/$G.parquet',columns=['ID','LOG10P','N'])
d=d.merge(m[['ID','rs']],on='ID')
d['P']=10**(-d.LOG10P)
d[['rs','P','N']].rename(columns={'rs':'SNP'}).to_csv(W+'magma_out/$G.pval.txt',sep='\t',index=False)
print('$G snps for magma',len(d))
P
[ -f $O/genes_w0.genes.annot ] || $M --annotate --snp-loc $O/snploc.txt --gene-loc $W/magma/NCBI37.3.gene.loc --out $O/genes_w0
$M --bfile $W/magma/g1000_eur synonyms=$W/magma/g1000_eur.synonyms --pval $O/$G.pval.txt ncol=N --gene-annot $O/genes_w0.genes.annot --out $O/$G
for s in h.all c2.cp.reactome c2.cp.kegg_medicus c5.go.bp c5.go.cc c5.go.mf c8.all c7.immunesigdb project_sets; do
  f=$W/gs/$s.sets; [ $s = project_sets ] && f=$W/gs/project_sets.entrez.txt
  $M --gene-results $O/$G.genes.raw --set-annot $f --out $O/$G.gs.$s > /dev/null
done
for c in gtex54 hpa_sc monaco_immune brain_sn; do
  $M --gene-results $O/$G.genes.raw --gene-covar $W/ref/covar_$c.txt --model condition-hide=Average direction=greater --out $O/$G.gp.$c > /dev/null
done
echo DONE $G
