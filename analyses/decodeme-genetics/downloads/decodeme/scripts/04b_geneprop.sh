#!/bin/bash
# gene-property (tissue/cell-type) analyses only; usage 04b_geneprop.sh <gwas>
G=$1; W=/tmp/work_decodeme; M=$W/magma/magma; O=$W/magma_out
for c in gtex54 hpa_sc monaco_immune brain_sn; do
  $M --gene-results $O/$G.genes.raw --gene-covar $W/ref/covar_$c.txt --model condition-hide=Average direction=greater --out $O/$G.gp.$c > /dev/null
done
echo DONE-GP $G
