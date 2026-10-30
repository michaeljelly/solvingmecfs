export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /mnt/project-files/biomap/deep/coupling
for a in "metab species all 30 30" "immune species all 30 30" "metab species Short 30 30" "metab kegg all 20 20" "immune kegg all 20 20"; do
  set -- $a
  python3 01_predictability.py $a > log_pred_$1_$2_$3.txt 2>&1
done
echo DONE > run1.done
