export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /mnt/project-files/biomap/deep/coupling
while [ ! -f run2.done ]; do sleep 20; done
python3 01_predictability.py immune kegg all 10 10 > log_pred_immune_kegg_all.txt 2>&1
echo DONE > run3.done
