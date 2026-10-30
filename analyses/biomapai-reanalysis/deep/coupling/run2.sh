export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /mnt/project-files/biomap/deep/coupling
while [ ! -f run1.done ]; do sleep 20; done
python3 01_predictability.py metab species Short 30 30 > log_pred_metab_species_Short.txt 2>&1
echo DONE > run2.done
