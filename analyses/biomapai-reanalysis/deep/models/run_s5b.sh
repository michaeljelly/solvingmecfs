cd /mnt/project-files/biomap/deep/models
python3 s5_importance.py short metabolome lgbm
python3 s5_importance.py short specie rf
while [ ! -f results/s2/pred_short_metabqc_lgbm.csv ]; do sleep 30; done
python3 s5_importance.py short metabqc lgbm
