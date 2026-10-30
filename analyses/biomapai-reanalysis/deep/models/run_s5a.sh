cd /mnt/project-files/biomap/deep/models
while pgrep -f "s3_external.py" > /dev/null; do sleep 30; done
python3 s5_importance.py all metabqc lgbm
python3 s5_importance.py all metabqc enet
python3 s5_importance.py all specie rf
python3 s5_importance.py all metabolome lgbm
python3 s5_importance.py all immune rf
