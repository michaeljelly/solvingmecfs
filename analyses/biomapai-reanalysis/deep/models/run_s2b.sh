cd /mnt/project-files/biomap/deep/models
while pgrep -f "run_s2.sh" > /dev/null; do sleep 30; done
python3 s2_classify.py short metabqc enet,lgbm
python3 s2_classify.py long metabqc enet,lgbm
python3 s2_fusion.py all,short,long
