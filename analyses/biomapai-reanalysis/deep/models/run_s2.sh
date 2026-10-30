#!/bin/bash
cd /mnt/project-files/biomap/deep/models
for T in all short long; do
  python3 s2_classify.py $T demo enet
  python3 s2_classify.py $T immune,specie,metabolome,quest enet,rf,lgbm
  python3 s2_classify.py $T kegg enet,rf,lgbm
done
python3 s2_fusion.py all,short,long
