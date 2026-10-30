# MHC extracts (remote tabix) for negative-control serology traits (same Butler-Laporte 2020 UKB platform) and other traits
import sys; sys.path.insert(0,'/mnt/project-files/mecfs-research/deep2/genetics_mr/scripts'); import tabixlib
T={'BKVP1_BL2020':'GCST90006885','JCVP1_BL2020':'GCST90006922','MCVVP1_BL2020':'GCST90006924','HSV1mgG_BL2020':'GCST90006918','VZVgEgI_BL2020':'GCST90006929',
   'HHV6IE1A_BL2020':'GCST90006904','HHV7U14_BL2020':'GCST90006909','HpyloriOMP_BL2020':'GCST90006914','ToxoSag1_BL2020':'GCST90006927','CtrachPGP3_BL2020':'GCST90006889',
   'HHV6IE1B_BL2020':'GCST90006906','HpyloriGroEL_BL2020':'GCST90006913'}
for eid,acc in T.items():
    try:
        u=tabixlib.acc_url(acc); d=tabixlib.std(tabixlib.region(u,6,25000000,35000000))
        d.to_csv(f'/tmp/work_gmr/exp/{eid}.mhc.tsv',sep='\t',index=False); print(eid,acc,len(d),flush=True)
    except Exception as e: print(eid,'fail',e,flush=True)
