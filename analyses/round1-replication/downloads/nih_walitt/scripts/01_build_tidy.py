"""Build tidy person-level tables from Walitt et al. 2024 (Nat Commun 15:907) source data.

Inputs: raw/nature_esm/sourcedata (Source Data zip, MOESM6) and raw/nature_esm/suppdata (Supplementary Data 1-24, MOESM4).
Outputs: results/tidy/*.csv
No participant IDs exist in most figure tables; each table is analysed on its own.
"""
import pandas as pd, numpy as np, os, re
B = '/mnt/project-files/mecfs-research/downloads/nih_walitt'
SRC = f'{B}/raw/nature_esm/sourcedata'
SUP = f'{B}/raw/nature_esm/suppdata'
OUT = f'{B}/results/tidy'
os.makedirs(OUT, exist_ok=True)

def grp(x):
    x = str(x)
    if x.startswith('HV') or x == 'Control': return 'HV'
    if x.startswith('PI') or x.startswith('ME') or x == 'Patient': return 'PI-ME/CFS'
    return np.nan

# 1. simple two-column group tables -> one long table of scalar measures
simple = {
 'Figure 2B.xlsx': 'HRV_SDNNi_ms', 'Figure 2C.xlsx': 'HRV_rMSSD_ms', 'Figure 2D.xlsx': 'HRV_pNN50_pct',
 'Figure 2H.xlsx': 'Baroreflex_slope_ms_per_mmHg', 'Figure 3C.xlsx': 'Grip_max_force_kg', 'Figure 3D.xlsx': 'Grip_time_to_failure_s',
 'Figure 5A.xlsx': 'CPET_peak_power_W', 'Figure 5C.xlsx': 'CPET_peak_HR', 'Figure 5D.xlsx': 'CPET_peak_VO2', 'Figure 5E.xlsx': 'CPET_pct_pred_VO2',
 'Figure 5F.xlsx': 'CPET_HR_reserve', 'Figure 5H.xlsx': 'CPET_VO2_at_AT', 'Figure S8F.xlsx': 'CPET_HR_at_AT',
 'Figure 6A.xlsx': 'CSF_DOPA', 'Figure 6B.xlsx': 'CSF_DOPAC', 'Figure 6C.xlsx': 'CSF_DHPG',
 'Figure 7A.xlsx': 'Blood_B_naive_pct', 'Figure 7B.xlsx': 'Blood_B_switched_memory_pct', 'Figure 7C.xlsx': 'CSF_CD8_PD1_pct',
 'Figure 7D.xlsx': 'Blood_CD8_CD226_pct', 'Figure 7E.xlsx': 'CSF_CD8_CXCR5_pct', 'Figure 7F.xlsx': 'Blood_CD8_naive_pct',
 'Figure S1A.xlsx': 'IENFD_leg', 'Figure S1B.xlsx': 'IENFD_thigh', 'Figure S1C.xlsx': 'SGNFD_leg', 'Figure S1D.xlsx': 'SGNFD_thigh',
 'Figure S1E.xlsx': 'Plasma_Tau', 'Figure S1F.xlsx': 'Plasma_NfL', 'Figure S1G.xlsx': 'Plasma_GFAP', 'Figure S1H.xlsx': 'CSF_Tau',
 'Figure S1I.xlsx': 'CSF_NfL', 'Figure S1J.xlsx': 'CSF_GFAP', 'Figure S1K.xlsx': 'CSF_UCHL1',
 'Figure S2A.xlsx': 'DEXA_lean_mass', 'Figure S2B.xlsx': 'DEXA_fat_mass', 'Figure S2D.xlsx': 'DEXA_body_fat_pct', 'Figure S2E.xlsx': 'DEXA_visceral_fat',
 'Figure S3A.xlsx': 'PBMC_Seahorse_basal', 'Figure S3B.xlsx': 'PBMC_Seahorse_max_resp', 'Figure S3D.xlsx': 'Muscle_type2_type1_ratio',
 'Figure S20A.xlsx': 'Stool_observed_features', 'Figure S20B.xlsx': 'Stool_inverse_simpson',
}
rows = []
for f, name in simple.items():
    d = pd.read_excel(f'{SRC}/{f}')
    d = d.iloc[:, :2].dropna()
    for g, v in zip(d.iloc[:, 0], d.iloc[:, 1]):
        try: v = float(v)
        except: continue
        rows.append(dict(measure=name, group=grp(g), value=v, source=f))
scal = pd.DataFrame(rows).dropna(subset=['group'])
scal.to_csv(f'{OUT}/scalar_measures_long.csv', index=False)

# 2. paired tables (catechols vs effort/grip/cognition). Keep as-is.
paired = ['Figure 3E.xlsx', 'Figure 6D.xlsx', 'Figure 6E.xlsx', 'Figure 6F.xlsx', 'Figure 6G.xlsx', 'Figure S11B.xlsx', 'Figure S11C.xlsx'] + \
         [f'Figure S12{c}.xlsx' for c in 'ABCDEFGHIJKLMNO'] + ['Figure S6B.xlsx', 'Figure S8B.xlsx']
prow = []
for f in paired:
    d = pd.read_excel(f'{SRC}/{f}').iloc[:, :3].dropna()
    for _, r in d.iterrows():
        prow.append(dict(source=f, group=grp(r.iloc[0]), x_name=d.columns[1].strip(), x=float(r.iloc[1]), y_name=d.columns[2].strip(), y=float(r.iloc[2])))
pd.DataFrame(prow).to_csv(f'{OUT}/paired_measures_long.csv', index=False)

# 3. tilt plasma catecholamines (Supp Fig 4), nmol/L presumably; row i in NE file = row i in EPI file (same order, file S4C has both)
t = pd.read_csv(f'{SRC}/Figure S4C.csv')
t.insert(0, 'row', range(len(t)))
t['group'] = t['Group'].map(grp)
t.to_csv(f'{OUT}/tilt_plasma_catecholamines_wide.csv', index=False)

# 4. CSF metabolomics, batch-normalised imputed Metabolon data (Supp Fig 13A = all 38)
m = pd.read_excel(f'{SRC}/Figure S13A.xlsx')
m = m[m['ID'].notna()].copy()
m.insert(0, 'group', m['ID'].map(grp)); m = m.drop(columns='ID')
m.rename(columns={'Birth Sex': 'female'}, inplace=True)  # 1=female: IDs 301,312 are female in GEO GSE251790; counts 7 M PI, 10 M HV match paper
# attach NINR sample ids and male participant numbers from Fig 6J (same row order, one blank row)
j = pd.read_excel(f'{SRC}/Figure 6J.xlsx')
j = j[j['CLIENT_SAMPLE_ID'].notna()].reset_index(drop=True)
assert np.allclose(j.iloc[:, 3].values.astype(float), m.iloc[:, 2].values.astype(float))
m.insert(1, 'sample_id', j['CLIENT_SAMPLE_ID'].values)
m.insert(2, 'participant_no_if_female', [x if str(x)[0].isdigit() else '' for x in j['ID'].astype(str)])
m.to_csv(f'{OUT}/csf_metabolomics_wide.csv', index=False)
# annotation (pathways) from Fig 6I (318 metabolites) and SD14
ann = pd.read_excel(f'{SRC}/Figure 6I.xlsx')[['CHEMICAL_NAME', 'SUPER_PATHWAY', 'SUB_PATHWAY', 'HMDB', 'KEGG', 'pval_combined_cfs', 'fc']]
ann.to_csv(f'{OUT}/csf_metabolomics_paper_stats_fig6I.csv', index=False)

# 5. Flow cytometry (Supplementary Data 15): rows = populations, columns = people (21 HV, 17 PI; no IDs)
fc = pd.read_excel(f'{SUP}/Supplementary Data 15.xlsx', header=None)
out = []; section = None; hdr = None
for i, r in fc.iterrows():
    lab = r[0]
    if isinstance(lab, str) and lab.startswith('Supplementary Data 15') and ':' in lab:
        section = 'blood' if 'blood' in lab else 'csf'; sub = lab.split(':')[0][-3:]; continue
    if r[1] in ('HV',):
        hdr = [grp(x) for x in r.values]; continue
    if isinstance(lab, str) and hdr is not None:
        for k in range(1, len(r)):
            if hdr[k] in ('HV', 'PI-ME/CFS'):
                try: v = float(r[k])
                except: v = np.nan
                out.append(dict(tissue=section, table=sub, population=lab.strip(), person_col=k, group=hdr[k], value=v))
fl = pd.DataFrame(out).dropna(subset=['value'])
fl.to_csv(f'{OUT}/flow_cytometry_long.csv', index=False)

# 6. Plasma lipidomics (Supp Fig 19A, per person, lipid species) and class totals (SD20)
l = pd.read_excel(f'{SRC}/Figure S19A.xlsx')
l = l[l['Group'].notna()].copy(); l['group'] = l['Group'].map(grp)
l.to_csv(f'{OUT}/plasma_lipidomics_species_wide.csv', index=False)
def person_matrix(file, header_row):
    d = pd.read_excel(file, header=None)
    h = d.iloc[header_row]
    body = d.iloc[header_row + 1:]
    recs = []
    for _, r in body.iterrows():
        if not isinstance(r[0], str): continue
        for k in range(1, len(r)):
            g = grp(h[k])
            if g in ('HV', 'PI-ME/CFS'):
                try: recs.append(dict(feature=r[0].strip(), person_col=k, group=g, value=float(r[k])))
                except: pass
    return pd.DataFrame(recs)
person_matrix(f'{SUP}/Supplementary Data 20.xlsx', 2).to_csv(f'{OUT}/plasma_lipid_classes_long.csv', index=False)
person_matrix(f'{SUP}/Supplementary Data 21.xlsx', 3).to_csv(f'{OUT}/stool_nmr_metabolites_long.csv', index=False)

# 7. HRV 5-min segment tables (collapsed) from the GitHub RData were not re-used; person-level HRV from Fig 2B-D above.
# 8. EEfRT trial data with mapIDs (GitHub repo)
e = pd.read_excel('/mnt/project-files/mecfs-research/datasets/repos/docwalitt_National-Institutes-of-Health-Myalgic-Encephalomyelitis-Chronic-Fatigue-Syndrome-Code-Repository/EEfRT/EEfRT_Full_Data_SAS_mapID.xlsx')
e.to_csv(f'{OUT}/eefrt_trials.csv', index=False)
print(scal.groupby(['measure', 'group']).size().unstack().to_string())
print('tilt', t.shape, 'csf metab', m.shape, 'flow', fl.shape, 'lipid species', l.shape)
