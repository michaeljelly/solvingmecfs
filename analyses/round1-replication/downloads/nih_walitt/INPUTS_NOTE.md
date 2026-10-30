# INPUTS_NOTE: NIH intramural PI-ME/CFS deep phenotyping (Walitt et al. 2024, Nat Commun 15:907)

All paths under `/mnt/project-files/mecfs-research/downloads/nih_walitt/`.

## 1. Nature Communications supplementary files (open, no login)
- Source: https://www.nature.com/articles/s41467-024-45107-3 , files `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-45107-3/MediaObjects/41467_2024_45107_MOESM{1..6}_ESM.*`
- Path: `raw/nature_esm/` (156 MB total). MOESM1 supp info PDF (32 MB), MOESM4 zip = Supplementary Data 1-24 (xlsx, unzipped to `raw/nature_esm/suppdata/`), MOESM6 zip = Source Data, 183 files xlsx/csv/txt/nii (unzipped to `raw/nature_esm/sourcedata/`).
- Person-level content, no participant IDs in most tables (one row per person, group label only):
  - Supp Fig 4 (`Figure S4C.csv`): plasma NE and EPI during 40-min head-up tilt, supine + every 4 min, 16 PI / 17 HV (1 HV all missing). NA after the last draw marks tilt stopped early.
  - Fig 6A-C, 6D-G, S11, S12: CSF DOPA, DOPAC, DHPG (21 HV / 16 PI) and CSF NE, DA, cys-DOPA, cys-DA paired with grip time to failure, hard-task choice and cognitive scores.
  - Fig 6H/6J/6K, S13A-E: CSF Metabolon metabolomics, 445 metabolites x 38 people (21 HV, 17 PI), batch-normalised imputed. 6J gives NINR sample IDs and participant numbers for women only.
  - Supplementary Data 15: blood and CSF flow cytometry, ~80 populations x 38 columns (people).
  - Fig S19A-C: plasma lipidomics, 792 lipid species x 33 people; SD20 lipid-class totals; SD21 stool NMR, 32 metabolites x 33.
  - Fig 2 (HRV, baroreflex), 3 (grip), 5 (CPET, n=9/8), 7 (B/T cells), S1 (nerve fibre, NfL, GFAP, tau), S2 (DEXA), S3 (Seahorse, fibre type).
  - Fig 3A/3B: EEfRT trial-level data (1,565 trials, 31 people) with IDs; same as GitHub `EEfRT_Full_Data_SAS_mapID.xlsx`.
- Tidy versions in `results/tidy/` (built by `scripts/01_build_tidy.py`):
  - `tilt_plasma_catecholamines_wide.csv` 33 rows x 25 cols
  - `csf_metabolomics_wide.csv` 38 rows x 489 cols (`female` 1=female; verified against GEO sex of participants 301/312)
  - `flow_cytometry_long.csv` 2,670 rows; `scalar_measures_long.csv` 1,179 rows (42 measures); `paired_measures_long.csv` 714 rows
  - `plasma_lipidomics_species_wide.csv` 33 x 795; `plasma_lipid_classes_long.csv` 429 rows; `stool_nmr_metabolites_long.csv` 1,056 rows; `eefrt_trials.csv` 1,565 x 18

Copied data (first rows):
```
row,Group,NE.SUPINE,NE.4,NE.8,...,EPI.SUPINE,EPI.4,EPI.8,...
0,PI-ME/CFS,0.59576,1.52975,1.57985,...,0.12081,0.25674,0.33509,...
1,PI-ME/CFS,1.08636,1.88577,2.13857,...,0.42259,0.55545,0.78393,...
group,sample_id,participant_no_if_female,female,S-1-pyrroline-5-carboxylate,homocarnosine,...
HV,NINR-00417,,0.0,1.1265,0.9186,...
HV,NINR-00418,103,1.0,3.6635,3.6393,...
tissue,table,population,person_col,group,value
blood,15A,B cell switched memory (%),1,HV,20.5
```

Surprises:
- "Birth Sex" in the metabolomics sheets is coded 1 = female (not male).
- Supplementary Data 14A: the "adjusted p-value" column is sorted independently of the rows (monotone increasing down the sheet); 70 of 142 rows show adjusted p smaller than raw p, which BH cannot produce. Treat that column as scrambled.
- CSF metabolomics: 94% of the 445 metabolites are lower in patients; HV samples were numbered NINR-00417..437 and PI 00438..454 (group blocks), and within patients the overall level falls with sample number (rho -0.68, p=0.002).
- Fig 6J/6K put one patient sample (ID 301) in the HV block with a blank group label; S13A labels it correctly.

## 2. GitHub code repo (already cloned): docwalitt/National-Institutes-of-Health-...-Code-Repository
- `Nature Comms MECFS Source Data.RData` is an RDS (gzip), readable with python `rdata`; 120 tables that duplicate the Source Data zip (fewer tables). EEfRT data + SAS scripts. Microbiome PPM table (`data/Walitt2023_Microbiome_...xlsx`). No tilt HR/BP data.

## 3. GEO (open)
- `raw/geo/GSE251790_adat.txt.gz` (CSF SomaScan 1.3k, 217 KB) and `raw/geo/GSE254030_adat.txt.gz` (serum SomaScan 1.3k, 251 KB), plus series matrices for GSE251790, GSE254030, GSE245661 (muscle RNA-seq). ADAT: 42 study samples each (21 Patient / 21 Control), i.e. 4 more patients than the 17 adjudicated in the paper; no adjudication flag.
- Not downloaded: GSE245661_RAW.tar and GSE251872_RAW.tar (RNA-seq counts; PBMC and muscle count matrices are already in the Source Data as Fig 8/9 csv), SRA PRJNA954397 / SRP467038 (stool shotgun reads).

## 4. mapMECFS (RTI CKAN)
- Anonymous `package_search` works: `raw/mapmecfs/package_search.json` (119 datasets, full metadata), `raw/mapmecfs/mapmecfs_catalog.csv` (one row per dataset), `raw/mapmecfs/nih_resources.tsv` (all 17 "Post-Infectious MECFS at the NIH" datasets and their resource URLs).
- Every file download, `package_list`, `resource_show` and `datastore_search` returns 403 "Authentication is required". Registration is free (CC BY 4.0) but has a Cloudflare Turnstile check, so it needs a person to register and create an API token.
- The login-gated NIH files that matter most: Head-up Tilt Table Testing (datafile-tilt-table-dataset.csv, 2 KB), Tilt Catecholamines (23 KB), Valsalva (3 KB), Lumbar Puncture CSF Catecholamines (4 KB), Clinical Master Labs (16 KB), CSF/Blood flow cytometry databases, Metabolomics CSF (208 KB), Lipidomics, Stool NMR, EEfRT summary, CPET, accelerometry, all questionnaires. These carry participant IDs, so they would let tilt HR/BP be joined to CSF catechols, metabolomics and immune data.
