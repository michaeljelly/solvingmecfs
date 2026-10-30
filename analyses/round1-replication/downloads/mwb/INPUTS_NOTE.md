# Inputs note: ME/CFS metabolomics downloads (folder `downloads/mwb/raw/`)

Everything below was downloaded on 24 Sep 2026. Total is about 16 MB. The loaders and parsing code are `mwparse.py` and `01_build.py`.

| Dataset | Path (under raw/) | Size | Format | Rows x cols | Notes |
|---|---|---|---|---|---|
| Naviaux 2016 PNAS, MW ST000450 (plasma, 45 CFS / 39 ctrl) | AN000705.mwtab.txt, AN000706.mwtab.txt, ST000450_*.json | 250 KB + 232 KB | mwTab text | 84 samples x 226 (pos) and 223 (neg) metabolites | Sex is the only per-sample covariate. KEGG and PubChem IDs are in METABOLITES block. No IPA and no 3-hydroxybutyrate. |
| Naviaux 2017 re-run, MW ST000617 | AN000946/947.mwtab.txt | 97 + 102 KB | mwTab | 38 x 189/194 | The same male participants as ST000450, re-measured on another instrument. **Not an independent cohort**, so it is excluded from the meta-analysis. HMDB IDs given. |
| Nagy-Szakal 2018 Sci Rep, MW ST000800/801/802/813 (plasma, 50/50) | ST000800_AN001271/2_Results.txt (CSH lipidomics), ST000801_AN001273_Results.txt (HILIC), AN001274.mwtab.txt (oxylipins), ST000813_AN001286_Results.txt (GC-TOF) | 1.0 MB, 0.36 MB, 1.6 MB, 0.5 MB | separate tab-delimited results files | about 100 samples x 1,900 / 660 / 2,700 / 264 features (mostly unannotated) | The CSH files label samples 1..100. They map to ESP IDs through the `_NNN` suffix in the mwTab SUBJECT_SAMPLE_FACTORS. Covariates: sex, age (<=50 / >50), BMI, IBS. 7 "Biorec" QC samples. No IPA, kynurenine or bile acids (except GUDCA). |
| CSF, MW ST000910-913 (32 ME/CFS, 40 MS, 19 no-neuro-disease) | AN001480-1483.mwtab.txt, ST000913_AN001483_Results.txt | 0.5 MB etc | mwTab + results file | 91 x 309/321/102/116 | Sample names differ between assays, so they are joined on subject (DPCSF IDs). Age and sex given. |
| NIH Walitt 2024 Nat Commun, CSF, MW ST003178 (17 PI-ME/CFS / 21 HV) | AN005217-5220.mwtab.txt | 111 KB etc | mwTab (Metabolon) | 38 x 195/36/131/83 | Age, sex, BMI, race. HMDB, KEGG and PubChem IDs given. **No plasma metabolomics was deposited.** |
| NCNP Japan serum, MW ST004940 (28 ME/CFS / 10 HC) | AN008373.mwtab.txt | 24 KB | mwTab | 38 x 28 targeted tryptophan/indole metabolites (nmol/L) | "N.D." = not detected. No age or sex per sample. Released July 2026 and probably unpublished. |
| NCNP Japan stool, MW ST004941 (15/15) | AN008374.mwtab.txt | 24 KB | mwTab | 30 x 28 (pmol/g dry faeces) | Controls were age- and sex-matched, but no per-sample covariates are given. |
| UConn stool indoles, MW ST003346 (10 ME/CFS / 10 HC) | AN005481.mwtab.txt | 60 KB | mwTab | 180 rows = 20 people x several aliquots x triplicate; 5 indoles (uM) | Averaged per person. Sex is given. |
| Long COVID vs recovered, MW ST003103 (40/40, Spain) | AN005077-5081.mwtab.txt | 0.5 MB etc | mwTab | 92 x 422/107/76/125/90 | Downloaded but **not analysed**. It is not an ME/CFS study, and names are m/z-RT style for most features. |
| Armstrong 2015, MetaboLights MTBLS161 (serum + urine NMR, 34 CFS / 25 ctrl) | MTBLS161/*.txt, *_maf.tsv | 28 KB MAF | ISA-Tab | 41 metabolites x 118 samples (serum and urine) | Only tryptophan overlaps the leads. ChEBI IDs are given. |
| Hoel 2021 JCI Insight (83 ME/CFS / 35 HC, Metabolon + Lipidyzer) | Hoel2021/sdd1-4.xlsx, sd1-3.pdf | 0.28 / 0.54 / 0.07 / 0.24 MB | Excel | 610 metabolites x fold change and Welch p; 18 lipid classes | **Summary statistics only**, no per-sample data. Effects are rebuilt from p and fold change. |
| Germain 2018 Metabolites (Hanson; 32 ME/CFS / 19 HC) | Germain2018_Metabolites/Supplementary_File_1.xlsx | 471 KB | Excel, "OrigScale" | 832 metabolites x 51 samples (C1-19, P1-32) | Per-sample Metabolon data, **with no covariates**. HMDB, KEGG and PubChem IDs given. Obtained via Europe PMC because MDPI returns 403. |
| Germain 2020 Metabolites (Hanson/Barupal; 26/26) | Germain2020_Metabolites/Supplementary File 1.xlsx (HD4), Supplementary File 2.xlsx (CLP lipids) | 463 + 540 KB | Excel | 768 x 52; 1,022 lipids x 52 | No covariates. It may share participants with Germain 2018 or 2022, since all are from the Hanson lab (unverified). |
| Baraniuk 2025 Sci Rep CSF (60 ME/CFS / 32 sedentary controls) | Baraniuk2025_CSF/MOESM2 (per sample) + MOESM1/3/4 | 869 KB | Excel | 92 x 948 (metabolites, lipids, questionnaires) | Values are already standardised. There are exercise and non-exercise groups (different people). Symptom and pain questionnaires are included. |
| Walitt 2024 Nat Commun supplements | Walitt2024_NatCommun/*.xlsx | small | Excel | SD20: 13 plasma lipid classes x 33 people; SD21: 32 stool NMR metabolites x 32; Fig 6J/K: CSF 445 metabolites x 39 | These are figure source data. |
| Guo 2023 Cell Host Microbe tables | Guo2023_CHM/*.xlsx | small | Excel | model tables only | No per-sample plasma data. Not used. |

Sample rows (copied data, ST004940 serum, nmol/L):
```
Samples                 HC-01  HC-02  HC-03  HC-04  HC-05
Indole                  75.6   96.1   62.0   76.7   72.7
Indole-3-acetaldehyde   N.D.   N.D.   1.7    6.9    1.5
```
Sample rows (copied data, Germain 2018 OrigScale header):
```
PATHWAY SORTORDER  BIOCHEMICAL          SUPER PATHWAY  SUB PATHWAY                          ... HMDB       C1
398                1-methylguanidine    Amino Acid     Guanidino and Acetamido Metabolism   ... HMDB01522  63132
```

Surprises:
- ST000617 is a re-measurement of Naviaux 2016 men, not a validation cohort.
- NIH ST003178 holds CSF only.
- The Japanese ST004940/4941 serum and faecal indole panels (2026) are new, and so far are only on Metabolomics Workbench.
- The Nagy-Szakal CSH results files use bare column numbers 1-100.
- MetaboLights has only one ME/CFS study (MTBLS161).
