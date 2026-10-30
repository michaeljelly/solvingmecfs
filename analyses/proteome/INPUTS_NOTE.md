# Inputs note: proteome / EV stream (deep2/proteome)

Local working copies are in /tmp/work_proteome/raw (645 MB) and were not copied here. The harmonised matrices (log2, samples x features) are in /tmp/work_proteome/data/<ds>_X.parquet with _meta.csv and _feat.csv files. Scripts to rebuild everything are in scripts/, starting with 01_load.py.

| dataset | source | local file | size | format | rows x cols | groups |
|---|---|---|---|---|---|---|
| Hoel/Fluge SomaScan 7k serum | Europe PMC PMC13006441 supplement mmc2.xlsx (sheet 2_log2Matrix, 3_Aptamer_annotations); PAD000026 | raw/supp/PMC13006441/ | 40 MB | xlsx | 79 x 7326 aptamers | 50 ME, 29 HC (24 HC with BMI) |
| Germain 2025 SomaScan 7k plasma, 2-day CPET | PMC12757482 mmc1.xlsx sheet Somalogic_dataset | raw/supp/PMC12757482/ | 78 MB | xlsx | 500 x 7288 | 4 timepoints; D1PRE 78 ME, 51 HC |
| Walitt NIH SomaScan 1.3k CSF | GEO GSE251790 (adat) | downloads/nih_walitt/raw/geo/GSE251790_adat.txt.gz (existing) | - | adat | 42 x 1299 | 21/21 |
| Walitt NIH SomaScan 1.3k serum | GEO GSE254030 (adat) | same folder | - | adat | 42 x 1299 | 21/21 |
| Rydland 2026 plasma EV (PXD072203) | PMC13312568 mmc1.xlsx (the PRIDE 7z has no txt folder) | raw/supp/PMC13312568/ | 0.7 MB | xlsx | 98 x 158 genes | 48 ME, 50 HC |
| Giloteaux 2024 EV, TMT, 0h/15min/24h after CPET | PMC10764978 JEV2-13-12403-s002.xlsx | raw/supp/PMC10764978/ | 1.7 MB | xlsx | 105 x 301 | 54 ME / 51 HC samples, 3 times |
| Bragee 2026 CSF (PXD076216) | PMC13194927 MOESM1 sheet VSN_ProteinQuant, plus PRIDE proteinGroups.txt and sdrf | raw/supp/PMC13194927, raw/pride/ | 0.4 MB + 20 MB | xlsx, tsv | 31 x 902 (1586 protein groups in MaxQuant) | patients only: 10 POTS, 21 no POTS |
| IgG-bound serum proteins (PXD065439) | PRIDE txt/proteinGroups.txt | raw/pride/pxd065439/txt | 11 MB | MaxQuant tsv | 106 x 689 | 67 cases (ME + post-COVID ME), 39 HC |
| UK Biobank Olink ME/CFS summary stats (Beentjes 2025) | paper supplement | raw/supp/ | - | xlsx | ~2900 proteins | summary statistics only (z) |
| Long COVID UK Olink Explore | Zenodo 14772494 plate1-3_olink.csv, meta_UKsamples.csv | raw/zen/ | 89 MB | long csv | CA 51 vs CO 51 | long COVID vs convalescent |
| Long COVID IRIS 3 months | Zenodo 13363117 IRIS_LongCOVID_olink_data.csv | raw/zen/ | 1 MB | csv | 36 vs 9 | |
| Long COVID Bansal | Zenodo 15237088 baseline57.xlsx | raw/zen/ | 1.2 MB | xlsx | 12 LC vs 44 | |
| Jonsjo 2020 Olink inflammation | Kopia_av_MEFCS_Inflammation_Shared_dataset.xlsx | raw/zen/ | 15 KB | xlsx | patients only | no controls, not used |
| Baraniuk/Schutzer 2023 CSF TMT | raw/baraniuk2023_csf.xlsx | | | xlsx | ME/CFS +/- FM | no controls, not used |
| Nunes 2024, Milivojevic 2020, Dara 2025, Glass 2025 | supplements | raw/*.xls(x) | small | | significant-protein lists only | not per-sample |
| Blood RNA: GSE227375, GSE214283, GSE98139, GSE293840 | GEO | reused from other streams / downloaded | | | 23/25, 30/28, 28/19, 93/75 | |
| KEGG REST lists; HPA tissue/cell enriched lists | rest.kegg.jp; proteinatlas.org | raw/kegg, raw/hpa | small | tsv | | |

## Sample rows (copied data, first 3 samples x 4 features, log2)
```
[copied data] hoel_serum        seq.10000.28 seq.10001.7 seq.10003.15 seq.10006.25
hoel0  10.54  9.04  8.41  9.78
hoel1  10.39  8.75  8.82  9.83
[copied data] germain_plasma    10000-28 10001-7 10003-15 10006-25
COR-1002_D1PRE   9.37 9.36 7.43 9.39
COR-1002_D1POST  9.34 8.88 7.15 9.32
[copied data] walitt_csf        10336-3_3 10337-83_3 10339-48_3 10342-55_3
Fatigue-126  4.19 4.04 14.75 3.64
[copied data] bragee_csf        SSC5D IGLON5 SHISA7 VSTM2B
P01_severe_POTS    NaN   20.71 20.61 24.78
P02_severe_NoPOTS  19.90 22.34 NaN   22.94
[copied data] igg_serum         pg0 pg1 pg2 pg3
1336_PT002_Healthy-Control_01  24.95 NaN NaN NaN
```

## Surprises
- Hoel serum: 982 case-control FDR hits, but they are dominated by intracellular proteins leaking from cells (histones, VIM, ALDOA). A 100-aptamer leak score has d = -1.16 in Hoel but -0.02 in Germain, +0.13 in Walitt serum and 0.0 in Walitt CSF, so the artefact belongs to the Hoel samples and not to the disease. The Sex code is inverted relative to its obvious reading (Sex==0 is male).
- Walitt CSF adat has no 'Group' column; group and sex had to be found from the values.
- Excel turned gene symbols into dates in several supplements (e.g. SEPT/MARCH genes).
- PXD072203 (7 GB combined_final.7z) has no MaxQuant txt folder; PXD073644 and PXD068504 are raw files only; figshare 30195544 is embargoed; PXD016622 (Milivojevic) is Spectronaut .sne only (3.2 GB), not reprocessed.
- Long COVID UK IDs differ between files (CA009 vs CA9).
- GSE227375 row index mixes symbols and Entrez IDs.
