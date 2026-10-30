# INPUTS_NOTE: GEO / PRIDE / figshare downloads (geo worker, 2026-09-24)

All files are under `/mnt/project-files/mecfs-research/downloads/geo/data/`. Catalog of every ME/CFS omics dataset found: `results/DATASET_CATALOG.csv`. Raw search hits: `data/gds_hits.json`, `data/gds_hits2.json`, `data/otherrepos.txt`.

Not copied here (too big): GSE214283_RAW.tar (1.3 GB, 116 10x `filtered_feature_bc_matrix.h5` files) from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214283/suppl/GSE214283_RAW.tar . The pseudobulk made from it is in `processed/`.

## 1. GSE214283 pseudobulk (made here from the 10x h5 files)
- Path: `processed/GSE214283_pseudobulk_{ALL,CD4T_other,CD8T,NK,B,cMono,ncMono}.parquet` (4-7 MB each), `processed/GSE214283_celltype_counts.csv`, `processed/GSE214283_meta.csv`
- Format: parquet, genes x samples raw UMI sums; 36,591 genes x 116 samples (58 people x D1 pre-CPET / D2 24 h post). 30 ME/CFS, 28 controls.
- Cell types come from simple marker scoring (CD3/CD8A/NKG7/MS4A1/LYZ/FCGR3A), not the authors' labels. Cells kept: >=200 genes, <20% mito, >=500 UMI (about 4,000 cells per sample).
```
[copied data] pseudobulk_ALL
          GSM6603110  GSM6603111  GSM6603112  GSM6603113
A1BG             950         808         967         834
A1BG-AS1          75          96         119          82
```
- Surprise: one sample (GSM6603131, COR-5218-D1) has only 881 cells.

## 2. GSE227375 (Van Booven 2023, PBMC RNA-seq, exercise T0/T1/T2)
- Path: `raw_small/GSE227375_ncbi_counts.tsv.gz` (6.8 MB), `processed/GSE227375_counts_symbol.parquet`, `processed/GSE227375_meta.csv`
- Source: NCBI-generated counts, https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE227375&format=file&file=GSE227375_raw_counts_GRCh38.p13_NCBI.tsv.gz
- 39,375 Entrez genes x 162 samples (GEO lists 187; 25 have no NCBI counts). Authors gave no supplementary matrix.
```
[copied data]
GeneID	GSM7099227	GSM7099228	GSM7099230
100287102	130	548	45
653635	1125	1537	531
```
- Surprise: GEO does not give subject IDs. Sample numbers run in consecutive triplets (e.g. CF0_2653/CF1_2654/CF2_2655). I used that to assign subjects and checked it with genotype-like genes (GSTM1, ERAP2, HLA-DRB4/5, UGT2B17...): within-subject r = 0.84 (median), between-subject r = -0.02. Result: 31 case and 36 control subjects. Column `subject` in the meta file.

## 3. GSE293840 (plasma cell-free RNA, 93 ME/CFS / 75 controls, 2025)
- Path: `raw_small/GSE293840_raw_counts_all.csv.gz` (5 MB), `processed/GSE293840_meta.csv`
- CSV, Ensembl gene (versioned) x 168 samples; columns cfs_cfrna_N mapped to GSM via series-matrix description.
```
[copied data]
,cfs_cfrna_1,cfs_cfrna_2,cfs_cfrna_3
ENSG00000000005.6,0,0,0
```
- Surprise: almost every gene set is about +0.17 SD higher in cases (global composition shift; global score d = +0.36). Set-level results must be read after adjusting for this (the `_globadj` files).

## 4. PAD000026 (PRIDE affinity archive; Hoel/Fluge SomaScan 7k serum, 50 ME/CFS / 29 controls)
- Path: `raw_small/PAD000026_Log2Matrix.tsv` (7 MB), `raw_small/PAD000026.adat` (6 MB, has aptamer to gene annotation), `raw_small/PAD000026_PAD000026_community_annotated.sdrf.tsv`
- Source: https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/01/PAD000026/
- TSV, 79 rows (samples) x 7,335 columns (9 metadata + 7,326 aptamers). Decimal commas.
```
[copied data]
Group	Metabotype	Sex	Age	BMI_Cat	Fasting_baseline	SF36PF	Daily_Step_Count	CycloME	seq.10000.28
ME	Metabotype 1	1	31–35	Normal	0	20	1394		10,54167745
```
- Surprise: 1,108 of 7,326 aptamers differ at q<0.05, led by intracellular proteins (histones, ALDOA, LDHB, VIM) all lower in patients (leak score d = -1.42). This looks like a sample-handling difference between patient and control serum, not disease biology. Any single-protein result needs adjusting for it.

## 5. Smaller sets used
- `raw_small/GSE128078_ncbi_counts.tsv.gz` (3.6 MB; Bouquet 2019 whole blood; 25 people x days 1/2/3/7; 99 samples) plus `GSE128078_FES_isoforms_FPKM.txt.gz`.
- `raw_small/GSE98139_All_normalized_counts.csv.gz` (1.8 MB; semicolon-separated; Ensembl x 47 adolescents; 28 CFS / 19 controls).
- `raw_small/GSE16059_series_matrix.txt.gz` (7.3 MB; GPL570 RMA; 88 arrays = 44 discordant MZ twin pairs). Probe annotation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=data (80 MB, not copied).
- Not analysed but saved: GSE268212 CD8 counts, GSE270045 counts, GSE275334 NanoString xlsx, series matrices for GSE130353, GSE14577, GSE269047.

## 6. figshare 30278434 (Esteban 2025, stool metabolomics + AHR assay, 19 ME/CFS / 23 controls)
- Path: `raw_small/figshare_AHR/` (about 50 small csv/R/qza files). DOI 10.6084/m9.figshare.30278434.v2
```
[copied data] volcano.csv
"","FC","log2(FC)","p.ajusted","-log10(p)"
"indole",3.7455,1.9052,0.0013785,2.8606
"tryptophan",2.4245,1.2777,0.027243,1.5647
```
- Relevant to the IPA lead: fecal indole-3-propionate was NOT different (median 3.15 vs 3.44, Mann-Whitney p = 0.83); fecal indole (p = 0.0004), tryptophan (p = 0.006) and indoleacetate (p = 0.008) were higher in patients. This is stool, not plasma.
