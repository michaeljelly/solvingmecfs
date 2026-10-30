# INPUTS_NOTE: DecodeME GWAS and supporting data (downloaded 2026-09-24)

## 1. DecodeME GWAS summary statistics (primary dataset)
- Source: OSF project https://osf.io/rgqs3/ , folder "DecodeME Summary Statistics" (linked from https://institute-genetics-cancer.ed.ac.uk/decodeme). Open download, no application needed.
- Path: `sumstats/` (raw REGENIE v4.1 output, space-separated, gz). Build GRCh38. Only MAF >= 0.01 filter applied; README says keep only variants in `gwas_qced.var.gz` (MAF >= 0.01, INFO >= 0.9, DENTIST-like test).

| file | size | rows (raw) | rows after QC list | cases / controls |
|---|---|---|---|---|
| gwas_1.regenie.gz (primary, GWAS-1) | 318 MB | 8,941,712 | 8,902,782 | 15,579 / 259,909 |
| gwas_2.regenie.gz (all cases vs UKB subset) | 317 MB | - | 8,895,648 | 15,579 / 155,790 |
| gwas_1_female.regenie.gz | 450 MB | - | 8,902,782 | 12,833 / 218,949 |
| gwas_1_male.regenie.gz | 305 MB | - | 8,880,709 | 2,746 / 40,960 |
| gwas_1_infectious_onset.regenie.gz | 318 MB | - | 8,902,782 | 9,738 / 259,909 |
| gwas_1_non_infectious_onset.regenie.gz | 316 MB | - | 8,900,939 | 5,841 / 259,909 |
| gwas_qced.var.gz (QC-passed IDs) | 39 MB | 8,902,782 | | |
| imputed.info.gz (INFO per variant) | 145 MB | | | |
| README_shared_sum_stats.txt, DecodeME-GWAS-Analysis-Plan-v1.pdf | small | | | |

17 columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS N N_CASES N_CONTROLS TEST BETA SE CHISQ LOG10P EXTRA. ALLELE1 is the effect allele; BETA is log-OR (Firth where needed; EXTRA all NA after QC).

```text
# copied data: gwas_1.regenie.gz, first 2 rows
CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS N N_CASES N_CONTROLS TEST BETA SE CHISQ LOG10P EXTRA
5 90000126 5:90000126:G:A G A 0.0154907 0.0161114 0.0164192 275488 15579 259909 ADD 0.0388003 0.0514265 0.569241 0.346247 NA
5 90002604 5:90002604:T:TG T TG 0.16074 0.164805 0.170375 275488 15579 259909 ADD -0.000160003 0.0172172 8.63641e-05 0.00323221 NA
```
Surprises:
- Files are not sorted by chromosome (start at chr5).
- README warns A1FREQ_CONTROLS is miscalculated by REGENIE v4.1; do not use it.
- Preprint Table 3 caption lists 6:26329103:A:T and 20:47529913:C:T as "GRCh38"; they are really GRCh37 coordinates (in GRCh38 they are 6:26328875:T:A and 20:48913376:C:T, the GWAS-Infection leads).
- The 12q24 SUDS3/TAOK3 lead is a poly-T indel (ID 12:118202773:CTTTTTTTTTTTTT:C).
- gwas_2 and gwas_1 share all cases, so they are not independent.

## 2. Paper files
- Preprint PDF (45 pages): `paper/DecodeME_Preprint_2025-08.pdf` from https://institute-genetics-cancer.ed.ac.uk/sites/default/files/2026-05/2025-08-03%20DecodeME%20Preprint.pdf (medRxiv 10.1101/2025.08.06.25333109). The supplementary xlsx files on research.ed.ac.uk (Pure) are behind a Cloudflare challenge and could not be fetched with curl: https://www.research.ed.ac.uk/files/533352488/Supplemental_Tables_S1-6_8-9_11.xlsx, .../533352486/Supplemental_Table_S7.xlsx, .../533352482/Supplement.pdf, .../533352484/Candidate_Genes.pdf. You need a browser to get them.

## 3. Reference data used for gene and enrichment analysis (in /tmp/work_decodeme, not copied: large or re-downloadable)
- MAGMA v1.10 static Linux binary, NCBI37.3/NCBI38 gene locations, 1000 Genomes EUR reference (g1000_eur, 503 people, 22.7M SNPs): https://cncr.nl/research/magma/ (SURFdrive links vu.data.surfsara.nl/index.php/s/lxDgt2dNdNr6DYt, Pj2orwuF2JYyKxq, VZNByNwpD8qqINe).
- MSigDB 2024.1.Hs Entrez GMTs (hallmark, Reactome, KEGG MEDICUS, GO BP/CC/MF, C8 cell type, C7 ImmuneSigDB): https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/
- GTEx v8 gene median TPM (54 tissues): https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
- Human Protein Atlas: rna_single_cell_type.tsv (154 cell types), rna_immune_cell_monaco.tsv (29 immune cell types incl. Th17 and Th1/Th17 memory CD4), rna_single_nuclei_cluster_type.tsv (34 human brain clusters): https://www.proteinatlas.org/download/tsv/
- NCBI Homo_sapiens.gene_info.gz (Ensembl to Entrez): https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/
- UCSC hg38ToHg19.over.chain.gz (for liftover to match the MAGMA 1000G reference).
- eQTLGen cis-eQTL full results (whole blood, n up to 31,684; 4.6 GB): https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTLs_full_20180905.txt.gz . Region extract saved as `results/eqtlgen_region_extract.tsv.gz`.
- GTEx v8 single-tissue eQTL/sQTL via API https://gtexportal.org/api/v2/association/singleTissueEqtl
- LDSC European LD scores (eur_w_ld_chr, HapMap3): https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/ (the Broad bucket now returns 404).

## 4. Other GWAS used for cross-trait checks (in /tmp/work_decodeme; not copied)
- Long COVID, HGI (Lammi et al. 2025 Nat Genet, PMID 40399555), GWAS Catalog GCST90454541 harmonised, GRCh38, ~6,000 cases and ~1.09M controls: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90454001-GCST90455000/GCST90454541/harmonised/GCST90454541.h.tsv.gz (318 MB)
- Multisite chronic pain (Johnston 2019, UKB, n=387,649), GCST008512: .../GCST008001-GCST009000/GCST008512/chronic_pain-bgen.stats.gz (BOLT-LMM, rsIDs, GRCh37)
- Major depressive disorder GCST90726344 (PMID 41233554; 151,582 cases / 362,873 controls, European), harmonised GRCh38.
- Not usable: Howard 2019 depression GCST007342 (the catalog only has the top 10,000 SNPs because of 23andMe terms).

## 5. Also found (not downloaded)
- Zenodo 20204356 (Maccallini 2026): DecodeME + Million Veteran Program ME/CFS meta-analysis sumstats (19,470 cases / 699,111 controls; GRCh37 and GRCh38, 312 MB each) plus the preprint. This is the largest ME/CFS GWAS available and would be the next thing to add.
- Zenodo 18198370 (Wielscher 2026): genetic correlation/MR of DecodeME with 22 traits (record, no files listed).
- Zenodo 22885464 (Huang et al. 2026): transcript-QTL colocalisation package for DecodeME loci (63 MB zip).
