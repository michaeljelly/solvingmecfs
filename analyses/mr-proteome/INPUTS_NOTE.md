# INPUTS_NOTE: genetics_mr stream (downloaded 2026-09-24)

Working copies live in `/tmp/work_gmr/` (local disk). Compact derived tables were copied to
`/mnt/project-files/mecfs-research/deep2/genetics_mr/data/`. Reused without re-download: DecodeME QC parquet
(`downloads/decodeme/qc_parquet/`), DecodeME+MVP METAL meta (Zenodo 20204356, in `/tmp/work_decodeme/meta/`),
Johnston 2019 chronic pain and GCST90726344 depression sumstats, eQTLGen, 1000G EUR (MAGMA reference), LDSC scores.

## 1. FinnGen R12 G6_POSTVIRFAT (postviral fatigue syndrome, ICD-10 G93.3)
- URL: https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_G6_POSTVIRFAT.gz (780 MB, open)
- Kept: `/tmp/work_gmr/dl/fg_G6_POSTVIRFAT.tsv.gz`, 299 MB, 11,047,960 variants x 9 columns (af 0.005-0.995), GRCh38. 283 cases / 463,029 controls.
```text
# copied data: fg_G6_POSTVIRFAT.tsv.gz first rows
#chrom	pos	ref	alt	rsids	pval	beta	sebeta	af_alt
1	13668	G	A	rs2691328	0.569987	0.555467	0.977814	0.00600458
1	14773	C	T	rs878915777	0.859634	-0.108451	0.613274	0.0142398
```
- Surprise: only 283 cases; no ME/CFS endpoint in the FinnGen+MVP+UKBB meta (meta_analysis/mvp_ukbb manifest has no POSTVIRFAT row).

## 2. FinnGen R12 AB1_EBV (infectious mononucleosis)
- URL: https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_AB1_EBV.gz (811 MB). 3,838 cases / 486,087 controls.
- Kept `/tmp/work_gmr/dl/fg_AB1_EBV.tsv.gz` (304 MB, 11,047,347 variants, same 9 columns). Only 2 independent genome-wide hits after clumping (both MHC).

## 3. Pan-UKB 20002_1482 (self-reported chronic fatigue syndrome), EUR
- URL: https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20002-both_sexes-1482.tsv.bgz (543 MB). EUR 1,882 cases / 418,591 controls (phenotype_manifest.tsv.bgz).
- Kept `/tmp/work_gmr/dl/panukb_20002_1482.tsv.gz`, 350 MB, 14,438,871 variants (low_confidence_EUR == false), GRCh37, columns chr pos ref alt af_cases_EUR af_controls_EUR beta_EUR se_EUR neglog10_pval_EUR.
```text
# copied data
chr	pos	ref	alt	af_cases_EUR	af_controls_EUR	beta_EUR	se_EUR	neglog10_pval_EUR
1	17641	G	A	1.415e-03	8.285e-04	7.906e-01	6.071e-01	7.149e-01
```
- Surprise: not independent of DecodeME. DecodeME controls are UK Biobank participants, so this set shares controls with DecodeME GWAS-1.

## 4. MVP-only ME/CFS estimates (derived, no download)
- Built from the DecodeME-1+MVP METAL file (sample-size scheme). The meta file's FRQ and Z refer to column A2, not A1 (checked against DecodeME at all leads). MVP-only Z = (Z_meta*sqrt(N) - sqrt(58,792)*Z_DME)/sqrt(N-58,792); MVP effective N 15,427 (3,891 cases / 439,202 controls).
- Output: `data/mvp_only_derived.parquet` (8,506,106 variants; rs chr pos38 ea oa eaf z beta se Neff). Validated against the MVP_Z column of the paper's S1 table at 8 SNPs (for example rs4810894 1.237 vs 1.2408; rs12579722 3.128 vs 3.14).

## 5. GWAS Catalog tables
- Studies: https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-download-studies-v1.0.3.1.txt (129 MB).
- Associations: .../gwas-catalog-associations_ontology-annotated-full.zip (74 MB zip, 741 MB tsv, 1,192,604 rows). Used for lead-SNP instruments: Chen 2023 metabolites (PMID 36635386; 1,400 European accessions GCST90199621-GCST90201020, of which 672 non-ratio metabolites had a p<5e-8 lead), Surendran 2022 (PMID 36357675, 913 studies, catalog has top hits only, no full sumstats), Orru 2020 immune cells (PMID 32929287, 731 traits, 600 with a hit), and risk factors (see results/mr_all_results.csv `source`).
- Surprise: catalog lead SNPs give only the risk allele; the other allele was filled from DecodeME. Some catalog entries are unusable (Revez vitamin D "?" alleles, Nordeidet heart rate given as chr:pos-?), so full files were streamed for those.

## 6. Streamed full summary statistics (kept only p<5e-6 rows, all MHC rows chr6:25-35 Mb, and HapMap3 rows)
Extracts in `/tmp/work_gmr/exp/` (`<id>.tsv`, `<id>.mhc.tsv`, `<id>.hm3.tsv.gz`). Source files (harmonised GRCh38 unless noted):
| id | accession / URL | n |
|---|---|---|
| RHR_Zhu2019 | GCST007609 ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz (GRCh37, space-separated) | 458,969 |
| VitD_Revez2020 | GCST90000618 | 496,946 |
| SleepDur_Dashti2019 | GCST007561 zip (columns BETA_SLEEPDURATION etc.) | 446,118 |
| EBVload_Kamitaki2026 | GCST90809803 (EBV DNA load in blood from WGS reads) | 654,938 |
| HHV6Bload_Kamitaki2026 | GCST90809811 | 645,249 |
| EBV/CMV serology (Butler-Laporte 2020, PMID 33204752, UKB) | GCST90006893/-895/-896/-897/-898/-899/-900/-901 | 5,136-8,735 |
| 12 control-pathogen serologies (BK, JC, MCV, HSV-1, VZV, HHV-6 IE1A/B, HHV-7, H. pylori OMP/GroEL, Toxoplasma, C. trachomatis) | GCST90006885, -889, -904, -906, -909, -913, -914, -918, -922, -924, -927, -929 (MHC only, remote tabix) | 1,000-8,700 |
| Chen 2023 metabolites (21 pre-specified) | GCST90199686 IPA, -650 indolelactate, -683 indoleacetate, -636 kynurenine, GCST90200441 tryptophan, -454 kynurenate, -348 quinolinate, GCST90199626 xanthurenate, GCST90200367 3-hydroxybutyrate, -318 acetoacetate, GCST90199621 carnitine, -669 C2, -666 C16, GCST90200387 GCA, -321 cholate, -386 DCA, -342 GCDCA, -333 TCA, GCST90199689 GLCA, GCST90200422 serotonin, -408 lactate | ~5,800-8,300 |
| Reticulocyte count / immature reticulocyte fraction (Vuckovic 2020) | GCST90002405, GCST90002387 (GRCh37-lifted harmonised; no .tbi) | 408,112 |
| Haemoglobin (Chen MH 2020, EUR) | GCST90002310 | 563,946 |
| EBV read positivity (Schmidt 2026) | GCST90809298 (non-harmonised only) | 56,180 / 304,103 |
Example row (EBVload_Kamitaki2026 harmonised):
```text
# copied data
chromosome	base_pair_location	effect_allele	other_allele	beta	standard_error	effect_allele_frequency	p_value	variant_id	rsid
1	701203	T	G	-0.0068	0.0072	0.0091	0.3434	1_701203_T_G	rs1012520212
```
- Surprise: EBI harmonised files with a `.tbi` can be queried remotely with pysam (about 17 s for the whole MHC), which avoids downloading 300-400 MB files. Vuckovic and Revez harmonised files have no index.

## 7. pQTL instruments
- UKB-PPP (Sun et al. 2023 Nature, n=34,557 discovery, Olink 2,923 proteins): Supplementary Table 9, https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx (34 MB). 1,954 proteins with a cis pQTL (p<1.7e-11); strongest cis variant per protein used.
- deCODE (Ferkingstad 2021 Nat Genet, n=35,559, SomaScan 4,907 aptamers): Supplementary Table 2, https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00978-w/MediaObjects/41588_2021_978_MOESM4_ESM.xlsx (44 MB). 1,688 aptamers with an rsID cis pQTL; no SE given, so SE = |beta|/z from the p-value.

## 8. APIs
- Open Targets Platform GraphQL (api.platform.opentargets.org/api/v4/graphql): `knownDrugs` no longer exists; the field is now `drugAndClinicalCandidates`.
- DGIdb GraphQL (dgidb.org/api/graphql), GTEx v2 API (single-tissue eQTL), GWAS Catalog REST v2 (slow, timed out at page size 500; the release TSVs are faster).

## Not obtained
- OpenGWAS (requires a JWT token since 2024); not needed.
- Full UKB-PPP and deCODE pQTL summary statistics (Synapse login / deCODE download form), needed for formal colocalisation of cis-pQTL hits.
- Surendran 2022 full summary statistics (not deposited in the catalog; only top hits), Eldjarn 2023 immunoglobulin GWAS (not in the catalog; Jonsson 2017 top hits used instead).
- HLA allele dosages for DecodeME (not public), needed to test DQA1*05:01 directly.
