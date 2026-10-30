# Mendelian randomisation, MVP/FinnGen replication and cis-pQTL protein screen

Catalogue id: `analysis-mr-proteome`

**What we found (one line):** Genetically higher IPA, indolelactate, kynurenine, ketones, carnitine or bile acids do not change ME/CFS risk; plasma EPHB1 and CLEC3B are new protein leads that replicate in the Million Veteran Program; DecodeME's genome-wide signal replicates in MVP (rg 0.73) though single loci mostly do not yet.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Mendelian randomisation, MVP/FinnGen replication and cis-pQTL protein screen. https://data.solvingmecfs.org/analysis-mr-proteome/
- **Source:** https://data.solvingmecfs.org/analysis-mr-proteome/README.md
- **Sample size:** see README
- **Data type:** our analysis: Mendelian randomisation, pQTL, HLA; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Includes MVP-only statistics we derived from the CC BY Zenodo 20204356 meta-analysis and MHC-region extracts of GWAS Catalog summary statistics (CC0) and DecodeME. Excluded: the FinnGen MHC extract (no explicit licence) and the pQTL instrument table built from the deCODE supplement (not openly licensed).

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `decodeme-gwas`: https://data.solvingmecfs.org/decodeme-gwas/README.md (if mirrored) or see the catalogue
- `decodeme-mvp-meta`: https://data.solvingmecfs.org/decodeme-mvp-meta/README.md (if mirrored) or see the catalogue
- `finngen-r12`: https://data.solvingmecfs.org/finngen-r12/README.md (if mirrored) or see the catalogue
- `gwas-catalog`: https://data.solvingmecfs.org/gwas-catalog/README.md (if mirrored) or see the catalogue
- `panukb-20002-1482`: https://data.solvingmecfs.org/panukb-20002-1482/README.md (if mirrored) or see the catalogue
- `pqtl-ukbppp-decode`: https://data.solvingmecfs.org/pqtl-ukbppp-decode/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-mr-proteome/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-mr-proteome/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `INPUTS_NOTE.md` | 7.4 KB |
| `data/mhc_extracts/BKVP1_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/CMVpp150_BL2020.mhc.tsv` | 4.7 MB |
| `data/mhc_extracts/CMVpp52_BL2020.mhc.tsv` | 4.7 MB |
| `data/mhc_extracts/CMVsero_BL2020.mhc.tsv` | 4.7 MB |
| `data/mhc_extracts/CtrachPGP3_BL2020.mhc.tsv` | 4.3 MB |
| `data/mhc_extracts/EAD_BL2020.mhc.tsv` | 4.7 MB |
| `data/mhc_extracts/EBNA1_BL2020.mhc.tsv` | 4.8 MB |
| `data/mhc_extracts/EBVload_Kamitaki2026.mhc.tsv` | 4.6 MB |
| `data/mhc_extracts/EBVsero_BL2020.mhc.tsv` | 4.6 MB |
| `data/mhc_extracts/HHV6Bload_Kamitaki2026.mhc.tsv` | 3.1 MB |
| `data/mhc_extracts/HHV6IE1A_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/HHV6IE1B_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/HHV7U14_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/HSV1mgG_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/HpyloriGroEL_BL2020.mhc.tsv` | 4.3 MB |
| `data/mhc_extracts/HpyloriOMP_BL2020.mhc.tsv` | 4.3 MB |
| `data/mhc_extracts/JCVP1_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/MCVVP1_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/ToxoSag1_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/VCAp18_BL2020.mhc.tsv` | 4.8 MB |
| `data/mhc_extracts/VZVgEgI_BL2020.mhc.tsv` | 4.4 MB |
| `data/mhc_extracts/ZEBRA_BL2020.mhc.tsv` | 4.8 MB |
| `data/mhc_extracts/mhc_dme_gwas_1.parquet` | 3.1 MB |
| `data/mhc_extracts/mhc_dme_gwas_1_infectious_onset.parquet` | 3.1 MB |
| `data/mhc_extracts/mhc_mvp.parquet` | 2.9 MB |
| `data/mvp_only_derived.parquet` | 262.9 MB |
| `results/cis_mr_top_ld_check.csv` | 1.3 KB |
| `results/drug_target_mapping.csv` | 26.2 KB |
| `results/hla_ebv_cmv_overlap.csv` | 25.2 KB |
| `results/hla_tag_lookup.csv` | 1.7 KB |
| `results/ldsc_rg_replication_infection_haem.csv` | 4.8 KB |
| `results/mr_all_results.csv` | 4.5 MB |
| `results/mr_discovery_replication_meta.csv` | 1.1 MB |
| `results/mr_forest_selected.png` | 126.2 KB |
| `results/mr_harmonised_instruments_GWAS1.csv.gz` | 304.3 KB |
| `results/rep_8loci.csv` | 2.4 KB |
| `results/rep_sign_concordance.csv` | 2.2 KB |
| `results/reverse_mr_instruments.csv` | 423 B |
| `results/reverse_mr_results.csv` | 4.2 KB |
| `results/smr_direction_magma_genes.csv` | 23.6 KB |
| `scripts/01_build_outcomes.py` | 2.6 KB |
| `scripts/02_replication.py` | 4.2 KB |
| `scripts/03_fetch_exposure.py` | 2.8 KB |
| `scripts/04_instruments.py` | 6.3 KB |
| `scripts/05_mr.py` | 7.0 KB |
| `scripts/06_ldsc_rg.py` | 6.0 KB |
| `scripts/07_hla.py` | 4.2 KB |
| `scripts/07b_fetch_mhc_controls.py` | 928 B |
| `scripts/08_smr_direction.py` | 4.0 KB |
| `scripts/09_reverse_mr.py` | 3.8 KB |
| `scripts/10_drug_targets.py` | 4.9 KB |
| `scripts/11_figures.py` | 2.2 KB |
| `scripts/mrlib.py` | 5.8 KB |
| `scripts/tabixlib.py` | 1.7 KB |
