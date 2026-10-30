# INPUTS_NOTE — microbiome_function stream

## Newly downloaded

### 1. Esteban 2025/2026 "Ramsay" ME/CFS tryptophan/AHR study (figshare)
- Source: https://api.figshare.com/v2/articles/30278434 (DOI 10.6084/m9.figshare.30278434.v2,
  "All data and analysis files MECFS AHR study", v2 2025-12-29). Paper: Europe PMC PMC13284739,
  "Tryptophan Metabolism and Aryl-Hydrocarbon Receptor Agonists in the Gut Microbiome of People With
  Myalgic Encephalomyelitis/Chronic Fatigue Syndrome" (2026).
- Path: `/mnt/project-files/mecfs-research/deep2/microbiome_function/data_esteban/` (~1.4 MB total)
- Key files:
  - `metabolites_no_outliers.csv` — 42 rows x 81 cols. **STOOL** metabolites (normalised to stool wet
    weight, per the authors' own R script), 19 patients / 23 controls, 79 metabolites including
    indole-3-propionate, indole-3-lactate, phenylpropionate, phenyllactate, 4-hydroxyphenyllactate,
    4-hydroxyphenylpropionate, cinnamoylglycine (100% missing), indole, SCFAs, bile acids.
  - `Metabolites_for_corr.csv` — 42 x 79, same metabolites, different normalisation.
  - `Relative_abudance_with_OTUID.csv` — 63 x 1191 (16S ASVs). **Surprise / caveat:** this is NOT the
    whole community. Row sums ~0.92, only 25 genera, and the abundant genera are largely
    environmental/contaminant (Mycobacterium, Streptomyces, Virgibacillus, Haloferax, a marine
    archaeon). It looks like the "tryptophan-correlated subcommunity" subset. Not usable for
    ordinary taxon inference. The full feature table is not deposited (only .qza distance matrices).
  - `taxonomy.csv` — 4223 x 9, ASV -> Silva taxonomy.
  - `Ramsay_All_Data%20-%20Demographics.csv`, `SF12_and_karnofsky_data.csv`, symptom tables.
- Sample of copied data (first 2 data columns of `metabolites_no_outliers.csv`):
```
Sample,Group,2-methylbutyrate,3-hydroxybenzoate,...
510,patient,517,0.178,...
591,patient,418,0.472,...
```

### 2. Guo et al. 2023 Cell Host Microbe supplementary tables
- Source: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10183837/supplementaryFiles
  (PMID 36758522, PMC10183837, "Deficient butyrate-producing capacity in the gut microbiome is
  associated with bacterial network disturbances and fatigue symptoms in ME/CFS", 106 ME/CFS + 91 HC,
  Columbia/Jackson Lab).
- Path: `data_esteban/NIHMS1884352-supplement-Table_S*.xlsx` (~180 KB)
- Table S2 = MaAsLin2 species results, 4 models; Model 1 (all ME/CFS vs HC, adjusted for sr-IBS, site,
  sex, BMI, race/ethnicity, age, antibiotics, pro/prebiotics) gives 42 species rows. **Caveat:** only
  species with nominally significant unadjusted p are listed, so the table is a selected subset.
  Coefficient sign is for `value = Control`, i.e. a negative coefficient means higher in ME/CFS.
- Extracted to `guo2023_species_maaslin.csv` (42 rows).
- Table S3 = GOmixer metabolic process / module results (butyrate production -53.3, p=0.019).
- Raijmakers 2020 supplement (PMC7690002) was also downloaded but is PDF/DOCX only, no usable tables.

## Reused local data (not re-downloaded)
- `/mnt/project-files/biomap/omics/Omics_Dataset/` — Metabolomics.csv (730 plasma metabolites x 414
  samples), `Specie Abundance.csv` (135 species x 479), KEGG.csv (4429 KOs x 479), Immune.csv
  (311 x 489), Quest_Lab_Measurement.csv (48 labs), Metadata.csv. 392 samples / 231 people have
  paired stool metagenome + plasma metabolome.
  - Surprise: the deposited KO table has only 4429 KOs and **does not contain the fld operon KOs**
    (K23109 fldH, K23110 fldB, K23111 fldC, K24016 acdA, K13607 fldA), so direct gene-level testing of
    the IPA pathway is impossible from these tables. K10797 (enr, 2-enoate reductase) is present and
    was used as the closest functional proxy.
  - The 135-species table is prevalence-filtered and contains no Peptostreptococcus or C. sporogenes.
- `/mnt/project-files/mecfs-research/datasets/repos/Microbe2991_.../ME&CFS_Microbiome_Results/combined_bracken_output.txt`
  — 12726 taxa x 58 samples, Kraken2/Bracken. **Not an independent cohort**: sample names are
  SAM000xxx, i.e. the same Jackson Laboratory cohort as BioMapAI (BioProject PRJNA878603). Valuable
  because it is unfiltered (includes C. sporogenes, Peptostreptococcus, Methanobrevibacter) and all
  58 samples have paired plasma metabolome. DB is bacteria/archaea only: 0 Candida, 0 crAssphage,
  essentially no phage at prevalence >= 50%.
- `/mnt/project-files/mecfs-research/datasets/repos/ohlab_MECFS_2021/` — also the same cohort
  (SAM000xxx ids), MaAsLin confounder outputs only.

## Could not obtain
- Guo 2023 per-sample species/KO tables and plasma metabolome (only summary supplements are public;
  raw reads are in SRA under the Columbia/JAX BioProjects and would need assembly + profiling).
- Raijmakers 2020 Q fever fatigue multi-omics per-sample tables (supplements are PDF/DOCX).
- Nagy-Szakal 2017, NIH Walitt PRJNA954397, UConn, Japanese NCNP: reads only, no processed tables found.
