# Gut microbiome function behind low IPA and indolelactate

Catalogue id: `analysis-microbiome-ipa`

**What we found (one line):** In 231 people with paired stool and blood, gut bacteria predict plasma IPA and account for its patient-control gap but do not predict indolelactate; low IPA goes with an overgrowth of R. gnavus, F. plautii and R. lactatiformans that also appears in an independent US cohort.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Gut microbiome function behind low IPA and indolelactate. https://data.solvingmecfs.org/analysis-microbiome-ipa/
- **Source:** https://data.solvingmecfs.org/analysis-microbiome-ipa/README.md
- **Sample size:** see README
- **Data type:** our analysis: microbiome function, metabolomics; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Per-person BioMapAI panel values, scores and Bracken guild abundances are excluded (sources have no redistribution licence). Copies of Esteban and Guo files are under figshare-esteban-2026 and guo-2023.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `biomapai-2025`: https://data.solvingmecfs.org/biomapai-2025/README.md (if mirrored) or see the catalogue
- `figshare-esteban-2026`: https://data.solvingmecfs.org/figshare-esteban-2026/README.md (if mirrored) or see the catalogue
- `github-microbe2991-myalgic-encephalomyelitis-chronic-fatigue`: https://data.solvingmecfs.org/github-microbe2991-myalgic-encephalomyelitis-chronic-fatigue/README.md (if mirrored) or see the catalogue
- `guo-2023`: https://data.solvingmecfs.org/guo-2023/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-microbiome-ipa/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-microbiome-ipa/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `INPUTS_NOTE.md` | 4.8 KB |
| `biomap_arldh_gene_panel.csv` | 11.0 KB |
| `biomap_aromatic_panel.csv` | 2.2 KB |
| `biomap_cv_prediction.csv` | 385 B |
| `biomap_immune_vs_modules.csv` | 57.4 KB |
| `biomap_kegg_assoc.csv` | 1.3 MB |
| `biomap_labs_vs_modules.csv` | 7.7 KB |
| `biomap_lactate_vs_parent_aa.csv` | 1.4 KB |
| `biomap_panel_corr_controls.csv` | 7.5 KB |
| `biomap_renal_adjustment.csv` | 832 B |
| `biomap_species_assoc.csv` | 36.8 KB |
| `biomap_symptom_assoc_patients.csv` | 19.0 KB |
| `biomap_vs_guo_species_concordance.csv` | 6.7 KB |
| `bracken_guilds.csv` | 1.2 KB |
| `bracken_taxon_scan.csv` | 499.0 KB |
| `esteban_family_metab_corr.csv` | 9.7 KB |
| `esteban_genus_metab_corr.csv` | 10.1 KB |
| `esteban_module_corr.csv` | 1.1 KB |
| `esteban_stool_metabolites.csv` | 1.4 KB |
| `guo2023_species_maaslin.csv` | 2.1 KB |
| `scripts/a1_panel.py` | 3.4 KB |
| `scripts/a2_ratio.py` | 3.2 KB |
| `scripts/a3_micro.py` | 3.3 KB |
| `scripts/a4_symptoms.py` | 3.6 KB |
| `scripts/a5_predict.py` | 3.0 KB |
| `scripts/b1_esteban.py` | 2.3 KB |
| `scripts/b2_esteban_taxa.py` | 2.1 KB |
| `scripts/c1_concord.py` | 1.4 KB |
| `scripts/d1_bracken.py` | 5.0 KB |
| `scripts/d2_bracken_scan.py` | 3.1 KB |
| `scripts/e1_arldh_genes.py` | 4.3 KB |
| `scripts/load.py` | 771 B |
