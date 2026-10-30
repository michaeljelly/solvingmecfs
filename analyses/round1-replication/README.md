# Round 1: public ME/CFS cohorts re-analysed (metabolomics, NIH Walitt, transcriptomics)

Catalogue id: `analysis-round1-replication`

**What we found (one line):** Indolelactate is lower in patients in 7 of 7 blood cohorts (pooled -0.39 SD), IPA modestly lower, kynurenine lower rather than higher; bile acids, VMA, ceramides and Th17 do not replicate, and no blood gene set passes FDR.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Round 1: public ME/CFS cohorts re-analysed (metabolomics, NIH Walitt, transcriptomics). https://data.solvingmecfs.org/analysis-round1-replication/
- **Source:** https://data.solvingmecfs.org/analysis-round1-replication/README.md
- **Sample size:** see README
- **Data type:** our analysis: metabolomics, proteomics, transcriptomics, microbiome (re-analysis); scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Excludes the pickled harmonised per-person Metabolomics Workbench table (CC BY-NC-ND source) and a mapMECFS catalogue table. Person-level tidy tables are included only for the NIH Walitt Source Data (CC BY 4.0).

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `baraniuk-2025-csf`: https://data.solvingmecfs.org/baraniuk-2025-csf/README.md (if mirrored) or see the catalogue
- `biomapai-2025`: https://data.solvingmecfs.org/biomapai-2025/README.md (if mirrored) or see the catalogue
- `figshare-esteban-2026`: https://data.solvingmecfs.org/figshare-esteban-2026/README.md (if mirrored) or see the catalogue
- `geo-gse128078`: https://data.solvingmecfs.org/geo-gse128078/README.md (if mirrored) or see the catalogue
- `geo-gse16059`: https://data.solvingmecfs.org/geo-gse16059/README.md (if mirrored) or see the catalogue
- `geo-gse214283`: https://data.solvingmecfs.org/geo-gse214283/README.md (if mirrored) or see the catalogue
- `geo-gse227375`: https://data.solvingmecfs.org/geo-gse227375/README.md (if mirrored) or see the catalogue
- `geo-gse251790`: https://data.solvingmecfs.org/geo-gse251790/README.md (if mirrored) or see the catalogue
- `geo-gse254030`: https://data.solvingmecfs.org/geo-gse254030/README.md (if mirrored) or see the catalogue
- `geo-gse293840`: https://data.solvingmecfs.org/geo-gse293840/README.md (if mirrored) or see the catalogue
- `geo-gse98139`: https://data.solvingmecfs.org/geo-gse98139/README.md (if mirrored) or see the catalogue
- `germain-2018`: https://data.solvingmecfs.org/germain-2018/README.md (if mirrored) or see the catalogue
- `germain-2020`: https://data.solvingmecfs.org/germain-2020/README.md (if mirrored) or see the catalogue
- `germain-2022`: https://data.solvingmecfs.org/germain-2022/README.md (if mirrored) or see the catalogue
- `github-nadjoubiii-myalgic-encephalomyelitis-project`: https://data.solvingmecfs.org/github-nadjoubiii-myalgic-encephalomyelitis-project/README.md (if mirrored) or see the catalogue
- `hoel-2021`: https://data.solvingmecfs.org/hoel-2021/README.md (if mirrored) or see the catalogue
- `metabolights-mtbls161`: https://data.solvingmecfs.org/metabolights-mtbls161/README.md (if mirrored) or see the catalogue
- `mwb-st000450`: https://data.solvingmecfs.org/mwb-st000450/README.md (if mirrored) or see the catalogue
- `mwb-st000617`: https://data.solvingmecfs.org/mwb-st000617/README.md (if mirrored) or see the catalogue
- `mwb-st000800`: https://data.solvingmecfs.org/mwb-st000800/README.md (if mirrored) or see the catalogue
- `mwb-st000910`: https://data.solvingmecfs.org/mwb-st000910/README.md (if mirrored) or see the catalogue
- `mwb-st002000`: https://data.solvingmecfs.org/mwb-st002000/README.md (if mirrored) or see the catalogue
- `mwb-st003103`: https://data.solvingmecfs.org/mwb-st003103/README.md (if mirrored) or see the catalogue
- `mwb-st003178`: https://data.solvingmecfs.org/mwb-st003178/README.md (if mirrored) or see the catalogue
- `mwb-st003346`: https://data.solvingmecfs.org/mwb-st003346/README.md (if mirrored) or see the catalogue
- `mwb-st004940`: https://data.solvingmecfs.org/mwb-st004940/README.md (if mirrored) or see the catalogue
- `nih-walitt-2024`: https://data.solvingmecfs.org/nih-walitt-2024/README.md (if mirrored) or see the catalogue
- `nih-walitt-2024-code`: https://data.solvingmecfs.org/nih-walitt-2024-code/README.md (if mirrored) or see the catalogue
- `pride-pad000026`: https://data.solvingmecfs.org/pride-pad000026/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-round1-replication/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-round1-replication/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `datasets/01_gut_nih_replication.py` | 7.7 KB |
| `datasets/02_plasma_replication.py` | 12.9 KB |
| `datasets/03_steroids_crosscheck.py` | 6.0 KB |
| `datasets/04_nih_proteomics_csf.py` | 3.9 KB |
| `datasets/05_mtbls161_nmr.py` | 1.5 KB |
| `datasets/06_figure.py` | 2.2 KB |
| `datasets/results/fig_ipa_gut_replication.png` | 132.9 KB |
| `datasets/results/germain2020_drug_label_check.csv` | 3.1 KB |
| `datasets/results/gut_3cohort_meta.csv` | 22.3 KB |
| `datasets/results/gut_nih_prespecified.csv` | 1.7 KB |
| `datasets/results/gut_signature_auc.txt` | 218 B |
| `datasets/results/mtbls161_nmr.csv` | 5.6 KB |
| `datasets/results/nih_csf_catecholamines.csv` | 625 B |
| `datasets/results/nih_somascan_plasma_csf.csv` | 277.4 KB |
| `datasets/results/plasma_allshared_G2020_vs_G2022.csv` | 116.7 KB |
| `datasets/results/plasma_targets_G2020_G2022_meta.csv` | 13.0 KB |
| `datasets/results/plasma_transfer_auc.csv` | 113 B |
| `datasets/results/steroids_crosscheck.csv` | 4.2 KB |
| `downloads/geo/INPUTS_NOTE.md` | 5.0 KB |
| `downloads/geo/data/gds_hits.json` | 109.9 KB |
| `downloads/geo/data/gds_hits2.json` | 83.8 KB |
| `downloads/geo/data/otherrepos.txt` | 22.7 KB |
| `downloads/geo/results/DATASET_CATALOG.csv` | 5.5 KB |
| `downloads/geo/results/GSE214283_celltype_props.csv` | 2.7 KB |
| `downloads/geo/results/GSE214283_change_case_vs_ctrl.csv` | 153.6 KB |
| `downloads/geo/results/GSE214283_change_within_group.csv` | 237.8 KB |
| `downloads/geo/results/GSE214283_genes_case_control.csv` | 259.3 KB |
| `downloads/geo/results/GSE214283_sets_case_control.csv` | 37.2 KB |
| `downloads/geo/results/GSE227375_change_case_vs_ctrl.csv` | 66.4 KB |
| `downloads/geo/results/GSE227375_change_within_group.csv` | 102.6 KB |
| `downloads/geo/results/GSE227375_genes_case_control.csv` | 56.0 KB |
| `downloads/geo/results/GSE227375_sets_case_control.csv` | 7.9 KB |
| `downloads/geo/results/GSE293840_genes_case_control.csv` | 16.5 KB |
| `downloads/geo/results/GSE293840_genes_case_control_globadj.csv` | 19.4 KB |
| `downloads/geo/results/GSE293840_sets_case_control.csv` | 2.6 KB |
| `downloads/geo/results/GSE293840_sets_case_control_globadj.csv` | 3.0 KB |
| `downloads/geo/results/META_genes_blood_transcriptome.csv` | 8.7 KB |
| `downloads/geo/results/META_sets_blood_transcriptome.csv` | 1.3 KB |
| `downloads/geo/results/PAD000026_all_aptamers_case_control.csv` | 925.2 KB |
| `downloads/geo/results/PAD000026_proteins_case_control.csv` | 11.8 KB |
| `downloads/geo/results/PAD000026_proteins_case_control_leakadj.csv` | 13.6 KB |
| `downloads/geo/results/PAD000026_sets_case_control.csv` | 2.3 KB |
| `downloads/geo/results/PAD000026_sets_case_control_leakadj.csv` | 2.6 KB |
| `downloads/geo/results/PAD000026_sets_vs_SF36PF_cases.csv` | 1.1 KB |
| `downloads/geo/results/SUMMARY_genes_baseline_by_dataset.csv` | 3.0 KB |
| `downloads/geo/results/SUMMARY_sets_baseline_by_dataset.csv` | 1.4 KB |
| `downloads/geo/results/fig_sets_by_dataset.png` | 204.8 KB |
| `downloads/geo/results/small_genes_case_control.csv` | 37.5 KB |
| `downloads/geo/results/small_paired.csv` | 99.1 KB |
| `downloads/geo/results/small_sets_case_control.csv` | 5.2 KB |
| `downloads/geo/scripts/01_pseudobulk_GSE214283.py` | 2.8 KB |
| `downloads/geo/scripts/02_prep_GSE227375.py` | 2.4 KB |
| `downloads/geo/scripts/03_pairs_GSE227375.py` | 1.1 KB |
| `downloads/geo/scripts/04_analyze_GSE214283.py` | 2.8 KB |
| `downloads/geo/scripts/05_analyze_GSE227375.py` | 1.5 KB |
| `downloads/geo/scripts/06_analyze_GSE293840.py` | 2.6 KB |
| `downloads/geo/scripts/07_analyze_PAD000026.py` | 3.9 KB |
| `downloads/geo/scripts/08_analyze_small.py` | 4.3 KB |
| `downloads/geo/scripts/09_summary.py` | 4.4 KB |
| `downloads/geo/scripts/10_figure.py` | 1.1 KB |
| `downloads/geo/scripts/common.py` | 7.3 KB |
| `downloads/mwb/01_build.py` | 22.4 KB |
| `downloads/mwb/02_effects.py` | 9.9 KB |
| `downloads/mwb/03_sensitivity.py` | 5.7 KB |
| `downloads/mwb/04_figure.py` | 2.3 KB |
| `downloads/mwb/INPUTS_NOTE.md` | 5.5 KB |
| `downloads/mwb/exertion_csf.csv` | 625 B |
| `downloads/mwb/fig_forest_leads.png` | 325.6 KB |
| `downloads/mwb/meta_blood.csv` | 16.2 KB |
| `downloads/mwb/meta_csf_stool.csv` | 4.8 KB |
| `downloads/mwb/mwparse.py` | 2.2 KB |
| `downloads/mwb/name_mapping.csv` | 41.9 KB |
| `downloads/mwb/per_study_effects.csv` | 43.3 KB |
| `downloads/mwb/sensitivity.csv` | 15.6 KB |
| `downloads/nih_walitt/INPUTS_NOTE.md` | 5.4 KB |
| `downloads/nih_walitt/figures/fig1_csf_metabolome_global_shift.png` | 79.9 KB |
| `downloads/nih_walitt/figures/fig2_tilt_catecholamines.png` | 223.9 KB |
| `downloads/nih_walitt/figures/fig3_effect_sizes.png` | 131.8 KB |
| `downloads/nih_walitt/results/02_reproduce_log.txt` | 22.5 KB |
| `downloads/nih_walitt/results/03_metabolomics_log.txt` | 20.0 KB |
| `downloads/nih_walitt/results/04_somalogic_log.txt` | 3.6 KB |
| `downloads/nih_walitt/results/A_scalar_group_tests.csv` | 6.7 KB |
| `downloads/nih_walitt/results/B_catechol_behaviour_correlations.csv` | 9.3 KB |
| `downloads/nih_walitt/results/C_tilt_group_tests.csv` | 2.4 KB |
| `downloads/nih_walitt/results/C_tilt_person_level_derived.csv` | 9.0 KB |
| `downloads/nih_walitt/results/D_flow_cytometry_tests.csv` | 13.5 KB |
| `downloads/nih_walitt/results/E_eefrt_person_level.csv` | 2.2 KB |
| `downloads/nih_walitt/results/F_SD14A_check.csv` | 22.8 KB |
| `downloads/nih_walitt/results/F_csf_lead_metabolites.csv` | 15.8 KB |
| `downloads/nih_walitt/results/F_csf_metabolomics_all_tests.csv` | 182.4 KB |
| `downloads/nih_walitt/results/G_plasma_lipid_classes.csv` | 1.2 KB |
| `downloads/nih_walitt/results/G_plasma_lipid_species.csv` | 52.8 KB |
| `downloads/nih_walitt/results/H_stool_nmr.csv` | 2.3 KB |
| `downloads/nih_walitt/results/I_somalogic_all_tests.csv` | 203.6 KB |
| `downloads/nih_walitt/results/I_somalogic_lead_proteins.csv` | 2.8 KB |
| `downloads/nih_walitt/results/tidy/csf_metabolomics_paper_stats_fig6I.csv` | 33.3 KB |
| `downloads/nih_walitt/results/tidy/csf_metabolomics_wide.csv` | 310.5 KB |
| `downloads/nih_walitt/results/tidy/eefrt_trials.csv` | 114.2 KB |
| `downloads/nih_walitt/results/tidy/flow_cytometry_long.csv` | 133.8 KB |
| `downloads/nih_walitt/results/tidy/paired_measures_long.csv` | 48.0 KB |
| `downloads/nih_walitt/results/tidy/plasma_lipid_classes_long.csv` | 8.8 KB |
| `downloads/nih_walitt/results/tidy/plasma_lipidomics_species_wide.csv` | 167.6 KB |
| `downloads/nih_walitt/results/tidy/scalar_measures_long.csv` | 52.0 KB |
| `downloads/nih_walitt/results/tidy/stool_nmr_metabolites_long.csv` | 28.0 KB |
| `downloads/nih_walitt/results/tidy/tilt_plasma_catecholamines_wide.csv` | 5.2 KB |
| `downloads/nih_walitt/scripts/01_build_tidy.py` | 7.3 KB |
| `downloads/nih_walitt/scripts/02_reproduce.py` | 9.2 KB |
| `downloads/nih_walitt/scripts/03_metabolomics_leads.py` | 8.5 KB |
| `downloads/nih_walitt/scripts/04_somalogic.py` | 2.8 KB |
| `downloads/nih_walitt/scripts/05_figures.py` | 3.6 KB |
