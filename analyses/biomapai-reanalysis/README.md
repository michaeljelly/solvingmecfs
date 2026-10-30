# Re-analysis of the BioMapAI (JAX) multi-omics cohort

Catalogue id: `analysis-biomapai-reanalysis`

**What we found (one line):** The gut butyrate-producer signal replicates in two external stool cohorts, but BioMapAI's reported accuracy is inflated by visit-level cross-validation, no data layer predicts individual symptom scores, and patients do not split into real subtypes.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Re-analysis of the BioMapAI (JAX) multi-omics cohort. https://data.solvingmecfs.org/analysis-biomapai-reanalysis/
- **Source:** https://data.solvingmecfs.org/analysis-biomapai-reanalysis/README.md
- **Sample size:** see README
- **Data type:** our analysis: multi-omics re-analysis; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. BioMapAI data have no redistribution licence, so all per-person or per-visit outputs (out-of-fold predictions, cluster assignments, trajectories, axis values) are excluded; only aggregate statistics, scripts and figures are published.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `biomapai-2025`: https://data.solvingmecfs.org/biomapai-2025/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-biomapai-reanalysis/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-biomapai-reanalysis/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `analysis/baseline_predictors_of_change.csv` | 398.6 KB |
| `analysis/butyrate.py` | 1.7 KB |
| `analysis/common.py` | 410 B |
| `analysis/cross_sectional_results.csv` | 6.2 MB |
| `analysis/crosssec.py` | 2.7 KB |
| `analysis/improvement.py` | 2.2 KB |
| `analysis/leakage.py` | 1.6 KB |
| `analysis/leakage_results.csv` | 1.3 KB |
| `analysis/longitudinal.py` | 3.2 KB |
| `analysis/targeted.py` | 3.4 KB |
| `analysis/targeted_composites.csv` | 9.1 KB |
| `analysis/within_person_associations.csv` | 4.1 MB |
| `deep/coupling/01_predictability.py` | 9.9 KB |
| `deep/coupling/02_diffcorr.py` | 8.1 KB |
| `deep/coupling/03_levels_symptoms.py` | 5.1 KB |
| `deep/coupling/04_immune_kyn.py` | 4.4 KB |
| `deep/coupling/04b_within_check.py` | 1.6 KB |
| `deep/coupling/05_levels_sensitivity.py` | 3.3 KB |
| `deep/coupling/06_network_counts.py` | 2.9 KB |
| `deep/coupling/07_external_replication.py` | 3.5 KB |
| `deep/coupling/08_figures.py` | 3.4 KB |
| `deep/coupling/cc.py` | 2.8 KB |
| `deep/coupling/diffcorr_global_perm.csv` | 90 B |
| `deep/coupling/diffcorr_prespecified.csv` | 18.7 KB |
| `deep/coupling/diffcorr_screen.csv` | 1.3 MB |
| `deep/coupling/external_replication.csv` | 3.7 KB |
| `deep/coupling/fig1_predictability_species_metabolome.png` | 103.5 KB |
| `deep/coupling/fig2_microbial_metabolites_levels_3cohorts.png` | 97.1 KB |
| `deep/coupling/immune_kyn_between.csv` | 93.4 KB |
| `deep/coupling/immune_kyn_within_check.csv` | 8.2 KB |
| `deep/coupling/immune_kyn_within_patients.csv` | 21.5 KB |
| `deep/coupling/levels_mixed.csv` | 10.3 KB |
| `deep/coupling/levels_sensitivity.csv` | 4.3 KB |
| `deep/coupling/log_pred_immune_kegg_all.txt` | 2.2 KB |
| `deep/coupling/log_pred_immune_species_all.txt` | 2.4 KB |
| `deep/coupling/log_pred_metab_kegg_all.txt` | 2.2 KB |
| `deep/coupling/log_pred_metab_species_Short.txt` | 2.5 KB |
| `deep/coupling/log_pred_metab_species_all.txt` | 2.4 KB |
| `deep/coupling/network_counts_matched.csv` | 1.3 KB |
| `deep/coupling/pred_immune_kegg_all.csv` | 104.7 KB |
| `deep/coupling/pred_immune_kegg_all_global.csv` | 536 B |
| `deep/coupling/pred_immune_species_all.csv` | 104.7 KB |
| `deep/coupling/pred_immune_species_all_global.csv` | 540 B |
| `deep/coupling/pred_metab_kegg_all.csv` | 299.9 KB |
| `deep/coupling/pred_metab_kegg_all_global.csv` | 577 B |
| `deep/coupling/pred_metab_species_Short.csv` | 299.6 KB |
| `deep/coupling/pred_metab_species_Short_global.csv` | 573 B |
| `deep/coupling/pred_metab_species_all.csv` | 300.3 KB |
| `deep/coupling/pred_metab_species_all_global.csv` | 563 B |
| `deep/coupling/pred_microbial_set_test.csv` | 563 B |
| `deep/coupling/run1.sh` | 352 B |
| `deep/coupling/run2.sh` | 270 B |
| `deep/coupling/run3.sh` | 262 B |
| `deep/coupling/symptoms_between_within.csv` | 113.9 KB |
| `deep/models/biomapai_dnn.py` | 2.1 KB |
| `deep/models/common.py` | 3.8 KB |
| `deep/models/fig_auc_comparison.png` | 67.9 KB |
| `deep/models/logs/s1.log` | 214 B |
| `deep/models/logs/s1_specie.log` | 214 B |
| `deep/models/logs/s1b.log` | 715 B |
| `deep/models/logs/s2.log` | 3.0 KB |
| `deep/models/logs/s2_metabqc.log` | 196 B |
| `deep/models/logs/s2b.log` | 12.2 KB |
| `deep/models/logs/s3.log` | 2.9 KB |
| `deep/models/logs/s3b.log` | 444 B |
| `deep/models/logs/s4b.log` | 1.0 KB |
| `deep/models/logs/s4perm.log` | 100 B |
| `deep/models/logs/s4w.log` | 1.4 KB |
| `deep/models/logs/s5a.log` | 13.5 KB |
| `deep/models/logs/s5b.log` | 7.7 KB |
| `deep/models/make_figure.py` | 2.6 KB |
| `deep/models/prep.py` | 5.0 KB |
| `deep/models/results/s1_summary.csv` | 2.4 KB |
| `deep/models/results/s2_summary_all_short_long.csv` | 15.6 KB |
| `deep/models/results/s3_external.csv` | 9.4 KB |
| `deep/models/results/s3_external_v1_oldpseudo.csv` | 9.4 KB |
| `deep/models/results/s3b_external_adjtable.csv` | 1.3 KB |
| `deep/models/results/s4_between.csv` | 21.1 KB |
| `deep/models/results/s4_between_metabqc.csv` | 6.5 KB |
| `deep/models/results/s4_perm_kegg_lgbm.csv` | 1.6 KB |
| `deep/models/results/s4_within.csv` | 19.0 KB |
| `deep/models/results/s4_within_metabqc.csv` | 6.3 KB |
| `deep/models/results/s5_importance_all_immune_rf.csv` | 30.9 KB |
| `deep/models/results/s5_importance_all_metabolome_lgbm.csv` | 70.8 KB |
| `deep/models/results/s5_importance_all_metabqc_enet.csv` | 88.0 KB |
| `deep/models/results/s5_importance_all_metabqc_lgbm.csv` | 85.7 KB |
| `deep/models/results/s5_importance_all_specie_rf.csv` | 11.1 KB |
| `deep/models/results/s5_importance_short_metabolome_lgbm.csv` | 70.7 KB |
| `deep/models/results/s5_importance_short_metabqc_lgbm.csv` | 91.1 KB |
| `deep/models/results/s5_importance_short_specie_rf.csv` | 11.1 KB |
| `deep/models/run_s2.sh` | 278 B |
| `deep/models/run_s2b.sh` | 228 B |
| `deep/models/run_s4.sh` | 148 B |
| `deep/models/run_s5a.sh` | 311 B |
| `deep/models/run_s5b.sh` | 246 B |
| `deep/models/s1_biomapai_cv.py` | 2.2 KB |
| `deep/models/s1_summary.py` | 1.9 KB |
| `deep/models/s2_classify.py` | 2.7 KB |
| `deep/models/s2_fusion.py` | 4.5 KB |
| `deep/models/s3_external.py` | 11.4 KB |
| `deep/models/s3b_external_adjtable.py` | 1.8 KB |
| `deep/models/s4_perm.py` | 1.6 KB |
| `deep/models/s4_symptoms.py` | 5.6 KB |
| `deep/models/s5_importance.py` | 2.1 KB |
| `deep/stats_replication/01_casecontrol.py` | 5.5 KB |
| `deep/stats_replication/02_symptoms.py` | 5.0 KB |
| `deep/stats_replication/03_replication.py` | 16.7 KB |
| `deep/stats_replication/04_exploratory_renal.py` | 4.2 KB |
| `deep/stats_replication/05_figures.py` | 4.9 KB |
| `deep/stats_replication/casecontrol_all_models.csv.gz` | 7.8 MB |
| `deep/stats_replication/casecontrol_n_significant.csv` | 917 B |
| `deep/stats_replication/casecontrol_sample_counts.csv` | 1.1 KB |
| `deep/stats_replication/data.py` | 4.2 KB |
| `deep/stats_replication/exploratory_filtration_markers.csv` | 2.4 KB |
| `deep/stats_replication/fig_replication_forest.png` | 122.8 KB |
| `deep/stats_replication/fig_symptom_between_within.png` | 60.6 KB |
| `deep/stats_replication/germain_bell_severity.csv` | 4.4 KB |
| `deep/stats_replication/lmm_vs_mwu_overlap.csv` | 813 B |
| `deep/stats_replication/lmmfast.py` | 4.2 KB |
| `deep/stats_replication/replication_per_cohort.csv` | 10.1 KB |
| `deep/stats_replication/replication_table.csv` | 25.4 KB |
| `deep/stats_replication/sex_concordance.csv` | 684 B |
| `deep/stats_replication/symptom_composites_by_sex.csv` | 1.6 KB |
| `deep/stats_replication/symptom_mixed_models.csv.gz` | 2.1 MB |
| `deep/stats_replication/symptom_targeted_composites.csv` | 17.0 KB |
| `deep/stats_replication/validate_lmm.py` | 1.5 KB |
| `deep/subtypes/axes.py` | 2.4 KB |
| `deep/subtypes/clust_tools.py` | 1.8 KB |
| `deep/subtypes/fig/fig1_symptom_clusters.png` | 135.3 KB |
| `deep/subtypes/fig/fig2_omics_clusters.png` | 119.4 KB |
| `deep/subtypes/fig/fig3_candidate_axes.png` | 103.1 KB |
| `deep/subtypes/fig/fig4_temporal.png` | 82.1 KB |
| `deep/subtypes/omics_space.py` | 1.5 KB |
| `deep/subtypes/results/axes_case_control.csv` | 2.7 KB |
| `deep/subtypes/results/axes_correlations.csv` | 18.1 KB |
| `deep/subtypes/results/axes_unimodality.csv` | 1.1 KB |
| `deep/subtypes/results/axes_vs_symptoms.csv` | 2.1 KB |
| `deep/subtypes/results/omics_cluster_selection.csv` | 5.3 KB |
| `deep/subtypes/results/omics_k2_clinical.csv` | 1.6 KB |
| `deep/subtypes/results/omics_k2_features.csv` | 130.8 KB |
| `deep/subtypes/results/omics_k3_clinical.csv` | 2.3 KB |
| `deep/subtypes/results/omics_k3_features.csv` | 154.7 KB |
| `deep/subtypes/results/symptom_cluster_copula_null.csv` | 1.2 KB |
| `deep/subtypes/results/symptom_cluster_selection.csv` | 3.6 KB |
| `deep/subtypes/results/symptom_partition_omics_tests.csv` | 829.7 KB |
| `deep/subtypes/results/symptom_partition_temporal.csv` | 371 B |
| `deep/subtypes/results/symptom_pca_loadings.csv` | 1.2 KB |
| `deep/subtypes/s1_symptom_clusters.py` | 1.8 KB |
| `deep/subtypes/s1b_copula_null.py` | 1.4 KB |
| `deep/subtypes/s2_symptom_omics.py` | 3.0 KB |
| `deep/subtypes/s3_omics_clusters.py` | 3.4 KB |
| `deep/subtypes/s3b_omics_characterise.py` | 5.0 KB |
| `deep/subtypes/s3c_gut_and_snf.py` | 3.2 KB |
| `deep/subtypes/s4a_symptom_temporal.py` | 4.2 KB |
| `deep/subtypes/s5_axes.py` | 4.4 KB |
| `deep/subtypes/s7_figures.py` | 7.7 KB |
| `deep/subtypes/sub_common.py` | 2.1 KB |
