# Blood and CSF proteomics, extracellular vesicles and cross-omics convergence

Catalogue id: `analysis-proteome`

**What we found (one line):** No single blood protein replicates strongly; muscle-derived proteins are modestly lower in three SomaScan cohorts and UK Biobank; the many Hoel/Fluge SomaScan hits are a sample-handling artefact; protein, RNA, metabolite and genetic layers do not converge.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Blood and CSF proteomics, extracellular vesicles and cross-omics convergence. https://data.solvingmecfs.org/analysis-proteome/
- **Source:** https://data.solvingmecfs.org/analysis-proteome/README.md
- **Sample size:** see README
- **Data type:** our analysis: proteomics, EVs, cross-omics; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Only per-cohort case-control statistics; no per-person rows. Raw inputs (several under CC BY-NC-ND) are not mirrored.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `bragee-2026-csf`: https://data.solvingmecfs.org/bragee-2026-csf/README.md (if mirrored) or see the catalogue
- `geo-gse214283`: https://data.solvingmecfs.org/geo-gse214283/README.md (if mirrored) or see the catalogue
- `geo-gse227375`: https://data.solvingmecfs.org/geo-gse227375/README.md (if mirrored) or see the catalogue
- `geo-gse251790`: https://data.solvingmecfs.org/geo-gse251790/README.md (if mirrored) or see the catalogue
- `geo-gse254030`: https://data.solvingmecfs.org/geo-gse254030/README.md (if mirrored) or see the catalogue
- `geo-gse293840`: https://data.solvingmecfs.org/geo-gse293840/README.md (if mirrored) or see the catalogue
- `geo-gse98139`: https://data.solvingmecfs.org/geo-gse98139/README.md (if mirrored) or see the catalogue
- `germain-2025-proteome`: https://data.solvingmecfs.org/germain-2025-proteome/README.md (if mirrored) or see the catalogue
- `giloteaux-2024-ev`: https://data.solvingmecfs.org/giloteaux-2024-ev/README.md (if mirrored) or see the catalogue
- `nih-walitt-2024`: https://data.solvingmecfs.org/nih-walitt-2024/README.md (if mirrored) or see the catalogue
- `pride-pad000026`: https://data.solvingmecfs.org/pride-pad000026/README.md (if mirrored) or see the catalogue
- `pride-pxd065439`: https://data.solvingmecfs.org/pride-pxd065439/README.md (if mirrored) or see the catalogue
- `rydland-2026-ev`: https://data.solvingmecfs.org/rydland-2026-ev/README.md (if mirrored) or see the catalogue
- `ukb-beentjes-2025`: https://data.solvingmecfs.org/ukb-beentjes-2025/README.md (if mirrored) or see the catalogue
- `zenodo-13363117`: https://data.solvingmecfs.org/zenodo-13363117/README.md (if mirrored) or see the catalogue
- `zenodo-14772494`: https://data.solvingmecfs.org/zenodo-14772494/README.md (if mirrored) or see the catalogue
- `zenodo-15237088`: https://data.solvingmecfs.org/zenodo-15237088/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-proteome/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-proteome/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `INPUTS_NOTE.md` | 4.6 KB |
| `fig_tissue_theme_heatmap.png` | 88.3 KB |
| `results/06b_log.txt` | 10.4 KB |
| `results/candidate_proteins.csv` | 24.8 KB |
| `results/cc_bragee_csf_POTS.csv` | 104.2 KB |
| `results/cc_germain_D1PRE.csv` | 743.6 KB |
| `results/cc_germain_D1PRE_nobmi.csv` | 777.4 KB |
| `results/cc_giloteaux_ev_0h.csv` | 31.0 KB |
| `results/cc_hoel_serum_leakadj.csv` | 830.4 KB |
| `results/cc_hoel_serum_raw.csv` | 810.7 KB |
| `results/cc_igg_serum.csv` | 14.0 KB |
| `results/cc_igg_serum_MEonly.csv` | 14.8 KB |
| `results/cc_lc_bansal.csv` | 20.2 KB |
| `results/cc_lc_iris_3mo.csv` | 20.0 KB |
| `results/cc_lc_uk_explore.csv` | 324.0 KB |
| `results/cc_rydland_ev.csv` | 19.3 KB |
| `results/cc_ukb_olink.csv` | 323.2 KB |
| `results/cc_walitt_csf.csv` | 153.3 KB |
| `results/cc_walitt_serum.csv` | 156.0 KB |
| `results/concordance_pairs.csv` | 2.5 KB |
| `results/convergence_kegg.csv` | 98.2 KB |
| `results/convergence_null_test.csv` | 752 B |
| `results/germain_exertion_acute_D1.csv` | 1.1 MB |
| `results/germain_exertion_acute_D2.csv` | 1.1 MB |
| `results/germain_exertion_post2_vs_post1.csv` | 1.1 MB |
| `results/germain_exertion_recovery_24h.csv` | 1.1 MB |
| `results/leak_aptamers_top100.csv` | 1.2 KB |
| `results/leak_score_other_cohorts.csv` | 196 B |
| `results/meta_blood_genes.csv` | 1.8 MB |
| `results/meta_blood_rep_null.csv` | 602 B |
| `results/plasmacell_rna_score.csv` | 300 B |
| `results/rna_GSE214283_PBMCsc.csv` | 1.5 MB |
| `results/rna_GSE227375_PBMC.csv` | 1.9 MB |
| `results/rna_GSE293840_cfRNA.csv` | 1.5 MB |
| `results/rna_GSE98139_WB.csv` | 1.5 MB |
| `results/rna_meta_genomewide.csv` | 1.9 MB |
| `results/themes_by_dataset.csv` | 35.7 KB |
| `results/themes_meta_blood.csv` | 3.6 KB |
| `results/tissue_and_theme_camera.csv` | 170.0 KB |
| `results/tissue_and_theme_camera_meta.csv` | 40.8 KB |
| `results/tissue_origin_by_dataset.csv` | 29.3 KB |
| `results/tissue_origin_by_location.csv` | 3.7 KB |
| `results/tissue_origin_drivers.csv` | 118.6 KB |
| `results/tissue_origin_meta.csv` | 9.3 KB |
| `results/tissue_scores_clinical.csv` | 12.0 KB |
| `scripts/01_load.py` | 8.3 KB |
| `scripts/02_leak.py` | 2.4 KB |
| `scripts/03_casecontrol.py` | 2.0 KB |
| `scripts/04_meta.py` | 4.0 KB |
| `scripts/05_themes.py` | 4.1 KB |
| `scripts/06_tissue.py` | 1.9 KB |
| `scripts/06b_tissue_camera.py` | 5.2 KB |
| `scripts/07_tissue_followup.py` | 2.6 KB |
| `scripts/08_scores_clinical.py` | 4.2 KB |
| `scripts/09_germain_exertion.py` | 1.4 KB |
| `scripts/10_longcovid.py` | 3.3 KB |
| `scripts/11_candidates.py` | 1.8 KB |
| `scripts/12_rna_genomewide.py` | 4.2 KB |
| `scripts/13_convergence.py` | 6.6 KB |
| `scripts/14_plasmacell_rna.py` | 2.6 KB |
| `scripts/common.py` | 3.5 KB |
| `scripts/httpfile.py` | 1.4 KB |
| `scripts/themes.py` | 3.3 KB |
