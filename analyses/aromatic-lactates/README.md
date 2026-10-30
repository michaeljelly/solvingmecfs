# Metabolome-wide meta-analysis across 10 blood cohorts (aromatic lactates)

Catalogue id: `analysis-aromatic-lactates`

**What we found (one line):** Four aromatic lactates (phenyllactate, 4-hydroxyphenyllactate, indolelactate, imidazole lactate) are lower in patients (-0.53 SD, 5 of 5 cohorts, p=4e-11) after adjusting for drugs, diet, caffeine and kidney function; collagen-breakdown markers are higher; nothing replicated tracks severity.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Metabolome-wide meta-analysis across 10 blood cohorts (aromatic lactates). https://data.solvingmecfs.org/analysis-aromatic-lactates/
- **Source:** https://data.solvingmecfs.org/analysis-aromatic-lactates/README.md
- **Sample size:** see README
- **Data type:** our analysis: metabolomics meta-analysis; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. All tables are per-study or pooled estimates; no per-person rows.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `biomapai-2025`: https://data.solvingmecfs.org/biomapai-2025/README.md (if mirrored) or see the catalogue
- `germain-2018`: https://data.solvingmecfs.org/germain-2018/README.md (if mirrored) or see the catalogue
- `germain-2020`: https://data.solvingmecfs.org/germain-2020/README.md (if mirrored) or see the catalogue
- `germain-2022`: https://data.solvingmecfs.org/germain-2022/README.md (if mirrored) or see the catalogue
- `hoel-2021`: https://data.solvingmecfs.org/hoel-2021/README.md (if mirrored) or see the catalogue
- `mwb-st000450`: https://data.solvingmecfs.org/mwb-st000450/README.md (if mirrored) or see the catalogue
- `mwb-st000617`: https://data.solvingmecfs.org/mwb-st000617/README.md (if mirrored) or see the catalogue
- `mwb-st000800`: https://data.solvingmecfs.org/mwb-st000800/README.md (if mirrored) or see the catalogue
- `mwb-st000910`: https://data.solvingmecfs.org/mwb-st000910/README.md (if mirrored) or see the catalogue
- `mwb-st002000`: https://data.solvingmecfs.org/mwb-st002000/README.md (if mirrored) or see the catalogue
- `mwb-st003103`: https://data.solvingmecfs.org/mwb-st003103/README.md (if mirrored) or see the catalogue
- `mwb-st003178`: https://data.solvingmecfs.org/mwb-st003178/README.md (if mirrored) or see the catalogue
- `mwb-st003346`: https://data.solvingmecfs.org/mwb-st003346/README.md (if mirrored) or see the catalogue
- `mwb-st004940`: https://data.solvingmecfs.org/mwb-st004940/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-aromatic-lactates/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-aromatic-lactates/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `10_build_all.py` | 12.9 KB |
| `11_harmonise.py` | 8.4 KB |
| `12_assoc_meta.py` | 4.6 KB |
| `13_confounders.py` | 6.1 KB |
| `14_ratios_modules.py` | 10.8 KB |
| `15_gut_species.py` | 9.9 KB |
| `16_key_sensitivity.py` | 4.0 KB |
| `17_severity.py` | 5.1 KB |
| `18_classifier.py` | 4.0 KB |
| `19_figure.py` | 2.2 KB |
| `INPUTS_NOTE.md` | 1.6 KB |
| `adjusted_meta.csv` | 881.6 KB |
| `classifier_transfer.csv` | 18.9 KB |
| `cometabolite_correlations.csv` | 8.0 KB |
| `drug_exposure_by_group.csv` | 4.7 KB |
| `fig_forest_key.png` | 157.0 KB |
| `harmonised_map.csv` | 597.9 KB |
| `key_items_per_study.csv` | 86.1 KB |
| `key_items_sensitivity.csv` | 13.0 KB |
| `meta_all.csv` | 1.3 MB |
| `mwlib.py` | 6.9 KB |
| `per_study_all.csv` | 1.0 MB |
| `ratios_modules_meta.csv` | 86.0 KB |
| `ratios_per_study.csv` | 213.1 KB |
| `replication_designs.csv` | 21.6 KB |
| `severity_key_items.csv` | 51.3 KB |
| `severity_metabolome_wide.csv` | 144.2 KB |
| `species_cv_r2.csv` | 1.3 KB |
| `species_metabolite_partial.csv` | 951.8 KB |
| `taxa_candidates.csv` | 39.8 KB |
| `taxa_mediation.csv` | 954 B |
