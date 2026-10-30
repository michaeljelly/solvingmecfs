# Post-exertional response re-analysis of the Cornell 2-day CPET metabolomics

Catalogue id: `analysis-pem-germain2022`

**What we found (one line):** Patients' plasma metabolome differs at rest, but within-person changes after exercise barely differ from controls; the paper's growing divergence over two days does not hold within people; lactate and glucose rises are blunted in patients.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Post-exertional response re-analysis of the Cornell 2-day CPET metabolomics. https://data.solvingmecfs.org/analysis-pem-germain2022/
- **Source:** https://data.solvingmecfs.org/analysis-pem-germain2022/README.md
- **Sample size:** see README
- **Data type:** our analysis: metabolomics, exercise response; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Pickled working copies of the per-person data are excluded (the source xlsx is under germain-2022).

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `germain-2022`: https://data.solvingmecfs.org/germain-2022/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-pem-germain2022/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-pem-germain2022/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `00_load.py` | 1.5 KB |
| `01_contrasts.py` | 3.9 KB |
| `01_output.txt` | 24.0 KB |
| `02_output.txt` | 9.0 KB |
| `02_perm_sex_bell.py` | 3.2 KB |
| `03_output.txt` | 62.4 KB |
| `03_sets_leads_sex_bell.py` | 6.1 KB |
| `04_output.txt` | 4.4 KB |
| `04_replicate_paper.py` | 3.1 KB |
| `05_mixed.py` | 2.0 KB |
| `05_output.txt` | 12.0 KB |
| `06_output.txt` | 843 B |
| `06_subgroups.py` | 4.7 KB |
| `07_effort_adjusted.py` | 2.2 KB |
| `07_output.txt` | 5.7 KB |
| `08_figure.py` | 2.3 KB |
| `contrast_results.xlsx` | 905.4 KB |
| `effort_adjusted_response.csv` | 363.3 KB |
| `mixed_model_results.csv` | 218.0 KB |
| `pem_figure.png` | 210.2 KB |
| `permutation_global.csv` | 548 B |
| `response_reproducibility_R1_vs_R2.csv` | 7.5 KB |
| `resting_mean_case_vs_ctrl.csv` | 94.7 KB |
| `resting_subpathway_shift.csv` | 5.3 KB |
| `set_scores_results.csv` | 309.6 KB |
| `wanted_urls.txt` | 1.7 KB |
