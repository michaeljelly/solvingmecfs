# Finnish registry enrichments, ClinicalTrials.gov results and literature map

Catalogue id: `analysis-epi-literature`

**What we found (one line):** Finnish registry data (crude observed/expected) show postviral fatigue syndrome enriched for hypermobility, IgG subclass deficiency, autonomic disorders and low reticulocytes; posted trial results are almost all null; full recovery after the first year is uncommon (median 5%).

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Finnish registry enrichments, ClinicalTrials.gov results and literature map. https://data.solvingmecfs.org/analysis-epi-literature/
- **Source:** https://data.solvingmecfs.org/analysis-epi-literature/README.md
- **Sample size:** see README
- **Data type:** our analysis: registry epidemiology, trial registry, literature mining; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Raw ClinicalTrials.gov and PubMed downloads, raw Risteys scrapes, and the sentence-level abstract excerpts are excluded; the paper-level table is published without its excerpt column. One curated table that referred to an individual patient is excluded.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `clinicaltrials-gov`: https://data.solvingmecfs.org/clinicaltrials-gov/README.md (if mirrored) or see the catalogue
- `pubmed-mecfs`: https://data.solvingmecfs.org/pubmed-mecfs/README.md (if mirrored) or see the catalogue
- `risteys-finngen`: https://data.solvingmecfs.org/risteys-finngen/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-epi-literature/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-epi-literature/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `INPUTS_NOTE.md` | 5.1 KB |
| `ctgov_all_interventional_trials.csv` | 381.5 KB |
| `ctgov_curated_effects.csv` | 178.1 KB |
| `ctgov_posted_results_outcomes.csv` | 1009.0 KB |
| `ctgov_primary_summary.csv` | 28.0 KB |
| `curated/natural_history.csv` | 1.1 KB |
| `curated/registry_studies.csv` | 890 B |
| `literature/biomarker_paper_calls.csv` | 126.0 KB |
| `literature/biomarker_replication_catalogue.csv` | 7.8 KB |
| `risteys/POSTVIRFAT_overlap_enrichment.csv` | 682.7 KB |
| `risteys/POSTVIRFAT_selected_FR_vs_FG_enrichment.csv` | 5.8 KB |
| `scripts/biomarker_dict.py` | 11.2 KB |
| `scripts/biomarker_mine.py` | 8.0 KB |
| `scripts/ctgov_effects.py` | 4.3 KB |
| `scripts/ctgov_extract.py` | 4.1 KB |
| `scripts/ctgov_pull.py` | 1.1 KB |
| `scripts/pubmed_generic.py` | 1.8 KB |
| `scripts/pubmed_pull.py` | 2.6 KB |
| `scripts/risteys_codewas.py` | 2.6 KB |
| `scripts/risteys_fr_counts.py` | 1.5 KB |
| `scripts/risteys_live.py` | 3.6 KB |
| `scripts/risteys_rel.py` | 2.4 KB |
