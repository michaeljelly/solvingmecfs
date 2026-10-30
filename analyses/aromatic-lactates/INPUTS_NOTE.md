# INPUTS_NOTE: metabolome_wide stream

No new datasets were downloaded in this stream. Everything was built from local copies:

- `/mnt/project-files/mecfs-research/downloads/mwb/raw/` (Metabolomics Workbench mwTab/Results files for ST000450, ST000800/801/813, ST004940; MetaboLights MTBLS161; Europe PMC supplements for Germain 2018 and 2020; Hoel 2021 JCI Insight sdd1.xlsx). These are described in `downloads/mwb/INPUTS_NOTE.md`.
- `/mnt/project-files/biomap/BioMapAI/data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv` (QC-normalised Metabolon, 515 visits), `data/metagenomics/metaphlan4_species_full.csv`, `data/indenpendent_cohorts/` (Che 2022, Germain 2022 JCI Insight).

New derived table worth indexing: `harmonised_map.csv` (CSV, 5,764 rows x 9 columns): one row per (harmonised metabolite group `hkey`, study, study-specific feature name), with the identifier types that linked it. 2,344 groups; 1,117 seen in >=2 cohorts, 480 in >=5.

```text (copied data, harmonised_map.csv)
hkey,label,study,feature,n_group,links,super,sub,k_studies
indolepropionate -> BioMapAI: indolepropionate | Che2022: indole_3_propionic_acid_PM_a1 | NCNP2026: Indole-3-propionic acid
1-palmitoyl-2-oleoyl-GPC (16:0/18:1) -> Che2022: PC(34:1) | NagySzakal2018: PC (34:1)   (sum-composition lipid match)
```

Surprises: in Hoel 2021 sdd1 the xenobiotics sit on a separate sheet (`11.HD4 Xenobiotics`); Che 2022 lists some compounds twice (GC and LC platforms, `_PM_a1` suffix), resolved by keeping the first listed; about 85-95% of people in every Metabolon cohort have at least one drug metabolite detected, mostly salicylate/ibuprofen/paracetamol traces.
