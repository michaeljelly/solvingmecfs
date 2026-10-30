# Single-cell PBMC, DNA methylation and B/T-cell repertoire re-analysis

Catalogue id: `analysis-singlecell-methylation`

**What we found (one line):** No immune-cell type, gene or methylation site separates ME/CFS from controls reproducibly; the de Vega 2017 methylation signal is explained by granulocyte contamination, Sato 2021 IGHV3-30/3-33 overuse does not replicate, and epigenetic ageing is normal.

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). Single-cell PBMC, DNA methylation and B/T-cell repertoire re-analysis. https://data.solvingmecfs.org/analysis-singlecell-methylation/
- **Source:** https://data.solvingmecfs.org/analysis-singlecell-methylation/README.md
- **Sample size:** see README
- **Data type:** our analysis: scRNA-seq, DNA methylation, immune repertoires; scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. Per-cell and pseudobulk tables are derived from GEO GSE214283 (open). Repertoire per-sample metrics derive from ENA/SRA PRJEB83659 and PRJNA940635.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `ena-prjeb83659`: https://data.solvingmecfs.org/ena-prjeb83659/README.md (if mirrored) or see the catalogue
- `geo-gse111183`: https://data.solvingmecfs.org/geo-gse111183/README.md (if mirrored) or see the catalogue
- `geo-gse156792`: https://data.solvingmecfs.org/geo-gse156792/README.md (if mirrored) or see the catalogue
- `geo-gse214283`: https://data.solvingmecfs.org/geo-gse214283/README.md (if mirrored) or see the catalogue
- `geo-gse268212`: https://data.solvingmecfs.org/geo-gse268212/README.md (if mirrored) or see the catalogue
- `geo-gse270045`: https://data.solvingmecfs.org/geo-gse270045/README.md (if mirrored) or see the catalogue
- `geo-gse327255`: https://data.solvingmecfs.org/geo-gse327255/README.md (if mirrored) or see the catalogue
- `geo-gse93266`: https://data.solvingmecfs.org/geo-gse93266/README.md (if mirrored) or see the catalogue
- `sra-prjna940635`: https://data.solvingmecfs.org/sra-prjna940635/README.md (if mirrored) or see the catalogue
- `vdjdb`: https://data.solvingmecfs.org/vdjdb/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-singlecell-methylation/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-singlecell-methylation/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `INPUTS_NOTE.md` | 5.4 KB |
| `data/cells_all.parquet` | 53.7 MB |
| `data/pb_B_memory.parquet` | 1.4 MB |
| `data/pb_B_naive.parquet` | 1.5 MB |
| `data/pb_CD4_Tem.parquet` | 1.6 MB |
| `data/pb_CD4_naive_CM.parquet` | 2.0 MB |
| `data/pb_CD8_Tem.parquet` | 1.5 MB |
| `data/pb_CD8_Temra.parquet` | 1.5 MB |
| `data/pb_CD8_naive_CM.parquet` | 1.4 MB |
| `data/pb_MAIT.parquet` | 1.2 MB |
| `data/pb_Mono_classical.parquet` | 2.1 MB |
| `data/pb_Mono_intermediate.parquet` | 1.7 MB |
| `data/pb_Mono_nonclassical.parquet` | 1.5 MB |
| `data/pb_NK_CD56dim.parquet` | 1.7 MB |
| `data/pb_Platelet.parquet` | 842.0 KB |
| `data/pb_Treg.parquet` | 1.3 MB |
| `data/pb_cDC2.parquet` | 1.1 MB |
| `figures/fig1_scRNA_prespecified_baseline.png` | 121.0 KB |
| `figures/fig2_methylation_confound_clocks.png` | 119.5 KB |
| `results/bcr_group_tests.csv` | 15.5 KB |
| `results/bcr_sample_metrics.csv` | 69.6 KB |
| `results/bulk_replication_scRNA_leads.csv` | 8.8 KB |
| `results/meth_GSE111183_sample_scores.csv` | 20.2 KB |
| `results/meth_GSE156792_sample_scores.csv` | 71.3 KB |
| `results/meth_GSE327255_sample_scores.csv` | 46.2 KB |
| `results/meth_GSE93266_sample_scores.csv` | 56.1 KB |
| `results/meth_clocks_composition_group_tests.csv` | 18.3 KB |
| `results/meth_ewas_GSE111183_top2000.csv` | 476.4 KB |
| `results/meth_ewas_GSE156792_top2000.csv` | 463.0 KB |
| `results/meth_ewas_GSE327255_top2000.csv` | 473.1 KB |
| `results/meth_ewas_GSE93266_top2000.csv` | 461.5 KB |
| `results/meth_ewas_meta_top5000_modelB.csv` | 1.2 MB |
| `results/meth_ewas_meta_top5000_modelC.csv` | 1.2 MB |
| `results/meth_ewas_replication_top1000_modelB.csv` | 976 B |
| `results/meth_ewas_replication_top1000_modelC.csv` | 1010 B |
| `results/meth_targeted_cpgs_modelB.csv` | 2.3 MB |
| `results/meth_targeted_cpgs_modelC.csv` | 2.3 MB |
| `results/sc_celltype_counts_per_sample.csv` | 9.9 KB |
| `results/sc_composition_scores_tests.csv` | 30.0 KB |
| `results/sc_features_long.csv` | 1.2 MB |
| `results/sc_pctmt_by_celltype_personmean.csv` | 857 B |
| `results/sc_pseudobulk_DE_q10.csv` | 160 B |
| `results/sc_pseudobulk_DE_summary.csv` | 1.2 KB |
| `results/sc_pseudobulk_DE_top30_per_celltype.csv` | 95.9 KB |
| `results/sc_sample_meta.csv` | 7.5 KB |
| `results/sc_sample_quality.csv` | 16.3 KB |
| `results/sc_sample_quality_tests.csv` | 1.2 KB |
| `results/tcr_group_tests.csv` | 21.5 KB |
| `results/tcr_sample_metrics.csv` | 80.8 KB |
| `scripts/01_sc_annotate.py` | 6.3 KB |
| `scripts/02_sc_collect.py` | 2.8 KB |
| `scripts/03_sc_composition_scores.py` | 6.0 KB |
| `scripts/04_sc_pseudobulk_de.py` | 3.6 KB |
| `scripts/05_sc_sample_quality.py` | 2.3 KB |
| `scripts/06_bulk_replication.py` | 5.2 KB |
| `scripts/10_meth_prepare.py` | 972 B |
| `scripts/11_meth_clocks_composition.py` | 7.4 KB |
| `scripts/12_meth_ewas.py` | 5.7 KB |
| `scripts/12b_meth_replication.py` | 3.9 KB |
| `scripts/13_meth_targeted.py` | 2.9 KB |
| `scripts/20_bcr_process.py` | 3.8 KB |
| `scripts/21_tcr_process.py` | 3.2 KB |
| `scripts/22_bcr_tcr_analysis.py` | 5.1 KB |
| `scripts/30_figures.py` | 3.5 KB |
