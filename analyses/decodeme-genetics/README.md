# DecodeME summary statistics re-analysed (loci, MAGMA, colocalisation, genetic correlation)

Catalogue id: `analysis-decodeme-genetics`

**What we found (one line):** All 8 DecodeME loci reproduce; genetic signal is enriched in brain tissue and neurons and in none of 29 immune cell types; infection and non-infection onset are genetically the same (rg 0.98); strong overlap with chronic pain (rg 0.45) and depression (rg 0.55).

This is an exploratory re-analysis of public data by the Solving ME/CFS project. It has not been peer reviewed. It is research reading, not medical advice.

- **Citation:** Solving ME/CFS research notes (2026). DecodeME summary statistics re-analysed (loci, MAGMA, colocalisation, genetic correlation). https://data.solvingmecfs.org/analysis-decodeme-genetics/
- **Source:** https://data.solvingmecfs.org/analysis-decodeme-genetics/README.md
- **Sample size:** see README
- **Data type:** our analysis: GWAS follow-up (MAGMA, coloc, LDSC); scripts, result tables, figures
- **Licence:** CC-BY-4.0 (results, figures, notes); MIT (scripts) (https://creativecommons.org/licenses/by/4.0/)
- **Attribution:** Cite 'Solving ME/CFS data catalogue' with the URL; credit the underlying datasets listed in the README.
- **How the licence was checked:** Our own work.

## Notes

Exploratory re-analysis, not peer reviewed. eQTLGen region extract excluded (licence not stated). QC-filtered parquet copies of the sumstats are not mirrored (duplicates of decodeme-gwas); the GRCh38-to-rsID map is included.

## Licences for this bundle

Scripts (`.py`, `.R`, `.sh`) are released under the MIT licence; result tables, figures and notes under CC BY 4.0. Derived tables inherit attribution requirements from the datasets they were computed from (listed below).

## Datasets used

- `decodeme-gwas`: https://data.solvingmecfs.org/decodeme-gwas/README.md (if mirrored) or see the catalogue
- `decodeme-mvp-meta`: https://data.solvingmecfs.org/decodeme-mvp-meta/README.md (if mirrored) or see the catalogue

## Files

Base URL: `https://data.solvingmecfs.org/analysis-decodeme-genetics/` (fallback `https://pub-befaaf5d33ee4324b51dc65bd8a75aa2.r2.dev/analysis-decodeme-genetics/`). Files ending `.gz` that were plain text at source were gzipped by us; content is unchanged.

| path | size |
|---|---|
| `downloads/decodeme/INPUTS_NOTE.md` | 5.9 KB |
| `downloads/decodeme/magma_results/gwas_1.genes.out` | 1.5 MB |
| `downloads/decodeme/magma_results/gwas_1.gp.brain_sn.gsa.out` | 4.5 KB |
| `downloads/decodeme/magma_results/gwas_1.gp.gtex54.gsa.out` | 6.8 KB |
| `downloads/decodeme/magma_results/gwas_1.gp.hpa_sc.gsa.out` | 18.9 KB |
| `downloads/decodeme/magma_results/gwas_1.gp.monaco_immune.gsa.out` | 3.8 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c2.cp.kegg_medicus.gsa.out` | 108.2 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c2.cp.reactome.gsa.out` | 254.3 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c5.go.bp.gsa.out` | 1.1 MB |
| `downloads/decodeme/magma_results/gwas_1.gs.c5.go.cc.gsa.genes.out` | 1.7 MB |
| `downloads/decodeme/magma_results/gwas_1.gs.c5.go.cc.gsa.out` | 134.6 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c5.go.cc.gsa.sets.genes.out` | 1.9 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c5.go.mf.gsa.out` | 259.1 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c7.immunesigdb.gsa.out` | 722.9 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c8.all.gsa.genes.out` | 1.7 MB |
| `downloads/decodeme/magma_results/gwas_1.gs.c8.all.gsa.out` | 119.6 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.c8.all.gsa.sets.genes.out` | 1.2 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.h.all.gsa.out` | 6.6 KB |
| `downloads/decodeme/magma_results/gwas_1.gs.project_sets.gsa.out` | 1.0 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.genes.out` | 1.5 MB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gp.brain_sn.gsa.out` | 4.5 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gp.gtex54.gsa.out` | 6.8 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gp.hpa_sc.gsa.out` | 18.9 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gp.monaco_immune.gsa.out` | 3.8 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c2.cp.kegg_medicus.gsa.out` | 108.2 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c2.cp.reactome.gsa.genes.out` | 1.7 MB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c2.cp.reactome.gsa.out` | 254.3 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c2.cp.reactome.gsa.sets.genes.out` | 10.9 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c5.go.bp.gsa.out` | 1.1 MB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c5.go.cc.gsa.out` | 134.6 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c5.go.mf.gsa.out` | 259.1 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c7.immunesigdb.gsa.out` | 722.9 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.c8.all.gsa.out` | 119.6 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.h.all.gsa.out` | 6.6 KB |
| `downloads/decodeme/magma_results/gwas_1_infectious_onset.gs.project_sets.gsa.out` | 1.0 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.genes.out` | 1.5 MB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gp.brain_sn.gsa.out` | 4.5 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gp.gtex54.gsa.out` | 6.8 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gp.hpa_sc.gsa.out` | 18.9 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gp.monaco_immune.gsa.out` | 3.8 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c2.cp.kegg_medicus.gsa.out` | 108.2 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c2.cp.reactome.gsa.out` | 254.3 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c5.go.bp.gsa.out` | 1.1 MB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c5.go.cc.gsa.out` | 134.6 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c5.go.mf.gsa.out` | 259.1 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c7.immunesigdb.gsa.out` | 722.9 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.c8.all.gsa.out` | 119.6 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.h.all.gsa.out` | 6.6 KB |
| `downloads/decodeme/magma_results/gwas_1_non_infectious_onset.gs.project_sets.gsa.out` | 1.0 KB |
| `downloads/decodeme/magma_results/gwas_2.genes.out` | 1.5 MB |
| `downloads/decodeme/magma_results/gwas_2.gp.brain_sn.gsa.out` | 4.5 KB |
| `downloads/decodeme/magma_results/gwas_2.gp.gtex54.gsa.out` | 6.8 KB |
| `downloads/decodeme/magma_results/gwas_2.gp.hpa_sc.gsa.out` | 18.9 KB |
| `downloads/decodeme/magma_results/gwas_2.gp.monaco_immune.gsa.out` | 3.8 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c2.cp.kegg_medicus.gsa.out` | 108.2 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c2.cp.reactome.gsa.out` | 254.3 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c5.go.bp.gsa.out` | 1.1 MB |
| `downloads/decodeme/magma_results/gwas_2.gs.c5.go.cc.gsa.genes.out` | 1.7 MB |
| `downloads/decodeme/magma_results/gwas_2.gs.c5.go.cc.gsa.out` | 134.6 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c5.go.cc.gsa.sets.genes.out` | 1.9 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c5.go.mf.gsa.out` | 259.1 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c7.immunesigdb.gsa.out` | 722.9 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c8.all.gsa.genes.out` | 1.7 MB |
| `downloads/decodeme/magma_results/gwas_2.gs.c8.all.gsa.out` | 119.6 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.c8.all.gsa.sets.genes.out` | 1.2 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.h.all.gsa.out` | 6.6 KB |
| `downloads/decodeme/magma_results/gwas_2.gs.project_sets.gsa.out` | 1.0 KB |
| `downloads/decodeme/qc_parquet/id_rs_map.tsv.gz` | 89.1 MB |
| `downloads/decodeme/results/coloc_eqtlgen_blood.csv` | 12.6 KB |
| `downloads/decodeme/results/coloc_ld_check.csv` | 1.8 KB |
| `downloads/decodeme/results/gtex_v8_lead_qtl_lookup.csv` | 29.3 KB |
| `downloads/decodeme/results/gtex_v8_lead_qtl_summary.csv` | 3.3 KB |
| `downloads/decodeme/results/ldsc_rg.csv` | 1.2 KB |
| `downloads/decodeme/results/lead_loci_all_gwas.csv` | 4.7 KB |
| `downloads/decodeme/results/longcovid_lookup_decodeme_suggestive.csv` | 5.1 KB |
| `downloads/decodeme/results/magma_gwas_1_genes.csv` | 1.5 MB |
| `downloads/decodeme/results/magma_gwas_1_genesets.csv` | 3.0 MB |
| `downloads/decodeme/results/magma_gwas_1_infectious_onset_genes.csv` | 1.5 MB |
| `downloads/decodeme/results/magma_gwas_1_infectious_onset_genesets.csv` | 3.0 MB |
| `downloads/decodeme/results/magma_gwas_1_infectious_onset_project_candidate_genes.csv` | 10.4 KB |
| `downloads/decodeme/results/magma_gwas_1_infectious_onset_tissue_celltype.csv` | 30.8 KB |
| `downloads/decodeme/results/magma_gwas_1_non_infectious_onset_genes.csv` | 1.5 MB |
| `downloads/decodeme/results/magma_gwas_1_non_infectious_onset_genesets.csv` | 2.8 MB |
| `downloads/decodeme/results/magma_gwas_1_non_infectious_onset_project_candidate_genes.csv` | 11.7 KB |
| `downloads/decodeme/results/magma_gwas_1_non_infectious_onset_tissue_celltype.csv` | 30.7 KB |
| `downloads/decodeme/results/magma_gwas_1_project_candidate_genes.csv` | 11.6 KB |
| `downloads/decodeme/results/magma_gwas_1_tissue_celltype.csv` | 31.0 KB |
| `downloads/decodeme/results/magma_gwas_2_genes.csv` | 1.4 MB |
| `downloads/decodeme/results/magma_gwas_2_genesets.csv` | 3.0 MB |
| `downloads/decodeme/results/magma_gwas_2_project_candidate_genes.csv` | 11.2 KB |
| `downloads/decodeme/results/magma_gwas_2_tissue_celltype.csv` | 30.8 KB |
| `downloads/decodeme/results/magma_tissue_celltype.png` | 210.4 KB |
| `downloads/decodeme/results/manhattan_gwas1_infection.png` | 329.3 KB |
| `downloads/decodeme/results/meta_mvp_gw_loci.csv` | 311 B |
| `downloads/decodeme/results/meta_mvp_lookup.csv` | 578 B |
| `downloads/decodeme/results/paper_table3_lookup.csv` | 1.5 KB |
| `downloads/decodeme/results/project_sets_symbols.tsv` | 2.0 KB |
| `downloads/decodeme/scripts/01_load_qc.py` | 963 B |
| `downloads/decodeme/scripts/01b_load_qc_arrow.py` | 902 B |
| `downloads/decodeme/scripts/02_liftover_rsid.py` | 1.3 KB |
| `downloads/decodeme/scripts/03_gene_covariates.py` | 1.9 KB |
| `downloads/decodeme/scripts/04_magma.sh` | 1.4 KB |
| `downloads/decodeme/scripts/04b_geneprop.sh` | 371 B |
| `downloads/decodeme/scripts/05_leads.py` | 3.3 KB |
| `downloads/decodeme/scripts/06_summarise_magma.py` | 2.7 KB |
| `downloads/decodeme/scripts/07_coloc_eqtlgen.py` | 3.5 KB |
| `downloads/decodeme/scripts/07b_ld_check.py` | 1.4 KB |
| `downloads/decodeme/scripts/08_gtex_lookup.py` | 1.4 KB |
| `downloads/decodeme/scripts/09_cross_trait.py` | 4.9 KB |
| `downloads/decodeme/scripts/10_longcovid_lookup.py` | 2.6 KB |
| `downloads/decodeme/scripts/11_figures.py` | 1.7 KB |
| `downloads/decodeme/scripts/12_meta_lookup.py` | 1.3 KB |
| `downloads/decodeme/scripts/custom_sets.py` | 1.9 KB |
| `genetics/candidate_gene_gtex_expression.csv` | 5.5 KB |
| `genetics/decodeme_lead_variants.csv` | 437 B |
| `genetics/gene_integration.csv` | 3.5 KB |
| `genetics/locus_genes_500kb.csv` | 3.0 KB |
| `genetics/locus_tissue_permutation.csv` | 3.6 KB |
| `genetics/scripts/01_locus_genes.py` | 2.3 KB |
| `genetics/scripts/02_gtex_expression.py` | 2.0 KB |
| `genetics/scripts/03_locus_tissue_permutation.py` | 2.5 KB |
| `genetics/wanted_urls.txt` | 3.2 KB |
