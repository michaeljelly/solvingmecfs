# Inputs obtained by the epi_literature stream (2026-09-24)

## 1. Risteys (FinnGen R13/R14 + FinRegistry), endpoint G6_POSTVIRFAT and 9 comparison endpoints
- URLs: https://risteys.finngen.fi/endpoints/G6_POSTVIRFAT (HTML + Phoenix LiveView websocket for the Relationships table); legacy API https://r10.risteys.finngen.fi/api/endpoints/G6_POSTVIRFAT/assocs.json and https://r9.risteys.finngen.fi/api/endpoints/<EP>/assocs.json
- Paths: risteys/<EP>_codewas.csv, <EP>_labwas.csv, <EP>_keyfigures.json, <EP>_assocs_r9/r10.json; risteys/G6_POSTVIRFAT.relationships.csv; flattened risteys/risteys_survival_assocs_r9_r10.csv
- Format CSV/JSON. G6_POSTVIRFAT_codewas 372 rows x 7; labwas 468 x 13; relationships 5014 x 8; survival table 9159 x 11.
- Endpoints: G6_POSTVIRFAT, G6_AUTONOM, M13_HYPERMOBILITY, AB1_EBV, R18_MALAI_FATIG, M13_FIBROMYALGIA, I9_HYPOTE, R18_SYNCO_COLLAP, D3_IMMUDEFPOSTEBV, AB1_CMV
- Sample rows (copied data):
```
code,vocab,description,OR,mlog10p,n_cases,n_controls
G93.3,ICD-10 Finland,Postviral fatigue syndrome,+∞,307.7,348,8
R53,ICD-10 Finland,Malaise and fatigue,10.8,108.0,245,627
concept,n_cases,n_controls,OR_measured,mlog10p_measured,...,mean_cases,mean_controls,unit,mlog10p_value,n_cases_v,n_controls_v
( 3036481 ) Cortisol [Moles/volume] in Serum or Plasma,204,295,15.30,152.93,...,386.22,376.44,nmol/l,0.29,195,278
{"FR_n": ["1622", "1089", "525"], "FG_n": ["334", "266", "68"], "FR_median_age": ["34.76", "35.80", "32.61"]}
```
- Surprising: the POSTVIRFAT survival HRs (before/after) are not computed in the current release (HR_FR column is all "-"); only r10 has them, and only for a subset. Relationships "cases_FR" is also "-" for FinRegistry; FG overlap counts ("N (Jaccard%)") are given. LabWAS mixes test-ordering OR (huge for cortisol, Borrelia) with value comparison (usually null).

## 2. FinnGen R12 manifest (case counts for O/E)
- URL: https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv (kept in /tmp/work_epi, 774 kB); used N=500,348.
- Derived: risteys/POSTVIRFAT_overlap_enrichment.csv (5014 rows: O, E, O/E, Poisson p, BH q) and risteys/POSTVIRFAT_selected_FR_vs_FG_enrichment.csv (37 rows; FinRegistry O/E using FR N=7,166,416)
```
endpoint,n_FR,longname,o_FR,o_FG,num_cases,OE_FG,q_FG,E_FR,OE_FR,p_FR
AB1_HERPES_SIMPLEX,37387.0,Herpesviral [herpes simplex] infections,31.0,13.0,4677.0,4.16,0.000146,8.46,3.66,1.96e-09
```
- Caveat: O/E is crude (no age/sex matching), and FR counts for the overlap came from the Risteys page, not individual data.

## 3. ClinicalTrials.gov API v2
- URL pattern: https://clinicaltrials.gov/api/v2/studies?query.cond=<cond>&filter.advanced=AREA[StudyType]INTERVENTIONAL&pageSize=100 (conds: "myalgic encephalomyelitis OR chronic fatigue syndrome", "long COVID OR post-COVID OR PASC", "postural orthostatic tachycardia")
- Raw: ctgov_raw_mecfs.json.gz (ME/CFS only, 1.1 MB gz); long COVID (17 MB) and POTS (5.9 MB) left in /tmp/work_epi/ctgov/.
- Derived: ctgov_all_interventional_trials.csv (1369 x 14; 282 ME/CFS, 915 LC, 308 POTS; 165 with results), ctgov_posted_results_outcomes.csv (1789 x 19), ctgov_curated_effects.csv (598 x 19, Hedges g from posted mean/SD), ctgov_primary_summary.csv
```
nct,title,n,otype,outcome,...,v1,v2,n1,n2,diff_arm1_minus_arm2,hedges_g,posted_p
NCT01071044,Efficacy and Safety of Lisdexamfetamine ...,26,PRIMARY,Change in BRIEF-A,...,21.38,3.36,15.0,11.0,18.02,1.34,.005
```
- Surprising: few ME/CFS trials post results (25/282); Hedges g from "change" arms uses SD of change, so g is inflated when SDs are small.

## 4. PubMed E-utilities
- esearch/efetch (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/). Query for main set: "chronic fatigue syndrome"[MeSH] OR myalgic encephalomyelitis OR postviral fatigue, all years.
- literature/raw/pubmed_mecfs.json.gz: 9404 records (7438 with abstracts), fields pmid, year, first, last, aff, journal, title, abstract, ptypes, mesh
- Extra sets: pm_lc_orth (538, long COVID orthostatic), pm_lc_rct (272), pm_recovery (2034, post-infective recovery/prognosis), pm_registry (370, registry/cohort)
```
{'pmid': '42773889', 'year': '2026', 'first': 'Ruseva', 'last': 'Martin', 'aff': 'QIMR Berghofer, Brisbane...', 'journal': 'Twin Res Hum Genet', 'title': 'Can Depression Polygenic Risk Score Predict Postinfective Fatigue Syndrome (PIFS...'}
```
- Derived: literature/biomarker_mentions_all.csv (854 sentence hits), biomarker_paper_calls.csv (659 paper-level calls, 439 papers), biomarker_replication_catalogue.csv (107 concepts). Regex precision on manual audit about 55-65%; recall is low (abstract-only).
- Curated by hand: curated/orthostatic_evidence.csv, curated/natural_history.csv, curated/registry_studies.csv

## Not obtained
- FinRegistry / FinnGen individual-level data: Findata permit (findata.fi) or FinnGen sandbox membership. Would allow age/sex-adjusted HRs of infections and drugs before G93.3 and a red-cell/reticulocyte test in POSTVIRFAT cases vs matched controls with proper covariates.
- Danish (DNPR via Danish Health Data Authority), Swedish (Socialstyrelsen patient register), UK Biobank (application, ~GBP 9k) and CPRD (ISAC protocol) for replication of the viral-infection and hypermobility signals.
