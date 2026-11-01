# Cohorts table (kb/cohorts.*)

Built 2026-11-01. Files:

- `cohorts.json` is the master copy: an array of 400 cohort objects.
- `cohorts.csv` has the same rows. List fields are joined with ` | ` and `severity_counts_raw` is a JSON string.
- `linkcheck.csv` holds the source-URL checks.
- `schema.md` defines the fields (see its `## cohorts` section).

The project treats long COVID with post-exertional malaise (PEM) as essentially the same disease as ME/CFS. Long COVID populations are included and labelled in `population`. PEM was not verified for most long COVID entries.

## What is in it

| Curation level | Rows | Meaning |
|---|---|---|
| `curated` | 101 | Every population behind the 110 catalogue datasets (where one exists), every cohort-type row of the 88-row INVENTORY, and the main named cohorts behind atlas findings. Severity was looked for in each one: in the paper (PMC full text), in our per-person downloads, or both. |
| `atlas_stub` | 299 | Primary studies cited in `atlas/build/research/*.json` that did not map to a curated cohort. Fields come from the atlas `n`/`design` text and the DOI. These rows were **not reviewed for severity**, apart from 3 that have explicit overrides. |

We left out the following:

- Funders, roadmaps, networks and portals: dhsc-plan, dekade, mecfs-research-foundation, zonmw, euromene, cordis, meruk-mea, lcrf, arpah-leap-czi, icancme, polybio, plrc, nih-roadmap, nih-crc-mapmecfs, immme-nkse and mecfs-research-review.
- Population biobanks with no published ME/CFS case series: hunt, moba, dbds, lifegene, estonian-biobank, ofh, generation-scotland, opensafely, bbj, clsa, ckb, korea-biobank and clalit.
- Inventory rows with no ME/CFS cohort data found: scripps-wearable, convalescence and simmaron. prime-uk is attached to decodeme.
- Catalogue entries that are data services or analyses rather than study populations: pubmed-mecfs, clinicaltrials-gov, gwas-catalog, vdjdb, pqtl-ukbppp-decode, mapmecfs, github-rti-mapmecfs-api, github-mecfs-bioinformatics, the `analysis-*` READMEs except where linked to a cohort, and several code-only GitHub repos (github-tijn-thesis, bs5ga, jasemurphy, k-roy, pkozakowski, prmunn, lgchaplin) whose cohort could not be identified.
- Atlas reviews, meta-analyses, guidelines, models, animal-only work and our own re-analyses. The re-analyses are mapped through `cites` onto the cohorts they used.

`cites` lists every matching atlas cite string with its DOI. `aliases` maps the cohort ids used in `effects.csv` (for example `germain-2022`) onto these cohorts. These link atlas findings to a `cohort_id`: 472 cite strings in all, 167 of them on curated cohorts.

## How severity was filled (rules)

1. The `n_mild`, `n_moderate`, `n_severe`, `n_very_severe`, `n_housebound` and `n_bedbound` fields hold a number only when a source gives a category with exactly that label. Intermediate categories go only in `severity_counts_raw`, verbatim. Examples are Bergen "mild/moderate" and "moderate/severe", and Mount Sinai "mild to moderate".
2. `n_severe` means functional severity: NICE or ICC "severe", or a group the authors defined as house- or bed-bound. Symptom-based or case-definition-based "severe" goes in raw only. Examples:
   - Hurwitz 2009 "severe (n=30)" was split on symptom scores.
   - Northwestern IM "severe ME/CFS" means meeting more than one case definition.
3. Bell scores are not converted to classes. The Cornell distribution is given as computed Bell counts in raw, and the anchor wording is quoted in the note.
4. `severity_denominator` is the number of people the counts refer to. For several cohorts this is a subset (Oslo 45 of 426, Bergen RituxME 111 of 151, Stockholm 31 of 229).
5. `severity_counts_basis` is either `reported` (taken from a paper) or `computed_from_per_person_data` (computed by us). Every number has a quote and its location (table, PMC id, file) in `severity_note`. Nothing is estimated. If a paper gives only percentages, the percentages are quoted and the count stays null.

### Computed from our downloads (per-person data)

- **Cornell 2-day CPET**: Germain 2022 supplement `biomap/BioMapAI/data/indenpendent_cohorts/jci.insight.2022.xlsx`, 60 cases. Bell median 37.5 (IQR 27.5–40). Bell ≤30: n=29. Bell ≤20: n=15. Bell ≤10: n=1. 75% female.
- **BioMapAI**: `codes/AI/input/metadata.csv`, one row per person. Duration was Short (<4 y since diagnosis) for 75 cases and Long (>10 y) for 78. 112 of 153 cases were female. Median age was 44. The public severity scores are rescaled, so no classes were made from them.
- **Bergen CycloME/PAD000026**: `PAD000026_Log2Matrix.tsv`, 50 cases. SF-36 PF median 30 (IQR 15–45). PF ≤20: n=20. Median 2,790 steps a day. This is in `notes`, not in the class fields.
- **UK ME/CFS Biobank repertoires**: BCR (PRJEB83659) has 37 severe, 31 mild/moderate, 24 healthy and 31 MS samples. TCR (PRJNA940635) has 304 samples, about 37–38 severe per cell type. Persons cannot be linked.
- **UAB family cohort**: 22 affected (21 F) and 7 unaffected, in 16 families.
- **Montreal GSE327255**: duration median 12 y (IQR 7–21.75). 2 of the 54 cases had illness for less than 2 years.
- **Cardiozorg**: the 2024 table was summed by us, giving 227/377/136.

## Totals by severity class

Fifteen cohorts give at least one count in a standard class. Three more give counts only in raw form: Cornell (computed Bell), Northwestern IM (definition-based severe) and Hurwitz (symptom-based severe).

| Class | Naive sum over cohorts | After overlap control* |
|---|---|---|
| Mild | 5,210 | 5,199 |
| Moderate | 10,783 | 10,734 |
| Severe | 2,705 | 2,681 |
| Very severe | 210 | 210 |
| Housebound (explicit) | 8 | 8 |
| Bedbound / bedridden (explicit) | 7 | 7 |

\*Within each `overlap_group` (for example bergen, stockholm, cardiozorg) we took the largest count per class instead of summing. This assumes trial and sub-study patients from the same clinic may be the same people.

Caveats:

- **DecodeME dominates.** It accounts for 2,073 severe and 131 very severe, all self-reported against NICE categories with no clinical verification. Mount Sinai CoRE adds 356 severe and 78 very severe, also from a self-selected online survey of ME/CFS respondents. Its long COVID column is kept in raw only. Without these two surveys, clinically assessed cohorts give about **252 severe** people (276 before overlap control), which is small.
- **UK ME/CFS Biobank**: the 54 severe (house- or bed-bound, mostly sampled at home) are the largest clinically verified, biosampled severe group. The many UKMEB sub-studies (Dibble, Bourgonje, Lee JS, Ryback & Cowan, BCR/TCR) reuse these same people and are not added again.
- Overlap between datasets that share participants is handled only by `overlap_group`. Overlap across groups cannot be ruled out. Examples are UKMEB volunteers who also joined DecodeME, and Bergen patients also in Oslo genotyping. Not all categories are equivalent either: NICE self-report, the clinician 5-level scale, ICC grading and DSQ item 79 measure different things.
- Housebound and bedbound counts are tiny because few papers report them as counts. Oslo (Lande 2020) gives 12.5% bedridden and 28.6% strictly housebound of a 426-person genotyped cohort, but only as percentages.

### Cohorts richest in severe patients

1. DecodeME: 2,073 severe and 131 very severe (self-report, saliva DNA only).
2. Mount Sinai CoRE survey: 356 severe and 78 very severe ME/CFS (survey only, no samples).
3. Stichting CardioZorg: 136 severe by ICC (clinic records and tilt; bedbound excluded).
4. UK ME/CFS Biobank: 54 severe, house- or bed-bound, with blood, PBMC and DNA, available by application.
5. OMF SIPS: 20 severe, all severely ill by design.
6. Bergen RituxME: 19 severe of 111. Bergen Sandvik 2023: 9 of 39. Stockholm/Bragée: 7 of 31 (CSF) and 9 of 34. Hsieh 2023 (Epsom): 9.
7. Oslo biobank: 6 bedridden and 8 strictly housebound in a 45-person subset. About 12.5% bedridden across the 426.

## Flags

| Flag | Meaning | Counts |
|---|---|---|
| `includes_severe` | true if the cohort includes house- or bed-bound or NICE-severe people (by count or by documented design, such as home visits); false if the design excludes them (for example two maximal CPETs) | true 19, false 5, null for the rest |
| `includes_early_onset_lt2y` | true if people within 2 years of onset are documented: DecodeME 0.5–1 y (n=254), post-infection prospective cohorts, RECOVER, post-COVID cohorts sampled at 6 months, and Montreal (2 people) | 12 true |
| `prospective_pre_illness` | true if samples or data exist from before onset for at least some cases: Northwestern IM, UK Biobank incident cases, RECOVER uninfected enrolees, birth and population cohorts, EHRs | 12 true |

## Sources used

- `site/data/catalogue.json`: 110 datasets. `dataset_ids` match catalogue ids, and each id was checked to exist.
- `global-resources/INVENTORY.csv` and `INVENTORY.md`. `inventory_ids` match inventory ids.
- `atlas/build/research/*.json`: 496 distinct studies, both findings[].studies and comprehensive_studies.
- Local data:
  - DecodeME preprint PDF (Table 1)
  - Germain 2022 xlsx
  - BioMapAI metadata
  - PAD000026 matrix
  - deep2/singlecell_methyl results (BCR/TCR and GSE327255 metadata)
  - UAB repo metadata
  - datasets/downloads/deep2 SUMMARY and INPUTS_NOTE files
- PubMed E-utilities and PMC: DOIs were resolved with the PMC ID converter and esearch (425 of 460 DOIs mapped to a PMID). 377 PMC full-text XMLs were fetched with efetch db=pmc and searched for severity tables. Europe PMC was not needed.

## Known gaps

- **Severity tables not reachable.** The following tables were in PDFs or images, or the full text was not in PMC: Montoya 2017 (PNAS, abstract only), Hoel 2021 Table 1, RituxME (Ann Intern Med; only the follow-up subset was used), Fluge 2016, Bragée 2020, Kedor 2022 and Legler 2023 (only Bell summaries), MECFS-R (only "median Bell 30"), Walitt 2024 (inclusion wording only), and OMF SIPS Chang 2021 (the PMC file was empty, so the number came from Hung 2025).
- **Stub rows were not reviewed.** Their n was parsed from atlas text only where it was unambiguous ("N vs M"), and `population` for many is "ME/CFS (assumed from atlas context)". Check the paper before using them.
- **Recruitment years, case definitions and settings** are filled only where the source text was seen. Fields marked "not re-checked" came from catalogue or inventory descriptions.
- **NCNED registry**: its 5-level table (n=519) appears to grade PEM, not overall illness, so it was not used.
- **Cohort naming**: the Cornell CRC row merges several papers with different n. The cohort n is the largest single paper (Keller 2024: 84/71).
- **Not researched**: Simmaron, Bateman Horne Center biobank beyond BioMapAI, and Solve M.E. You+ME severity.

## Link check (linkcheck.csv)

158 unique URLs were checked: every catalogue `source_url` plus every URL and domain named in INVENTORY.csv. We fetched each with `curl -sL`, looked for the expected organisation, study or accession name, and searched the visible text for gambling or spam terms. JS-rendered pages were confirmed through their APIs (PRIDE, ENA, OSF, MetaboLights, figshare).

| Status | Count |
|---|---|
| ok | 109 |
| blocked | 41 |
| dead | 4 |
| hijacked | 2 |
| unverified | 2 |

- **Hijacked**:
  - `decodeme.org.uk`, title "Non GamStop Casinos UK Rankings…"
  - `youandmeregistry.com`, which redirects to purpose4mypain.org, title "HOKI178: Link Situs Judi Slot Online…"
- **Dead**:
  - `solvecfs.org/you-me-registry/`, now a redirect to a PNG image. This is the catalogue source_url for youandme-registry.
  - `clalitresearch.org` (connection reset)
  - `ncnp.go.jp` and `mecfs-research-review.org` (proxy CONNECT 502)
- **Blocked**: 27 of these are GitHub, which returns 403 to this sandbox and says nothing about the repositories themselves. The rest are Cloudflare or Access Denied pages: MDPI, PNAS, Wiley, cureme.lshtm.ac.uk, omf.ngo, cprd.com, cdc.gov, ukbiobank.ac.uk, allofus and polybio.
- **Unverified**: two Elsevier DOIs resolve to linkinghub JavaScript redirects.

The earlier false positives (icancme.ca "slot" in CSS, lifelines.nl "judi" in a person's name, a Nature page) were checked and are ok.
