# Atlas coverage against the reference literature (proposal, 2026-10-29)

Reference set: 1378 unique papers from Europe PMC, in four groups:
- (a) reviews: 60 most-cited ME/CFS reviews from 2015-26, 40 recent ones and 25 long COVID (LC) mechanism reviews
- (b) the 300 most-cited ME/CFS primary biology papers
- (c) 696 ME/CFS biology papers from 2024-26
- (d) 391 LC papers on PEM, exercise and fatigue mechanisms

Scripts are in `scripts/` (01-09) and the raw Europe PMC pulls are in `downloads/`.

## Coverage

| Group | n | In atlas (DOI or author-year) | Proposed here | Out of scope | Not matched (covered by a broader finding, or not screened) |
|---|---|---|---|---|---|
| All | 1378 | 235 (17%) | 145 (11%) | 52 (4%) | 946 (69%) |
| a: ME/CFS reviews, most-cited | 60 | 13 | 0 | 0 | 47 |
| a: ME/CFS reviews, recent | 40 | 10 | 1 | 1 | 28 |
| a: LC reviews | 25 | 0 | 1 | 1 | 23 |
| b: most-cited primary | 300 | 75 (25%) | 29 | 3 | 193 |
| c: 2024-26 primary | 696 | 132 (19%) | 68 | 43 | 453 |
| d: LC PEM/exercise | 391 | 32 (8%) | 63 | 7 | 289 |

By domain:
- Rows were assigned to a domain by the matched finding, otherwise by the top TF-IDF suggestion, which is weak.
- "Proposed" is the share of each row's papers that this proposal adds.

| Domain | n | In atlas | Proposed |
|---|---|---|---|
| muscle-exertion | 144 | 19% | 12% |
| vascular | 79 | 20% | 18% |
| neuro | 230 | 10% | 16% |
| immune | 200 | 18% | 11% |
| infection | 182 | 9% | 9% |
| omics | 94 | 18% | 14% |
| endocrine | 82 | 17% | 7% |
| gut | 74 | 19% | 7% |
| genetics | 65 | 14% | 8% |
| trials | 184 | 9% | 5% |
| key reviews | 44 | 100% | - |

How to read these numbers:
- The "In atlas" share is a lower bound. Many unmatched papers are small replications of findings the atlas already makes at theme level.
- Title screening suggests most of the 946 unmatched papers fall in that group.
- The atlas is strongest on its named flagship cohorts (NIH, Cornell, Columbia, Charité, Bergen, Amsterdam). It is weakest on the LC exercise physiology literature and on pre-2010 neurophysiology.

## Biggest thematic gaps

1. **LC exercise physiology.** Covers peak VO2, invasive CPET (iCPET) oxygen-extraction failure, fat oxidation and lactate, 31P-MRS and NIRS oxidative capacity, myopathy on biopsy and EMG, and the Tryfonos 2024 crossover RCT that provoked no PEM. In group d, 8% of papers are matched.
2. **LC brain imaging.** Covers blood-brain barrier disruption (Greene 2024), subcortical and choroid plexus changes, and contested FDG-PET results. The atlas has none of these.
3. **Pain modulation.** Exercise-induced hyperalgesia and impaired conditioned pain modulation, a Meeus/Nijs body of work replicated across several studies. The atlas has none of it.
4. **Electrophysiology and network imaging.** Covers qEEG, resting-state fMRI (contested), serotonergic PET and basal ganglia reward fMRI.
5. **Negative and null trials.** Covers clonidine (NorCAPITAL), fluoxetine, and the 2026 LDN RCT. Also missing are small positives such as solriamfetol, saline, SCIG and vitamin D.
6. **Large post-COVID incidence registries and infection specificity.** Registry data show that bacterial pneumonia also raises ME/CFS risk.
7. **Null virus and persistence studies.** Covers blood virome studies, CATCH-EVseq (no enterovirus) and the plasma-antigen nulls. These should sit beside the positive persistence findings.
8. **Autonomic testing in LC.** Covers objective autonomic testing in LC and MSNA results that disagree with each other, plus the classic tilt papers (Bou-Holaigah 1995, Stewart 1999), which are missing from neuro-pots-prevalence.

## 15 most important additions (ranked)

1. lc-icpet-impaired-extraction (new, several): six iCPET/NIRS studies. It also justifies raising vascular-oxygen-extraction to several.
2. muscle-exertion-lc-peak-vo2 (new, several): Durstenfeld 2022 meta-analysis, plus contradicting cohorts.
3. lc-pem-provocation-response (new, contested): the Tryfonos 2024 RCT null against Berardi 2026.
4. neuro-pain-modulation-exercise-hyperalgesia (new, several).
5. lc-myopathy-histology-emg (new, several): Hejbøl, Agergaard and Tryfonos, against the null MRI and jitter studies.
6. lc-muscle-oxidative-capacity (new, several): Finnigan 2024 31P-MRS.
7. Grade change: vascular-muscle-capillary-basement-membrane from early to several (Aschman 2023 and Tryfonos 2025).
8. neuro-lc-bbb-disruption (new, several).
9. trials-clonidine (new, refuted). Clonidine is still used clinically.
10. infection-blood-virome-null (new, early), plus null contradicting entries for SARS-CoV-2 antigen and enterovirus.
11. Grade change: neuro-injury-markers from several to contested (Azcue 2024, NfL higher).
12. infection-viral-vs-bacterial-registry: grade change to contested, plus the Hsu 2024 and Chang 2025 registries.
13. omics-lc-serotonin (new, contested): Wong 2023 Cell against Mathé 2025.
14. neuro-lc-autonomic-objective (new, several), plus the classic tilt papers added to neuro-pots-prevalence.
15. trials-ldn: add the Naik 2026 RCT (n=160, null), which confirms the refuted grade.

## Caveats

- **Author-year matching is noisy.** Common surnames gave false matches, so an author-year match was counted only when the atlas row had no DOI.
- **Substance coverage was judged from titles.** The 69% unmatched figure is approximate and is not a list of gaps.
- **Preprints are included.** They are labelled in `cite`.
- **Two cited entries are verified=false** because their abstracts were not read: de Boer 2022 (a letter) and Rahmati 2025.
- **Some add_studies results are generic.** For roughly 15 entries (EV, methylation, GPCR-AAb, gut, Saurel, Hamrefors and others), the `result` text is a placeholder. Tighten it from the abstract before merging.
- **Subdomains for trials-ssri and trials-solriamfetol are placeholders.** New subdomains are suggested in `detail`.
