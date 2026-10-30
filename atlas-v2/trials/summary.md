# ME/CFS and long COVID treatment trials: database summary

Built 2026-10-29. Data file: `trials.json` (1028 trials). Scripts are in `scripts/`; raw downloads are in `raw/`.

Sources: ClinicalTrials.gov API v2 (all interventional ME/CFS records, plus long COVID records screened by target symptom), Europe PMC (published trials, matched to registry records), ISRCTN and EU CTIS. `verified: true` means the result was read in a published abstract or in ClinicalTrials.gov posted results. No result was taken from memory or from press releases.

## Counts

### By population
| Population | Trials |
|---|---:|
| long COVID | 618 |
| ME/CFS | 370 |
| both | 40 |

### By design
| Design | Trials |
|---|---:|
| RCT-placebo | 304 |
| single-arm | 232 |
| RCT-usual-care | 195 |
| RCT-active | 172 |
| crossover | 59 |
| open-label | 37 |
| RCT-unspecified | 21 |
| n-of-1 | 5 |
| unclear | 3 |

"RCT-unspecified" covers ISRCTN or CTIS records that are randomised but do not name the comparator.

### By status
| Status | Trials |
|---|---:|
| published only (not on CT.gov, ISRCTN or CTIS) | 441 |
| COMPLETED | 271 |
| UNKNOWN | 97 |
| RECRUITING | 76 |
| ACTIVE_NOT_RECRUITING | 30 |
| NOT_YET_RECRUITING | 29 |
| TERMINATED | 24 |
| ISRCTN (status not given) | 22 |
| WITHDRAWN | 22 |
| EU CTIS (status code not mapped) | 6 |
| ENROLLING_BY_INVITATION | 6 |
| SUSPENDED | 4 |

### By result
| Result | Trials |
|---|---:|
| positive | 332 |
| unreported-after-completion | 255 |
| pending | 186 |
| null | 100 |
| mixed | 86 |
| feasibility-only | 27 |
| withdrawn | 22 |
| not-assessable | 18 |
| negative/harm | 2 |

Two values were added to the requested set. **feasibility-only** means the trial reported only recruitment, adherence or safety. **not-assessable** means a paper or posted result exists but efficacy cannot be read from it, for example because there is no abstract or the trial was terminated with unusable data. **withdrawn** marks trials withdrawn before enrolment.

### By intervention class
| Class | Trials |
|---|---:|
| behavioural/rehab | 396 |
| other | 216 |
| supplement | 120 |
| device/neuromodulation | 88 |
| immunomodulator | 80 |
| metabolic/mitochondrial | 53 |
| autonomic/vasoactive | 28 |
| antiviral | 22 |
| microbiome | 13 |
| apheresis/immunoadsorption | 12 |

### By model step (one trial can count under several steps)

| Step | Trials |
|---|---:|
| 1 immune priming | 93 |
| 2 brain/endocrine control | 249 |
| 3 low blood volume | 11 |
| 4 vessels/capillaries | 32 |
| 5 tissue under-perfusion | 84 |
| 6 tissue injury/immune response | 139 |
| 7 sickness response | 91 |
| 8 sympathetic arousal/poor sleep | 70 |
| no step assigned (for example exercise, vitamins, acupuncture) | 477 |

### Verified results by population

| Population | positive | null | mixed | negative/harm | feasibility-only |
|---|---:|---:|---:|---:|---:|
| ME/CFS | 169 | 48 | 44 | 2 | 9 |
| long COVID | 156 | 51 | 39 | 0 | 18 |
| both | 7 | 1 | 3 | 0 | 0 |

### How often trials are positive, by design (verified efficacy results only)

| Group | Trials | % positive | % null |
|---|---:|---:|---:|
| All verified efficacy results | 520 | 64 | 19 |
| Uncontrolled (single-arm, open-label, n-of-1) | 158 | 87 | 2 |
| Randomised against usual care or wait list | 99 | 73 | 12 |
| Blinded placebo or sham (parallel or crossover) | 156 | 35 | 42 |
| Blinded placebo or sham with n of 100 or more | 36 | 17 | 56 |
| TCM, herbal, acupuncture or moxibustion | 56 | 80 | 5 |

The more tightly a trial is controlled, the less often it is positive. This is the strongest pattern in the database.

## Exclusions

ClinicalTrials.gov records screened out (ME/CFS and long COVID searches combined, after deduplication):

| Reason | Trials |
|---|---:|
| not an ME/CFS or long COVID population | 316 |
| primary outcome lung/breathlessness only | 68 |
| not a treatment trial (mechanistic, diagnostic, observational) | 61 |
| primary outcome psychiatric or sleep only | 31 |
| acute COVID, prevention or post-discharge recovery | 23 |
| primary outcome smell/taste only | 23 |
| primary outcome other organ (heart, skin, GI, eye, etc.) | 22 |
| primary outcome pain only | 18 |
| vaccine trial | 5 |
| dropped on manual review | 1 |

The long COVID screen kept trials aimed at fatigue, PEM, exercise tolerance, cognition, autonomic or orthostatic symptoms, or global function. It excluded trials of lung function or breathlessness only, smell or taste only, vaccines, acute COVID, and mental health only. Further exclusions after manual review: 3 diagnostic or mechanistic studies, 1 CT.gov record, and 56 Europe PMC publications that were not treatment trials or were out of scope (52 on abstract review plus 4 earlier). Also 6 duplicate publications or registry entries were merged into existing records, and the 2 platform parent records (RECOVER-AUTONOMIC, RECOVER-ENERGIZE) were split into their arms.

- ISRCTN: excluded (observational, diagnostic, mechanistic, breathlessness/mental-health only, or off-topic) 28; added as new record 22; merged into existing record 21
- CTIS: added as new record 6; merged into existing record 10; excluded (imaging, dyspnoea, acute COVID, prevention, off-topic) 12

## Overdue: completed more than 12 months ago with no results

255 trials have a primary completion date before 2025-10-29, have no results posted on ClinicalTrials.gov, and have no results paper that we could find in Europe PMC. They are 55 ME/CFS, 187 long COVID and 13 both. The 79 with UNKNOWN status are registry records that have not been updated past their completion date. Most of these probably finished, but some may never have run.

### Blinded or randomised trials, n of 50 or more (actual or planned), completed or terminated

| Trial | Intervention | Population | n | Primary completion | Sponsor |
|---|---|---|---:|---|---|
| NCT04950725 | breathing/respiratory training | long COVID | 2268 | 2022-12-31 | Mayo Clinic |
| NCT05642923 | vitamins/minerals | both | 528 | 2024-06-01 | ClinAmygate |
| NCT05513560 | pentoxifylline | long COVID | 460 | 2025-03-31 | University Health Network, Toronto |
| NCT05638633 | corticosteroid + vitamins/minerals | both | 321 | 2025-01-30 | Wuerzburg University Hospital |
| NCT05703074 | NADH/NAD+ precursor | long COVID | 310 | 2024-09-23 | University Hospital, Akershus |
| NCT00920777 | CBT/psychological therapy | ME/CFS | 234 | 2014-07 | Norwegian University of Science and Tech |
| NCT05823896 | nirmatrelvir/ritonavir | long COVID | 219 | 2024-11-28 | Karolinska Institutet |
| NCT04632719 | exercise/rehabilitation | long COVID | 200 | 2020-12-29 | University of Sao Paulo |
| NCT06568393 | exercise/rehabilitation | long COVID | 176 | 2024-04-01 | Bahçeşehir University |
| NCT02541890 | exercise/rehabilitation | ME/CFS | 175 | 2017-07 | Norwegian University of Science and Tech |
| NCT05846126 | immersive VR | long COVID | 172 | 2024-12-20 | Consorci Sanitari de Terrassa |
| NCT05689827 | meldonium/EMHPS | long COVID | 160 | 2022-11-10 | Promomed, LLC |
| NCT05356936 | vitamins/minerals | long COVID | 151 | 2024-04-01 | University Hospitals Cleveland Medical C |
| NCT05795816 | Testofen; Microcrystalline cellulose | long COVID | 150 | 2025-01-30 | RDC Clinical Pty Ltd |
| NCT05630339 | magnesium + vitamins/minerals | long COVID | 150 | 2023-05-01 | Coordinación de Investigación en Salud,  |
| NCT01765725 | Patient education program | ME/CFS | 146 | 2016-12 | Oslo University Hospital |
| NCT06063031 | exercise/rehabilitation + diet | long COVID | 146 | 2025-01-15 | Instituto de Investigación Sanitaria Ara |
| NCT06016192 | exercise/rehabilitation | long COVID | 139 | 2022-11-30 | University of Witten/Herdecke |
| NCT05298878 | exercise/rehabilitation | long COVID | 132 | 2023-09-22 | McGill University Health Centre/Research |
| NCT05848518 | exercise/rehabilitation | long COVID | 120 | 2023-09-30 | Campus docent Sant Joan de Déu-Universit |
| NCT05956405 | CBT/psychological therapy | long COVID | 120 | 2024-03-31 | Hospital Miguel Servet |
| NCT07184398 | cognitive rehabilitation | long COVID | 120 | 2024-09-30 | Universidad Libre seccional Cali |
| NCT07553897 | exercise/rehabilitation | long COVID | 120 | 2023-11-08 | Tri-Service General Hospital (TSGH) |
| NCT05911009 | BC007 (rovunaptabin) | long COVID | 119 | 2024-06-26 | Berlin Cures GmbH |
| NCT04944121 | RSLV-132 | long COVID | 112 | 2023-10-31 | Resolve Therapeutics |
| NCT04740736 | IV saline/fluids | ME/CFS | 106 | 2025-07-11 | Icahn School of Medicine at Mount Sinai |
| NCT04795557 | exercise/rehabilitation | long COVID | 100 | 2021-11-11 | Swedish Herbal Institute AB |
| NCT05676047 | cognitive rehabilitation | long COVID | 100 | 2025-01-21 | McMaster University |
| NCT05986422 | corticosteroid | long COVID | 96 | 2025-03-14 | Charite University, Berlin, Germany |
| NCT06124625 | exercise/rehabilitation | long COVID | 90 | 2025-09-08 | Schön Klinik Berchtesgadener Land |
| NCT07722000 | herbal/TCM | long COVID | 90 | 2023-12-01 | Beni-Suef University |
| NCT05798221 | herbal/TCM | long COVID | 86 | 2025-03-25 | Universität Duisburg-Essen |
| NCT05374668 | exercise/rehabilitation | long COVID | 85 | 2021-09-22 | Istinye University |
| NCT04718506 | exercise/rehabilitation + breathing/respiratory tr | long COVID | 83 | 2022-02-06 | Universidad de Murcia |
| NCT03691987 | microbiome therapy | ME/CFS | 80 | 2023-03-01 | University Hospital of North Norway |
| NCT05532904 | exercise/rehabilitation | long COVID | 80 | 2025-06-13 | Assistance Publique - Hôpitaux de Paris |
| NCT06492551 | exercise/rehabilitation | long COVID | 80 | 2022-09-23 | Trak Health Solutions S.L. |
| NCT07397130 | biofeedback | long COVID | 80 | 2025-08-01 | RSUP Persahabatan |
| NCT04983394 | exercise/rehabilitation | long COVID | 79 | 2022-02-27 | Cumhuriyet University |
| NCT01742013 | GCJBP Laennec Inj. | ME/CFS | 78 | 2013-08 | Ho Cheol Shin, M.D., Ph.D. |
| NCT05770505 | exercise/rehabilitation | long COVID | 77 | 2024-05-31 | Norwegian University of Science and Tech |
| NCT05435456 | exercise/rehabilitation | long COVID | 76 | 2021-03-31 | Halic University |
| NCT05901337 | acupuncture/manual therapy | long COVID | 76 | 2022-06-21 | Cairo University |
| NCT05216549 | exercise/rehabilitation | long COVID | 75 | 2022-06-05 | Józef Piłsudski University of Physical E |
| NCT04741841 | microbiome therapy | ME/CFS | 74 | 2021-03-30 | ImmuneBiotech Medical Sweden AB |
| NCT06492590 | exercise/rehabilitation | long COVID | 74 | 2025-09-30 | Consorci Sanitari de Terrassa |
| NCT01534130 | acupuncture/manual therapy | ME/CFS | 72 | 2012-08 | Chengdu University of Traditional Chines |
| NCT06840873 | exercise/rehabilitation | long COVID | 68 | 2024-06-30 | Shang-Lin Chiang |
| NCT04924881 | herbal/TCM | long COVID | 68 | 2022-10-03 | Chinese University of Hong Kong |
| NCT04301609 | sleep treatment | ME/CFS | 67 | 2022-12-30 | Vitae Health Innovation |
| NCT05710770 | immunoadsorption/apheresis | both | 66 | 2025-08-04 | Charite University, Berlin, Germany |
| NCT05204511 | exercise/rehabilitation | long COVID | 66 | 2023-05-31 | University of Vienna |
| NCT05254301 | CBT/psychological therapy + exercise/rehabilitatio | long COVID | 65 | 2024-01-30 | Universitair Ziekenhuis Brussel |
| NCT05954325 | immunoadsorption/apheresis | both | 63 | 2025-07-31 | Hannover Medical School |
| NCT05516836 | exercise/rehabilitation | long COVID | 62 | 2023-12-20 | Universidad Europea de Madrid |
| NCT06492577 | exercise/rehabilitation | long COVID | 61 | 2024-10-06 | Spitalul Clinic de Boli Infecțioase și P |
| NCT05323331 | exercise/rehabilitation | long COVID | 60 | 2022-10-01 | Riphah International University |
| NCT06315894 | exercise/rehabilitation | long COVID | 60 | 2024-09-05 | Fundación FLS de Lucha Contra el Sida, l |
| NCT06316843 | valacyclovir | long COVID | 59 | 2024-10-31 | Bateman Horne Center |
| NCT05631171 | CBT/psychological therapy | long COVID | 58 | 2023-06-30 | Adhera Health, Inc. |
| NCT02156128 | exercise/rehabilitation | ME/CFS | 54 | 2015-03 | Norwegian University of Science and Tech |
| NCT05477199 | CBT/psychological therapy | long COVID | 54 | 2022-09-15 | Universidad de Granada |
| NCT06814379 | exercise/rehabilitation + immersive VR | long COVID | 54 | 2025-07-31 | Universidad de Granada |
| NCT06419712 | vitamins/minerals | long COVID | 54 | 2023-04-20 | Escuela Militar de Graduados de Sanidad, |
| NCT05658536 | CBT/psychological therapy | long COVID | 51 | 2024-01-30 | University of Washington |
| NCT04794036 | exercise/rehabilitation | long COVID | 50 | 2022-06-13 | Universidad San Jorge |
| NCT05739032 | exercise/rehabilitation | long COVID | 50 | 2022-05-30 | Ukbe Sirayder |

### All other overdue trials

| Trial | Intervention | Status | n | Primary completion |
|---|---|---|---:|---|
| NCT01108549 | D-ribose | COMPLETED | 257 | 2009-09 |
| NCT01507701 | clonidine | COMPLETED | 5 | 2010-03 |
| NCT00468013 | exercise/rehabilitation | COMPLETED | 60 | 2010-03 |
| NCT00860236 | CBT/psychological therapy | COMPLETED | 40 | 2010-06 |
| NCT00540254 | CBT/psychological therapy | COMPLETED | 24 | 2010-07 |
| NCT01873482 | Movement to rhythm | COMPLETED | 17 | 2014-05 |
| NCT02082730 | sleep treatment + CBT/psychological therapy | COMPLETED | 20 | 2014-08-01 |
| NCT02336126 | CBT/psychological therapy | COMPLETED | 8 | 2015-01 |
| NCT02213679 | creatine/GAA | COMPLETED | 20 | 2015-06 |
| NCT01806246 | exercise/rehabilitation | COMPLETED | 13 | 2015-09 |
| NCT03562325 | CBT/psychological therapy | COMPLETED | 40 | 2017-10-31 |
| NCT02955420 | BC007 (rovunaptabin) | COMPLETED | 74 | 2018-04 |
| NCT05128292 | coenzyme Q10 | COMPLETED | 42 | 2018-09-01 |
| NCT04174300 | acupuncture/manual therapy | COMPLETED | 40 | 2020-03-10 |
| NCT04381793 | Recovery Factors | COMPLETED | 43 | 2020-07-04 |
| NCT04445376 | exercise/rehabilitation | COMPLETED | 20 | 2020-09-10 |
| NCT04435002 | herbal/TCM | COMPLETED | 39 | 2020-12-09 |
| NCT05013606 | molecular hydrogen | COMPLETED | 23 | 2021-04-01 |
| NCT04972864 | exercise/rehabilitation | COMPLETED | 24 | 2021-07-18 |
| NCT04980573 | Aromatherapy | COMPLETED | 47 | 2021-11-30 |
| NCT04841759 | exercise/rehabilitation | COMPLETED | 46 | 2021-12-03 |
| NCT04828668 | Targeted Wellness Formula C™ Sublingual Drops | COMPLETED | 32 | 2022-01-31 |
| NCT05311852 | polyphenol/omega-3/anti-inflammatory suppleme | COMPLETED | 34 | 2022-02-15 |
| NCT05571852 | cognitive rehabilitation + exercise/rehabilit | COMPLETED | 73 | 2022-02-28 |
| NCT05074771 | exercise/rehabilitation | COMPLETED | 31 | 2022-03-30 |
| NCT05107440 | CBT/psychological therapy | COMPLETED | 27 | 2022-06-20 |
| NCT05252481 | transcranial stimulation/neurofeedback | COMPLETED | 47 | 2022-06-26 |
| NCT05453188 | exercise/rehabilitation | COMPLETED | 30 | 2022-07-01 |
| NCT05973136 | exercise/rehabilitation | COMPLETED | 7 | 2022-08-08 |
| NCT04854772 | Mind Body Syndrome Therapy for Long Covid | COMPLETED | 23 | 2022-08-09 |
| NCT06161740 | exercise/rehabilitation | COMPLETED | 17 | 2022-09-20 |
| NCT04996212 | exercise/rehabilitation | COMPLETED | 70 | 2022-10-14 |
| NCT05119634 | exercise/rehabilitation | COMPLETED | 42 | 2022-10-21 |
| NCT07544186 | L-arginine/citrulline | COMPLETED | 43 | 2022-11-11 |
| NCT05543551 | breathing/respiratory training + exercise/reh | COMPLETED | 42 | 2022-12-01 |
| NCT05939622 | meldonium/EMHPS | COMPLETED | 30 | 2022-12-20 |
| NCT05302973 | exercise/rehabilitation | COMPLETED | 30 | 2022-12-31 |
| NCT05679505 | vagus nerve stimulation | COMPLETED | 40 | 2023-01-01 |
| NCT05225220 | vagus nerve stimulation | COMPLETED | 25 | 2023-01-06 |
| NCT04997395 | cannabinoid | COMPLETED | 12 | 2023-01-06 |
| NCT04604704 | low-dose naltrexone + NADH/NAD+ precursor | COMPLETED | 36 | 2023-01-23 |
| NCT05630378 | outpatient clinic with multimodal integrative | COMPLETED | 42 | 2023-02-01 |
| NCT05212831 | Inogen One® G4 | COMPLETED | 21 | 2023-03-02 |
| NCT05764070 | vagus nerve stimulation | COMPLETED | 44 | 2023-03-09 |
| NCT05172102 | breathing/respiratory training | COMPLETED | 40 | 2023-03-20 |
| NCT05566379 | CBT/psychological therapy | COMPLETED | 20 | 2023-07-29 |
| NCT05739552 | exercise/rehabilitation | COMPLETED | 30 | 2023-07-31 |
| NCT05228899 | IV saline/fluids | COMPLETED | 18 | 2023-10-01 |
| NCT05793736 | biofeedback | COMPLETED | 17 | 2023-10-30 |
| NCT05422924 | CBT/psychological therapy + diet | COMPLETED | 45 | 2023-10-31 |
| NCT05218174 | exercise/rehabilitation | COMPLETED | 21 | 2023-11-17 |
| NCT05338749 | cognitive rehabilitation | COMPLETED | 16 | 2023-11-30 |
| NCT07651371 | exercise/rehabilitation | COMPLETED | 74 | 2023-12-15 |
| NCT06437210 | Viusid Oral Solution | COMPLETED | 200 | 2023-12-15 |
| NCT06379737 | exercise/rehabilitation | COMPLETED | 72 | 2023-12-20 |
| NCT05664711 | stellate ganglion/sympathetic block | COMPLETED | 10 | 2023-12-21 |
| NCT05697042 | CBT/psychological therapy + exercise/rehabili | COMPLETED | 16 | 2023-12-30 |
| NCT06172803 | exercise/rehabilitation + breathing/respirato | COMPLETED | 4 | 2024-01-15 |
| NCT06726772 | CBT/psychological therapy | COMPLETED | 50 | 2024-01-31 |
| NCT05139979 | breathing/respiratory training | COMPLETED | 189 | 2024-02-07 |
| NCT05677932 | photobiomodulation | COMPLETED | 12 | 2024-02-28 |
| NCT06108297 | lithium | COMPLETED | 5 | 2024-03-26 |
| NCT06952127 | exercise/rehabilitation | COMPLETED | 22 | 2024-03-29 |
| NCT06054438 | herbal/TCM | COMPLETED | 110 | 2024-03-31 |
| NCT05675995 | exercise/rehabilitation | COMPLETED | 40 | 2024-05-15 |
| NCT06251011 | exercise/rehabilitation | COMPLETED | 32 | 2024-06-01 |
| NCT07492953 | EECP | COMPLETED | 10 | 2024-06-24 |
| NCT05507190 | CBT/psychological therapy | COMPLETED | 20 | 2024-06-30 |
| NCT06189066 | Splenic Ultrasound | COMPLETED | 15 | 2024-07-24 |
| NCT06004310 | exercise/rehabilitation | COMPLETED | 120 | 2024-07-25 |
| NCT05633472 | vitamins/minerals | COMPLETED | 33 | 2024-08-23 |
| NCT05212688 | acupuncture/manual therapy | COMPLETED | 119 | 2024-09-25 |
| NCT05630040 | vagus nerve stimulation | COMPLETED | 40 | 2024-10-03 |
| NCT06658340 | CBT/psychological therapy | COMPLETED | 34 | 2024-11-14 |
| NCT06341751 | CBT/psychological therapy | COMPLETED | 18 | 2025-01-15 |
| NCT07429903 | diet | COMPLETED | 31 | 2025-01-24 |
| NCT07363655 | CBT/psychological therapy + pacing/energy man | COMPLETED | 13 | 2025-02-23 |
| NCT05299333 | exercise/rehabilitation | COMPLETED | 30 | 2025-03-07 |
| NCT07745699 | diet | COMPLETED | 11 | 2025-06-01 |
| NCT07279272 | Cold Atmospheric Plasma; Atmospheric Air | COMPLETED | 40 | 2025-07-22 |
| NCT06620406 | microbiome therapy | COMPLETED | 36 | 2025-08-06 |
| NCT06330376 | breathing/respiratory training | COMPLETED | 16 | 2025-08-15 |
| NCT05732285 | CBT/psychological therapy | COMPLETED | 15 | 2025-09-01 |
| NCT06437223 | tonabersat | COMPLETED | 16 | 2025-09-30 |
| NCT01730495 | anti-TNF | TERMINATED | 4 | 2014-08 |
| NCT01368718 | sleep treatment | TERMINATED | 5 | 2015-12 |
| NCT01156922 | rituximab | TERMINATED | 8 | 2016-04 |
| NCT04177459 | exercise/rehabilitation | TERMINATED | 3 | 2020-04-30 |
| NCT06789952 | CBT/psychological therapy | TERMINATED | 41 | 2023-11-30 |
| NCT05813899 | microbiome therapy | TERMINATED | 39 | 2023-12-31 |
| NCT05317975 | exercise/rehabilitation | TERMINATED | 13 | 2024-02-17 |
| NCT06113679 | RDX-19 | TERMINATED | 41 | 2024-05-24 |
| NCT05274477 | fampridine | TERMINATED | 8 | 2024-06-30 |
| NCT07573709 | transcranial stimulation/neurofeedback | TERMINATED | 20 | 2025-05-01 |
| NCT05967052 | gabapentinoid + exercise/rehabilitation + CBT | TERMINATED | 23 | 2025-09-11 |
| NCT01793415 | IodoCarb (r) | UNKNOWN | 40 | 2013-12 |
| NCT01907711 | acupuncture/manual therapy | UNKNOWN | 60 | 2014-01 |
| NCT01966276 | stimulant | UNKNOWN | 134 | 2014-12 |
| NCT02964533 | acupuncture/manual therapy | UNKNOWN | 70 | 2017-02 |
| NCT02847845 | herbal/TCM | UNKNOWN | 50 | 2017-06 |
| NCT03502044 | Inactive Kinetic Oscillation Stimulation (KOS | UNKNOWN | 200 | 2019-12-31 |
| NCT05242094 | exercise/rehabilitation | UNKNOWN | 2000 | 2021-02-02 |
| NCT04599348 | herbal/TCM | UNKNOWN | 70 | 2021-04-30 |
| NCT04678700 | exercise/rehabilitation | UNKNOWN | 40 | 2021-06-01 |
| NCT02854683 | IV saline/fluids | UNKNOWN | 45 | 2021-12 |
| NCT04961333 | exercise/rehabilitation | UNKNOWN | 200 | 2021-12-31 |
| NCT04811859 | breathing/respiratory training | UNKNOWN | 40 | 2022-01-10 |
| NCT04898205 | exercise/rehabilitation | UNKNOWN | 24 | 2022-01-20 |
| NCT05289115 | transcranial stimulation/neurofeedback | UNKNOWN | 30 | 2022-04-14 |
| NCT04846010 | herbal/TCM + vagus nerve stimulation | UNKNOWN | 2000 | 2022-05-01 |
| NCT05289128 | transcranial stimulation/neurofeedback | UNKNOWN | 30 | 2022-06-25 |
| NCT05130736 | exercise/rehabilitation | UNKNOWN | 80 | 2022-07-31 |
| NCT05543408 | exercise/rehabilitation | UNKNOWN | 40 | 2022-08-02 |
| NCT04947488 | diet | UNKNOWN | 50 | 2022-09-18 |
| NCT05119608 | CBT/psychological therapy | UNKNOWN | 40 | 2022-09-30 |
| NCT05150782 | polyphenol/omega-3/anti-inflammatory suppleme | UNKNOWN | 32 | 2022-10-10 |
| NCT05217771 | High tone power therapy group | UNKNOWN | 56 | 2022-11-20 |
| NCT05582603 | exercise/rehabilitation + cognitive rehabilit | UNKNOWN | 20 | 2022-11-30 |
| NCT05516550 | exercise/rehabilitation | UNKNOWN | 412 | 2022-12 |
| NCT05667077 | amantadine | UNKNOWN | 83 | 2022-12-29 |
| NCT05665205 | exercise/rehabilitation | UNKNOWN | 20 | 2023-01-30 |
| NCT05722730 | exercise/rehabilitation | UNKNOWN | 60 | 2023-03-30 |
| NCT05397626 | biofeedback + molecular hydrogen | UNKNOWN | 39 | 2023-05-23 |
| NCT05419219 | exercise/rehabilitation | UNKNOWN | 200 | 2023-05-31 |
| NCT05012826 | exercise/rehabilitation | UNKNOWN | 104 | 2023-06-06 |
| NCT05799534 | exercise/rehabilitation | UNKNOWN | 400 | 2023-07 |
| NCT04532827 | exercise/rehabilitation | UNKNOWN | 105 | 2023-07-01 |
| NCT05934279 | exercise/rehabilitation | UNKNOWN | 40 | 2023-07-30 |
| NCT05749757 | acupuncture/manual therapy | UNKNOWN | 84 | 2023-07-31 |
| NCT05040893 | exercise/rehabilitation | UNKNOWN | 15 | 2023-08-09 |
| NCT05534997 | exercise/rehabilitation | UNKNOWN | 60 | 2023-09-01 |
| NCT05289154 | exercise/rehabilitation | UNKNOWN | 200 | 2023-09-01 |
| NCT05817032 | exercise/rehabilitation | UNKNOWN | 22 | 2023-10 |
| NCT05629897 | cognitive rehabilitation | UNKNOWN | 88 | 2023-11 |
| NCT06144320 | acupuncture/manual therapy | UNKNOWN | 40 | 2023-11 |
| NCT05629884 | exercise/rehabilitation | UNKNOWN | 56 | 2023-11-14 |
| NCT05501067 | exercise/rehabilitation | UNKNOWN | 60 | 2023-11-30 |
| NCT05360563 | exercise/rehabilitation | UNKNOWN | 94 | 2023-12 |
| NCT05975034 | microbiome therapy | UNKNOWN | 240 | 2023-12-11 |
| NCT04836351 | exercise/rehabilitation | UNKNOWN | 20 | 2023-12-30 |
| NCT05855356 | exercise/rehabilitation | UNKNOWN | 50 | 2023-12-30 |
| NCT05212467 | brain-retraining programme | UNKNOWN | 360 | 2023-12-31 |
| NCT05669261 | Adipose Tissue Harvest; ATCell | UNKNOWN | 20 | 2023-12-31 |
| NCT05077241 | breathing/respiratory training | UNKNOWN | 10 | 2024-01 |
| NCT05684952 | herbal/TCM | UNKNOWN | 152 | 2024-01-15 |
| NCT05715801 | hyperbaric oxygen | UNKNOWN | 80 | 2024-01-31 |
| NCT05725538 | exercise/rehabilitation | UNKNOWN | 100 | 2024-03-15 |
| NCT06128941 | exercise/rehabilitation | UNKNOWN | 60 | 2024-04 |
| NCT05668039 | EECP | UNKNOWN | 32 | 2024-04-01 |
| NCT06348186 | acupuncture/manual therapy | UNKNOWN | 30 | 2024-04-01 |
| NCT05851846 | brain-retraining programme | UNKNOWN | 130 | 2024-05-15 |
| NCT06042777 | herbal/TCM | UNKNOWN | 100 | 2024-05-31 |
| NCT05381779 | breathing/respiratory training + exercise/reh | UNKNOWN | 60 | 2024-06-05 |
| NCT06374446 | immersive VR | UNKNOWN | 56 | 2024-06-15 |
| NCT05616806 | iENDURE | UNKNOWN | 10 | 2024-07-01 |
| NCT06074848 | transcranial stimulation/neurofeedback + exer | UNKNOWN | 48 | 2024-07-01 |
| NCT05019963 | exercise/rehabilitation | UNKNOWN | 106 | 2024-08 |
| NCT06145867 | photobiomodulation | UNKNOWN | 10 | 2024-08-01 |
| NCT05231408 | exercise/rehabilitation | UNKNOWN | 30 | 2024-08-15 |
| NCT05481177 | ivabradine | UNKNOWN | 250 | 2024-09-01 |
| NCT05831839 | cognitive rehabilitation + exercise/rehabilit | UNKNOWN | 50 | 2024-09-01 |
| NCT05454683 | sleep treatment + magnesium | UNKNOWN | 106 | 2024-09-30 |
| NCT06091384 | breathing/respiratory training | UNKNOWN | 40 | 2024-09-30 |
| NCT06147050 | metformin | UNKNOWN | 16 | 2024-10 |
| NCT06227273 | molecular hydrogen | UNKNOWN | 50 | 2024-11-01 |
| NCT05932797 | CBT/psychological therapy | UNKNOWN | 60 | 2024-11-18 |
| NCT06211686 | exercise/rehabilitation | UNKNOWN | 50 | 2024-12-01 |
| NCT05918965 | vagus nerve stimulation | UNKNOWN | 45 | 2025-03-01 |
| NCT06142240 | exercise/rehabilitation | UNKNOWN | 40 | 2025-04-01 |
| NCT05947617 | VIX001 | UNKNOWN | 9 | 2025-04-01 |
| NCT06231225 | exercise/rehabilitation | UNKNOWN | 100 | 2025-06-30 |
| NCT06004362 | Complementary and Integrative Medicine online | UNKNOWN | 120 | 2025-07 |
| NCT06492798 | cell therapy (MSC/cord blood) | UNKNOWN | 76 | 2025-07-01 |
| NCT06348212 | microbiome therapy | UNKNOWN | 60 | 2025-08 |
| ISRCTN31455243 | Healthcare Evaluation and Assessment of Patie | see registry | 60 | 2003-07-01 |
| ISRCTN21288417 | CBT/psychological therapy | see registry |  | 2003-10-01 |
| ISRCTN37461058 | CBT/psychological therapy | see registry | 153 | 2004-01-31 |
| ISRCTN31632033 | CBT/psychological therapy | see registry | 130 | 2010-04-01 |
| ISRCTN07518149 | diet | see registry | 5 | 2011-12-31 |
| ISRCTN16045700 | CBT/psychological therapy | see registry | 20 | 2020-12-31 |
| ISRCTN58998320 | A novel online program for fatigue; a feasibi | see registry | 200 | 2022-03-31 |
| ISRCTN11868601 | CBT/psychological therapy | see registry | 94 | 2022-04-30 |
| ISRCTN11563542 | COVID VR | see registry | 6 | 2022-07-12 |
| ISRCTN36407216 | CBT/psychological therapy | see registry | 554 | 2023-07-31 |
| ISRCTN16012321 | CBT/psychological therapy + exercise/rehabili | see registry | 40 | 2023-12-31 |
| ISRCTN39864861 | cognitive rehabilitation | see registry | 1608 | 2025-05-30 |
| ISRCTN15068418 | CBT/psychological therapy | see registry | 269 | 2025-06-15 |
| ISRCTN12368131 | herbal/TCM | see registry | 105 | 2025-06-18 |

## Largest and best-designed trials

Blinded placebo- or sham-controlled trials with a verified result, sorted by size (n of 100 or more):

| Trial | Population | Intervention | n | Result | Finding |
|---|---|---|---:|---|---|
| NCT05595369 | long COVID | nirmatrelvir/ritonavir | 951 | null | RECOVER-VITAL: nirmatrelvir/ritonavir 15 or 25 d vs placebo (n=959; cognitive, autonomic, exercise phenotypes): no benefit in any phenotype |
| NCT05074888 | long COVID | Prospekta | 676 | positive | Prospekta vs placebo for post-COVID asthenia (n=676): FSS fell 16.5 vs 14.0 (p=0.0016); small absolute difference; industry-sponsored |
| PMID:15353532 | ME/CFS | galantamine | 434 | null | Galantamine (3 doses) vs placebo (n=434): no benefit on any outcome |
| PMID:18277062 | ME/CFS | distant healing | 409 | null | EUHEALS distant healing vs no healing, partially blinded (n=409): no effect on mental or physical health |
| NCT06128967 | both | fluvoxamine + metformin | 399 | mixed | Fluvoxamine and metformin vs placebo, adaptive (n=399): fluvoxamine reduced fatigue (day 60 and 90) and improved QoL; metformin no benefit |
| KCT0009342 | long COVID | metformin / ursodeoxycholic acid | 396 | null | Metformin or ursodeoxycholic acid vs placebo 2 weeks (n=396, Korea): no improvement in recovery |
| CTRI/2021/11/038234 | long COVID | colchicine | 346 | null | Colchicine vs placebo (n=346, India): no improvement in 6MWT, respiratory function or inflammation at 52 weeks |
| NCT05965739 | long COVID | cognitive rehabilitation + transcranial  | 328 | null | RECOVER-NEURO: BrainHQ, PASC-CoRE cognitive rehab and tDCS vs controls (n=328): no benefit of any intervention on cognition |
| NCT03186027 | ME/CFS | coenzyme Q10 + vitamins/minerals | 282 | mixed | CoQ10 200 mg + NADH 20 mg vs placebo for 12 weeks (n=207 analysed): within-group falls in cognitive fatigue and FIS-40 and better SF-36 and sleep on active; abstract repo |
| NCT05999435 | long COVID | LAU-7b (fenretinide) | 272 | null | LAU-7b (fenretinide) vs placebo (n=272): no SF-36 PCS improvement at week 12 |
| NCT05497089 | long COVID | temelimab | 203 | null | Temelimab (anti-HERV-W Env) vs placebo (n=203): no improvement in fatigue |
| NCT04652518 | long COVID | (deu)pirfenidone | 177 | null | Deupirfenidone (n=185): 6MWD improved in both arms; no difference vs placebo |
| NCT05430152 | both | low-dose naltrexone | 160 | null | LDN up to 4.5 mg vs placebo for 16 wk (n=160): no difference in fatigue (FSS), pain, symptoms, step count or QoL |
| NCT02807649 | ME/CFS | herbal/TCM | 159 | mixed | Cistanche+ginkgo product (n=159): symptom improvement vs placebo but effect size described as trivial |
| NCT01318356 | ME/CFS | CBT/psychological therapy + antibiotic | 155 | mixed | Qure: doxycycline or CBT vs placebo in Q-fever fatigue syndrome (n=155): posted CIS fatigue CBT 31.6, doxycycline 40.8, placebo 37.8; no between-group statistics posted;  |
| NCT05576662 | long COVID | nirmatrelvir/ritonavir | 155 | null | STOP-PASC: nirmatrelvir/ritonavir 15 d vs placebo/r (n=155): safe, no benefit on core symptoms |
| NCT05890534 | long COVID | polyphenol/omega-3/anti-inflammatory sup | 153 | null | PYCNOVID: Pycnogenol 200 mg vs placebo (n=153): no difference in EQ-VAS; both groups improved |
| NCT02229942 | ME/CFS | rituximab | 151 | null | RituxME (n=151): response 26% rituximab vs 35% placebo (P=0.22); no difference in fatigue over 24 months |
| PMID:23981369 | ME/CFS | acupuncture | 150 | positive | Acupuncture vs sham (n=~150, single-blind): moderate reductions in physical and mental fatigue |
| NCT05047952 | long COVID | vortioxetine | 149 | mixed | Vortioxetine vs placebo (n=149): DSST cognition change 0.33 vs 0.18, not significant unadjusted; depression and HRQoL improved in paper |
| PMID:28408193 | ME/CFS | imipramine | 138 | positive | STreSS-3: imipramine vs placebo in multiple functional somatic syndromes incl. CFS (n=138): better overall health |
| NCT05874037 | long COVID | fluvoxamine | 137 | positive | Fluvoxamine vs placebo (n=137): total symptom score fell 47.3 vs 31.1 (no p-value posted) |
| PMID:9828987 | ME/CFS | antidepressant + graded exercise therapy | 136 | mixed | Fluoxetine and graded exercise, factorial placebo-controlled (n=136): exercise reduced fatigue and improved work capacity; fluoxetine improved depression only |
| ISRCTN23930966 | ME/CFS | Sijunzi decoction (TCM) | 127 | positive | Sijunzi decoction vs placebo (n=127): larger fatigue improvement with good adherence; gut microbiota changes |
| NCT04960215 | long COVID | coenzyme Q10 | 121 | null | High-dose CoQ10 500 mg vs placebo, crossover (n=121): no reduction in PCC symptoms |
| NCT01040429 | ME/CFS | clonidine | 120 | negative/harm | NorCAPITAL clonidine (n=120 adolescents): lowered sympathetic activity and inflammation but reduced physical activity (primary steps/day); judged not clinically useful |
| NCT05922865 | long COVID | KNEESUP smart knee assistive device + KN | 120 | positive | Exercise telerehabilitation with smart sensor vs active control (n=120): VO2peak +2.1 and gait speed improved; QoL and sleep unchanged |
| PMID:40094666 | long COVID | cholinesterase inhibitor/memantine + CBT | 120 | null | Donepezil vs placebo (n=120): no benefit for fatigue or psychological symptoms |
| DOI:10.3310/nihropenres.14236.1 | long COVID | auricular vagus nerve stimulation | 114 | null | PAuSing-pCF: percutaneous auricular VNS vs sham vs placebo (n=114): no difference in fatigue at 8 weeks |
| PMID:8622391 | ME/CFS | fluoxetine | 107 | null | Fluoxetine 20 mg vs placebo 8 wk (n=107): no benefit on any CFS dimension, including depressed patients |
| PMID:15016577 | ME/CFS | homeopathy | 103 | mixed | Individualised homeopathy vs placebo (n=103): better on MFI general fatigue, not other scales; weak, equivocal |
| KCT0007515 | ME/CFS | Gongjin-Dan / Ssanghwa-Tang herbal medic | 103 | mixed | Gongjin-Dan and Ssanghwa-Tang vs placebo in chronic fatigue (n=103): some subscale improvements only |
| PMID:11150109 | ME/CFS | fludrocortisone | 100 | null | Fludrocortisone vs placebo for neurally mediated hypotension (n=100): no more effective than placebo |
| PMID:12829200 | ME/CFS | hydrocortisone + fludrocortisone | 100 | null | Hydrocortisone + fludrocortisone vs placebo, crossover (n=100): no improvement in fatigue or well-being |
| NCT05668091 | long COVID | nirmatrelvir/ritonavir | 100 | null | PAX LC: nirmatrelvir/ritonavir 15 d vs placebo/r (n=100, decentralised): no improvement at day 28 |
| PMID:35337143 | long COVID | ADAPT-232 adaptogen herbal | 100 | mixed | Adaptogens (Chisan/ADAPT-232) vs placebo: better physical performance and shorter fatigue duration reported; most outcomes not different |

Large unblinded RCTs (n of 200 or more) with a verified result:

| Trial | Population | Intervention | n | Result | Finding |
|---|---|---|---:|---|---|
| ISRCTN54285094 | ME/CFS | CBT/psychological therapy + pacing/energ | 641 | mixed | PACE: APT, CBT, GET vs specialist medical care (n=641): CBT and GET moderately improved fatigue and function vs SMC; APT no benefit. Outcome changes and effect sizes cont |
| ISRCTN11466448 | long COVID | online group rehabilitation | 585 | positive | REGAIN: online group physical and mental rehab vs usual care (n=585): better HRQoL at 3 and 12 months |
| PMID:35416437 | long COVID | Qingjin Yiqi granules (TCM) | 388 | mixed | Qingjin Yiqi granules vs control (n=388, PCC): breathlessness (mMRC, Borg) improved; overall symptom score not different |
| NTR4809 | ME/CFS | internet CBT stepped care | 360 | null | Stepped care with internet CBT vs face-to-face CBT (n=~360): noninferior; less therapist time |
| NCT05196451 | long COVID | CBT/psychological therapy | 314 | positive | SIPCOV brief outpatient rehab (n=314): effective on primary physical functioning vs usual care |
| NCT01650636 | ME/CFS | CBT/psychological therapy | 300 | mixed | Video CBSM vs health information: medium-large benefits only in high-PEM subgroup, none in low-PEM |
| PMID:22370004 | ME/CFS | CBT/psychological therapy + graded exerc | 300 | null | Counselling or graded exercise vs usual care in primary-care chronic fatigue (n=~300): no significant difference in fatigue |
| NCT06404073 | long COVID | pacing/energy management | 298 | positive | RECOVER-ENERGIZE structured pacing vs usual care (n=298): PEM frequency (-0.46, CI -0.65 to -0.27) and severity reduced |
| PMID:20418251 | ME/CFS | pragmatic rehabilitation / supportive li | 296 | null | FINE trial: nurse-led pragmatic rehabilitation or supportive listening vs GP care (n=296): small short-term fatigue gain, not significant at 70 weeks; listening ineffecti |
| PMID:11265953 | ME/CFS | CBT | 278 | positive | CBT vs support groups vs natural course (n=278): CBT better for fatigue and functional impairment at 14 months |
| ISRCTN23962803 | ME/CFS | graded exercise therapy + pacing/energy  | 241 | null | MAGENTA: GET vs activity management in paediatric ME/CFS: no evidence GET more effective; little improvement in either arm |
| ISRCTN18020851 | ME/CFS | CBT/psychological therapy + pacing/energ | 240 | mixed | FITNET-NHS: remote CBT vs activity management in paediatric ME/CFS (n=~240): better physical function and school attendance at 6 months; limited sustained effect, not cos |
| NTR4013 | ME/CFS | internet CBT | 240 | positive | Web-based CBT (two feedback formats) vs wait list (n=240): fatigue reduced in both iCBT arms |
| NCT06245642 | ME/CFS | herbal/TCM | 235 | positive | Compound Ciwujia granules vs placebo (TCM, n=235): fatigue scores improved significantly at 6 weeks |
| ISRCTN12595520 | long COVID | remote weight management (diet) | 234 | positive | ReDIRECT remote weight management vs wait list in overweight long COVID (n=234): small improvement in most troublesome symptom |
| NCT04900961 | long COVID | exercise/rehabilitation | 233 | positive | Personalised resistance exercise, 12 weeks: improved exercise capacity, HRQoL and grip strength; PEM not increased |
| ISRCTN22975026 | ME/CFS | graded exercise therapy + exercise/rehab | 211 | positive | GETSET: guided graded exercise self-help + SMC vs SMC (n=211): lower fatigue at 12 weeks; smaller effect on disability |
| PMID:26402868 | ME/CFS | group CBT | 204 | positive | Group CBT vs wait list (n=204): large improvements in fatigue (ES 1.1) and impairment |

Large trials that have finished or are about to report and are worth watching: RECOVER-AUTONOMIC IVIG (NCT06305793, n=200) and ivabradine (NCT06305806, n=181), the probiotic trial PURE-LC (NCT06643299, n=194), the tVNS trial NCT07190105 (n=546), the LIFT trial of pyridostigmine and LDN (NCT06366724, n=160), the Mount Sinai sirolimus trial (NCT06960928, n=80), the ResetME daratumumab RCT (CTIS 2024-520094-13-00, n=66), and the GLP-1 trial NCT07128082 (n=1,000).

## Registries

| Registry | Reached? | What was done |
|---|---|---|
| ClinicalTrials.gov (API v2) | yes | Pulled in full: condition and term searches for ME/CFS and long COVID |
| Europe PMC (REST, core) | yes | Searched for trial publications; abstracts were read for the judgements |
| ISRCTN | yes (XML API) | Searched by condition keyword only: 71 unique records; 21 merged, 22 added, 28 excluded |
| EU CTIS | yes (POST search API; GET returns 403) | 5 keyword searches: 28 records; 10 merged, 6 added, 12 excluded. Status codes are numbers and were not mapped |
| DRKS (Germany) | page loads, but no data | Search results are drawn by JavaScript and are not in the HTML; the guessed API URL returned 404. German trials appear here only when they are also on CT.gov, CTIS or in a paper |
| WHO ICTRP | page loads, but cannot be queried | The search is an ASP.NET form (__VIEWSTATE), so curl cannot run it. Chinese, Iranian, Indian and Japanese registries were reached only through IDs cited in papers |
| ANZCTR | no (HTTP 403) | Not searched |
| ChiCTR | no (HTTP 405) | Not searched; ChiCTR trials appear only when published |

## Coverage check against the treatments page

This compares the database with `site-src/data/treatments.json` (47 treatments, updated 2026-10-28), which was read but not edited. The same content is in `coverage_check.json`.

### Page entries whose evidence is out of date, wrong or incomplete

- **pacing-heart-rate**: The page says the only controlled test was PACE's APT arm. RECOVER-ENERGIZE structured pacing vs usual care (NCT06404073, n=298, long COVID) reduced how often PEM happened (-0.46, CI -0.65 to -0.27) and how severe it was. It is the first positive randomised trial of pacing.
- **low-dose-naltrexone**: The page calls the Nacul trial 'presented May 2026, unpublished'. It is now a medRxiv preprint, Oct 2026 (DOI 10.64898/2026.10.22.26363582; NCT05430152): no difference in fatigue (FSS), pain, symptoms, steps or quality of life at 16 weeks. The cited reference should be updated.
- **rintatolimod**: The page cites a 2024 company release for the long COVID trial. ClinicalTrials.gov now has posted results (NCT05592418, n=80): the change in PROMIS fatigue was essentially the same on drug and placebo. The 1994 placebo-controlled RCT (n=92, PMID 8148460, positive on Karnofsky score and treadmill work) is not cited.
- **tvns**: PAuSing-pCF is listed as a pending preprint. It is now published in NIHR Open Research (DOI 10.3310/nihropenres.14236.1, n=114): fatigue at 8 weeks did not differ from sham or placebo. It should move from 'pending' into the evidence line. Also pending: NCT07190105 (n=546, finished Dec 2025, not reported).
- **mito-supplements**: The page says 'creatine has only a 12-person long-COVID RCT'. That is out of date: a 3-arm placebo-controlled RCT (n=67, PMID 41958909) found that 6 g/day improved fatigue and muscle strength (18 g did not). There are also two other small long COVID creatine RCTs (PMID 38099612, 38684388) and an open-label ME/CFS MRS study (PMID 39408275, n=11). NCT05703074 (nicotinamide riboside plus mind-body therapy, n=310, completed Sep 2024) has not reported.
- **immunoadsorption**: The page says plasma exchange has only case reports. There is a sham-controlled plasma exchange RCT (NCT05445674, n=50, posted results): safe, but no difference on any efficacy outcome. That adds to the negative evidence. On BC007, the page says it failed a 114-person unpublished company trial. The database has the 30-person reCOVer crossover published as positive (CTIS 2022-001781-35) and the phase 2 RCT NCT05911009 (n=119, completed Jun 2024) with no posted results and no paper. The page should cite the registry record and make clear the failure is not yet in a publication or posted results.
- **stimulants**: Only methylphenidate is cited. Missing: lisdexamfetamine vs placebo (NCT01071044, n=26), which improved executive function and fatigue, and solriamfetol vs placebo (NCT04622293, n=44), which improved fatigue at week 8 but not week 6. Also NCT01966276 (methylphenidate plus nutrient formula, n=134 planned, completed 2014): no results ever reported.
- **cbt-supportive**: The evidence line rests on PACE alone. The database has about 40 other randomised CBT or psychological trials in ME/CFS: FITNET-NHS (n=~240), FINE (n=296, null at 70 weeks), Dutch group, internet and stepped-care CBT, and CBSM, which helped only the high-PEM subgroup. Most are unblinded and use self-rated outcomes, and many are positive. The page's conclusion can stand, but it should say that the wider body exists and why it is discounted (unblinded subjective outcomes, few objective measures).
- **graded-exercise**: The page could cite MAGENTA (ISRCTN23962803, paediatric, GET no better than activity management) and GETSET (ISRCTN22975026, guided self-help GET, a self-rated gain). PACE alone is not the full trial base.
- **ivig**: This matches the database. Add the terminated subcutaneous IgPro20 trial in post-COVID POTS (NCT06524739, n=16, uninterpretable). RECOVER-AUTONOMIC IVIG (NCT06305793, n=200) has finished, with no results posted or published.
- **ivabradine**: The page cites the RECOVER-AUTONOMIC ACC 2026 presentation. The database has no paper or posted results for NCT06305806 yet, so this could not be checked. NCT05481177 (ivabradine nested RCT, n=250, completed Sep 2024) has not reported.
- **steroids**: There are two unreported long COVID steroid RCTs: prednisolone plus B vitamins (NCT05638633, n=321, completed Jan 2025) and methylprednisolone vs placebo (NCT05986422, terminated, n=96).
- **fmt**: FMT matches the database. The probiotics entry (one 39-person RCT) is missing several placebo RCTs that have not reported: NCT04741841 (GutMagnific, ME/CFS, n=74, completed 2021), NCT05975034 (long COVID, n=240 planned, completion Dec 2023, status UNKNOWN), NCT05813899 (L. paracasei PS23, terminated). Also PURE-LC (NCT06643299, n=194), completed Jul 2026.
- **herpes-antivirals**: The valacyclovir plus celecoxib pilot RCT (NCT06316843, n=59) was completed Oct 2024 and has not reported. The page mentions only the planned trial NCT07597902.
- **sleep-medicines**: Add the sodium oxybate placebo RCTs (SAFFE NCT02055898, n=13, more slow-wave sleep but no daytime gain; NCT00498485, terminated, n=17). A second melatonin plus zinc RCT (NCT05454683, n=106 planned, completion Sep 2024) has not reported.
- **iv-saline**: The page says there is no randomised trial. Two have not reported: NCT04740736 (ME/CFS, saline infusion vs placebo in a PEM study, n=106, completed Jul 2025) and NCT02854683 (oral rehydration vs IV saline crossover, status UNKNOWN, completion 2021).
- **eecp**: This matches the database. NCT05668039 is still not reported (status UNKNOWN, primary completion Apr 2024).
- **stellate-ganglion-block**: The page cites STAR-CO NCT06253806, which is not in the database (it was not returned by the ME/CFS or long COVID searches; check it is the right ID). The database has two placebo-controlled RCTs that are recruiting (NCT07374562; a cervicothoracic sympathetic chain block vs saline, NCT07468604) and an open-label dual-block study, NCT05638620 (n=20, COMPASS-31 down 38%).

### Treatments with trial evidence and no page entry

- **fluvoxamine** (positive/mixed): Two placebo RCTs, both favourable: an adaptive fluvoxamine/metformin trial (NCT06128967, n=399), where fluvoxamine reduced fatigue at days 60 and 90, and NCT05874037 (n=137), with a bigger fall in symptom score but no p-value posted. Another RCT is recruiting (NCT07359482).
- **metformin (treatment of established long COVID)** (null): Null in the adaptive trial (NCT06128967) and in a Korean RCT (KCT0009342, n=396). NCT07280572 (n=1,250) is enrolling.
- **colchicine** (null/mixed): Null in a placebo RCT (CTRI/2021/11/038234, n=346). In STIMULATE-ICP (open-label) it reduced fatigue slightly. NCT07280572 is ongoing.
- **temelimab (anti-HERV-W)** (null): Placebo RCT, n=203: no improvement in fatigue.
- **efgartigimod (FcRn blocker)** (null): Placebo RCT in post-COVID POTS, n=53: no benefit (COMPASS-31 -14.4 vs -16.0). The extension study was terminated.
- **cell therapy (cord blood, MSC)** (positive (small)): REGENECYTE cord blood vs placebo (NCT05682560, n=30): positive on fatigue, but small and industry-run. NCT06492798 (n=76) is overdue. Five MSC trials were withdrawn.
- **Prospekta (Russian ultra-low-dose antibody product)** (positive (small effect)): Placebo RCT for post-COVID asthenia, n=676: FSS fell 16.5 vs 14.0. The difference is significant but small, and the trial was industry-sponsored. This is the largest positive blinded trial in the database.
- **vortioxetine** (mixed): Placebo RCT, n=149: no significant gain on the primary cognitive measure (DSST).
- **lithium (low dose)** (null): Placebo RCT, n=52: no effect on fatigue or brain fog.
- **Pycnogenol and other polyphenols** (mixed): PYCNOVID (n=153) was null. Smaller positive RCTs: quercetin phytosome (n=78), cocoa flavonoids (FLALOC, n=46), high-polyphenol chocolate (n=10). The page's fibre-polyphenols entry says there is no treatment trial.
- **cognitive rehabilitation and tDCS** (null): RECOVER-NEURO (NCT05965739, n=328): no benefit from BrainHQ, PASC-CoRE or tDCS. Single-centre tDCS trials are split (NCT04876417 null; NCT05389592 small gains).
- **amantadine** (mixed): Open-label RCT positive on FSS (PMID 38228731). The placebo RCT NCT06055244 (n=64) finished in Dec 2025 and has not reported. A 1997 ME/CFS RCT found it poorly tolerated.
- **galantamine** (null): Placebo RCT in ME/CFS, n=434: no benefit on any outcome.
- **LAU-7b (fenretinide), deupirfenidone** (null): Both null in placebo RCTs of long COVID (n=272 and n=177).
- **structured rehabilitation programmes for long COVID** (positive (unblinded)): REGAIN online group rehab (ISRCTN11466448, n=585) improved quality of life. SIPCOV (n=314) and personalised resistance exercise (n=233) were positive, and PEM did not increase. All were unblinded. The page has no long COVID rehab entry separate from GET.
- **plitidepsin** (pending): A phase 3 long COVID RCT is registered in CTIS (2023-504087-42-00, n=260 planned). No results yet.

These entries match the database and need no change: fludrocortisone, clonidine, rituximab, cyclophosphamide, daratumumab, covid-antivirals, oxaloacetate, coq10-nadh, hbot-clinical, red-light, mast-cell, microclot-therapy, rapamycin, pyridostigmine. Entries with no trials to check: upright-time, oral-fluids-salt, hbot-home, bpc-157, nicotine-patch, thyroid-testosterone, desmopressin, guanfacine-nac, midodrine (the database holds only a terminated droxidopa study), and glp1 (the trials are all ongoing).

## Data-quality caveats

- Screening used keyword rules on registry fields, followed by manual review of borderline and long COVID records. Some trials may be wrongly in or out, especially long COVID trials with mixed respiratory and fatigue aims.
- Results were judged from abstracts and posted results only, not full texts. A 'positive' label means the abstract reports benefit on its main outcome; it does not mean the trial was free of bias. Many positive trials are small, unblinded or industry-run.
- 56 verified efficacy results come from TCM, herbal, acupuncture or moxibustion trials. These are mostly Chinese, unblinded, and use 'total effective rate' outcomes, and they push up the positive count for the 'supplement' and 'other' classes.
- 26 records rest only on preprints (`preprint_only: true`), and 12 publication records have no abstract (result not-assessable).
- The 'overdue' count may overstate non-reporting. A results paper may exist that Europe PMC did not match to the trial (for example, no registry ID in the abstract, or a non-English paper). UNKNOWN-status CT.gov records past their completion date are counted as overdue.
- ISRCTN coverage is partial (condition keyword search only). CTIS status codes are not mapped. DRKS, ANZCTR, ChiCTR and ICTRP could not be searched, so non-US/EU trials appear only when published or cross-registered.
- Intervention names and classes are normalised with a keyword dictionary (`scripts/build_trials.py`). Model steps are assigned by intervention, not by trial outcome, and combination trials take the first named active drug's class.
- Population 'both' means the trial enrolled or targeted both ME/CFS and long COVID, including post-COVID ME/CFS.
- n_actual comes from CT.gov ACTUAL enrolment, or from the abstract for publication-only records; abstract n is sometimes the number analysed rather than randomised.
