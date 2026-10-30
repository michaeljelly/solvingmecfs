# INPUTS_NOTE: deep2/singlecell_methyl

Working copies live in /tmp/work_scm (local disk); results/, scripts/, figures/ are in this folder. Large raw files were not copied to the shared folder.

## 1. GSE214283 (Hanson/Grimson, Cornell) PBMC 10x scRNA
- Source: ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214283/suppl/GSE214283_RAW.tar (1.3 GB, 116 filtered_feature_bc_matrix.h5). Re-used the copy extracted by the previous round at /tmp/work_geo/sc (not re-downloaded).
- 116 samples = 58 people x 2 days (D1 pre-CPET, D2 24 h post); 30 ME/CFS, 28 controls; 18 chromium batches balanced by group.
- After QC (>=200 genes, <6000 genes, >=500 UMI, <20% mito): 503,739 cells. CellTypist Immune_All_Low (majority voting) + marker gating.
- Per-cell table (/tmp/work_scm/cells_all.parquet, 503,739 rows x ~80 cols) and counts per sample x type in results/sc_celltype_counts_per_sample.csv.
- Surprises: Scrublet returned no scores (doublets not removed). CellTypist labelled almost no Th17/Th1/gamma-delta cells at majority-vote level (<100 cells), so these were marker-gated instead. NK CD56bright cells were absorbed into CD56dim clusters by majority voting (113 cells).
```
(copied data) cells_all.parquet
gsm         type            ct_pred                 pct_mt    lib     ngenes
GSM6603110  Mono_classical  Classical monocytes     0.1125    9703    2524
GSM6603110  CD4_naive_CM    Tcm/Naive helper T cells 0.0633   5621    1454
GSM6603110  NK_CD56dim      CD16+ NK cells          0.1337    5100    1668
```

## 2. DNA methylation (GEO series-matrix betas)
| Accession | Tissue / array | n | Format | Rows x cols |
|---|---|---|---|---|
| GSE156792 (Herrera 2018, PMID 30516085) | CD3+ T cells, 450K | 61 ME / 48 HC | series matrix beta (minfi), 429 MB gz | 485,512 x 109 |
| GSE93266 (de Vega 2017, PMID 28231836) | "PBMC", 450K | 49 / 26 | series matrix beta, 185 MB gz | 485,512 x 75 |
| GSE111183 (Trivedi 2018, PMID 30036399) | PBMC, EPIC | 13 / 12, all female | series matrix beta, 66 MB gz | 866,895 x 25 |
| GSE327255 (2026 saliva, PMID 42010606) | saliva, EPIC | 54 / 21 | GSE327255_MatrixProcessedGEO_beta.csv.gz, 233 MB | 710,267 x 75 |
- Also downloaded but not analysed (small n / RRBS): GSE153667 (RRBS 10/10), GSE166592 (RRBS, 3 people x 5 time points), GSE304805 (RRBS PEM 0/24/48 h, 5 ME / 2 HC).
- Annotation: GPL13534 and GPL21145 manifests from GEO platform suppl.
- Clock coefficients: biolearn 0.x data folder (Horvath 2013, Horvath skin-blood 2018, Hannum 2013, Levine PhenoAge 2018, DunedinPACE 2022, Lu DNAmTL 2019, GrimAge v1 component CpGs, McCartney 2018 smoking/alcohol/BMI, Zhang 2017 mortality), Salas 12-cell and Reinius 6-cell reference matrices.
- Surprises: no age in GSE156792/GSE93266/GSE111183 metadata (only saliva has age). GSE93266 "PBMC" contains large amounts of granulocyte DNA (median 48% in controls vs 14% in cases), and sample G00199 appears twice (G00199, G00199_2). Saliva matrix has no X/Y probes. The saliva paper's lead CpG cg19803194 (PTPRN2) has a SNP (rs4909073) within 10 bp of the probe.
```
(copied data) GSE156792 betas
ID_REF      GSM4744443 GSM4744444 GSM4744445 GSM4744446
cg00000029  0.732      0.698      0.719      0.722
cg00000108  0.923      0.912      0.905      0.927
```

## 3. Replication / repertoire data
- GSE339049 (severe long COVID monocytes, 2026): GSE339049_Counts.tsv.gz is an ATAC peak count matrix (chr:start-end rows), 48 samples; not scRNA. Not used further.
- GSE235050 (long COVID scRNA, 12 samples): only per-cluster cell counts with vague labels (e.g. "ISG expressing immune cells"); not usable for replication.
- Bulk RNA re-used from previous round (/tmp/work_geo/raw): GSE227375, GSE268212 (sorted CD8, 7/7, same Cornell group), GSE98139, GSE128078, GSE270045 (long COVID whole blood 19/17).
- VDJdb 2026-06-03 release (Zenodo record 22104776, vdjdb-2026-06-03.zip, 41 MB): 120,575 human TRB entries (EBV 20,244; CMV 30,380).
- PRJEB83659 (ENA): BCR heavy-chain amplicons, UK ME/CFS Biobank, 123 samples (37 severe ME/CFS, 31 mild-moderate, 24 healthy, 31 MS), 495 MiSeq runs, 6.1 GB. Libraries A and B = replicate PCRs; reads start with an IgM (GGGCGGATGCACTCC) or IgG (GGCCCTTGGTGGAGG) constant-region primer. First 10,000 R1 reads of one run per library were streamed.
- PRJNA940635 (SRA/ENA): TCR-locus hybrid-capture on genomic DNA from sorted CD4 and CD8 T cells, UK ME/CFS Biobank, 304 samples (74 severe, 75 mild-moderate, 77 healthy, 78 MS; all female; age bands), 150 GB. First 150,000 R1 reads per run streamed; reads with a TRBJ core 20-mer kept.
- niaid_pi-mecfs-iggseq repo (local): this is IgG-bound gut microbiota 16S (IgG-seq), not a B-cell receptor repertoire; not used here.

## 4. Parsed repertoire examples
```
(copied data) BCR parsed reads (PRJEB83659)
sample      lib iso        v     j            cdr3
SBL046 SBL046-A IgG IGHV4-39 IGHJ5 ARDGRQWTPQGGFDP
SBL046 SBL046-A IgM IGHV3-33 IGHJ6            None
SBL046 SBL046-A IgM IGHV3-48 IGHJ5        ARYSWFDP
SBL046 SBL046-A IgM IGHV3-13 IGHJ4  ARGGRGSPQWAYDY
```
```
(copied data) TCR parsed reads (PRJNA940635)
      sample        v       j             cdr3
SAMN33477952 TRBV11-3 TRBJ2-7 CASRLDVGQSRCEQYF
SAMN33477952 TRBV12-4 TRBJ2-7  CASSLDRGLSSEQYF
SAMN33477952  TRBV7-6 TRBJ2-7  CIAVPEAVYEDEQYF
SAMN33477952   TRBV28 TRBJ2-7   CASSFFDGSYEQYF
```
- BCR: 246 runs (one per replicate library), ~8,900 of 10,000 reads assigned a V gene per run, ~62% with a CDR3. TCR: 304 runs, median ~1,450 productive TRB CDR3 reads per sample from 150k reads.
