# Weekly literature update (run from /mnt/project-files/mecfs-research)

1. `python3 -I lit/scripts/fetch_lit.py` (Europe PMC: PubMed, PMC, medRxiv, bioRxiv, Research Square). About 10 minutes.
2. `python3 -I lit/scripts/triage.py` (dedupe, kind, domains, cited flag; only curated atlas/trials/treatments evidence counts as cited).
3. `python3 -I lit/scripts/new_batches.py` (writes lit/batches/n<date>-NN.jsonl for papers never screened; prints none if nothing new).
4. Screen each new batch with an agent using lit/SCREEN_BRIEF.md, output to lit/screen/n<date>-NN.json. Use small batches (10 to 25 papers per agent).
5. `python3 -I lit/scripts/make_verify.py 10` writes lit/verify/v-NN.jsonl for links not yet verified (rename or clear old v-*.jsonl first). Run one Opus agent per batch with lit/VERIFY_BRIEF.md, output lit/verify/out-<date>-NN.json. integrate.py uses only verified links once any verdicts exist.
6. `python3 -I lit/scripts/integrate.py --apply` (rebuilds the "sweep" lists, writes grade_review.md and new_candidates.md).
7. Read new_candidates.md and grade_review.md; add or regrade findings by hand in atlas/build/research/*.json (add "added" and "regrade_note").
8. `bash site-src/build.sh deploy`; check "Current Version ID".
