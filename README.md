# Solving ME/CFS

Source for [solvingmecfs.org](https://solvingmecfs.org): an evidence atlas of ME/CFS research, a working model of the mechanism, open questions, proposed experiments, a register of treatment trials, and the re-analyses of public data behind several of the site's claims.

Written by [Michael Jelly](https://michaeljelly.com), with AI agents assisting with literature review and analysis.

## Layout

| Path | What it is |
|---|---|
| `site-src/data/graph.json` | The links between model steps, atlas findings, open questions, proposed experiments and treatments. Edit this to change the model. |
| `site-src/data/treatments.json` | 61 treatments: model step, rationale, evidence, risks, cost, sources. |
| `site-src/data/trials.json`, `funding.json` | Trial register and funding calls. |
| `site-src/pages/` | Page fragments. Blocks between `<!--gen:NAME-->` markers are generated; everything else is hand-written. |
| `site-src/*.py` | Generators: `make_graph.py`, `make_treatments.py`, `make_trials.py`, `make_funding.py`, `site_build.py`, `figs/gen_figs.py`. |
| `site-src/check_graph.py` | Fails the build if a link in the graph points at nothing, a treatment lacks a model rationale, or the same finding carries two grades. |
| `site-src/verify_refs.py` | Checks unverified references against Crossref, PubMed and ClinicalTrials.gov. |
| `atlas/build/research/*.json` | The evidence atlas: findings by body system, graded, with studies for and against, and the hypothesis register. |
| `atlas/build/build_atlas.py` | Builds the atlas page and `atlas-data.json`. |
| `atlas-v2/` | Pipelines for the trial register (`trials/`), atlas gap audit (`gaps/`), landscape of comparable efforts and funding calls. Raw downloads are not included. |
| `analyses/` | Re-analyses of public ME/CFS datasets (DecodeME genetics, proteomics, metabolomics, microbiome, single-cell, Mendelian randomisation, BioMapAI, registry epidemiology). Each folder has a README with inputs and results. Large data files are hosted at `data.solvingmecfs.org` and listed in the [data catalogue](https://solvingmecfs.org/data). |
| `data-catalogue/AUDIT.md` | Licence audit for the rehosted datasets. |
| `site/` | The built site, as deployed. |

## Build

```
bash site-src/build.sh          # check, regenerate every page into site/
bash site-src/build.sh deploy   # also deploy to Cloudflare Workers (needs credentials)
```

Python 3.10+ with the standard library only. The build order is: check the graph, build the atlas, render the linked blocks, draw the figures, then wrap every page.

## Corrections

If a grade in the atlas is wrong, a citation does not say what it is cited for, or a step in the model is contradicted by data we missed, open an issue or write to hello@solvingmecfs.org.
