#!/bin/bash
# Full rebuild of solvingmecfs.org: atlas fragment -> treatments fragment -> wrap all pages -> deploy.
set -e
R="${R:-$(cd "$(dirname "$0")/.." && pwd)}"
cd $R/site-src && python3 -I check_graph.py . $R/atlas/build/research
cd $R/atlas/build && python3 -I build_atlas.py research $R/site $R/site-src
cd $R/site-src && python3 -I make_graph.py . $R/atlas/build/research
cd $R/site-src && python3 -I figs/gen_figs.py svg
cd $R/site-src && python3 -I make_treatments.py . && python3 -I make_funding.py . && python3 -I make_trials.py . $R/site && python3 -I site_build.py . $R/site
if [ "$1" = "deploy" ]; then
  cd /tmp && CLOUDFLARE_ACCOUNT_ID=57807fa4814e54014ff99781d420a7b2 npx -y wrangler@4 deploy --assets $R/site --name mecfs-research --compatibility-date 2026-09-01
  rm -rf $R/site/.wrangler
fi
