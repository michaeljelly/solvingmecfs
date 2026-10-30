import json, sys
from collections import Counter
base = sys.argv[1]
R = [r for r in json.load(open(f"{base}/work/epmc_screen.json")) if r["screen"] == "KEEP"]
print(Counter((r["isme"], r["islc"], r["trialpt"]) for r in R))
R.sort(key=lambda r: (not r["isme"], r["year"] or ""))
with open(f"{base}/work/epmc_keep_titles.txt", "w") as f:
    for i, r in enumerate(R):
        r["idx"] = i
        f.write("%d|%s|%s|%s|%s|%s\n" % (i, r["pmid"] or r["doi"], r["year"], "T" if r["trialpt"] else "", "ME" if r["isme"] else "LC", r["title"][:150]))
json.dump(R, open(f"{base}/work/epmc_keep.json", "w"))
