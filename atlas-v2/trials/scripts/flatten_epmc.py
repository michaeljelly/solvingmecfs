import json, glob, sys, re
base = sys.argv[1]
recs = {}
for f in sorted(glob.glob(f"{base}/raw/epmc/*_*.json")):
    tag = f.split("/")[-1].split("_")[0]
    for r in json.load(open(f))["resultList"]["result"]:
        k = r.get("pmid") or r.get("doi") or r["id"]
        pt = r.get("pubTypeList", {}).get("pubType", [])
        ab = r.get("abstractText") or ""
        ab = re.sub(r"<[^>]+>", " ", ab)
        e = recs.setdefault(k, dict(pmid=r.get("pmid"), pmcid=r.get("pmcid"), doi=r.get("doi"), title=re.sub(r"<[^>]+>", "", r.get("title", "")), year=r.get("pubYear"), journal=(r.get("journalInfo", {}).get("journal", {}) or {}).get("title") or r.get("bookOrReportDetails", {}).get("publisher"), authors=r.get("authorString", "")[:200], pubtypes=pt, abstract=ab, source=r.get("source"), tags=[]))
        e["tags"].append(tag)
json.dump(list(recs.values()), open(f"{base}/work/epmc_flat.json", "w"))
print(len(recs))
