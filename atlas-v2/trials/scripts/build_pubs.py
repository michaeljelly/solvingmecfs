import json, sys, re, glob
base = sys.argv[1]
K = json.load(open(f"{base}/work/epmc_keep.json"))
byidx = {r["idx"]: r for r in K}
sel = []; extra = []
for line in open(f"{base}/work/pub_select.txt"):
    if line.startswith("#"): continue
    if line.startswith("EXTRA"): extra = line.split()[1:]; continue
    sel += [int(x) for x in line.split()]
out = []
for i in sel:
    r = byidx[i]; out.append({k: r[k] for k in ("pmid", "pmcid", "doi", "title", "year", "journal", "authors", "pubtypes", "abstract")} | {"isme": r["isme"], "islc": r["islc"]})
for pid in extra:
    d = json.load(open(f"{base}/raw/epmc_extra/{pid}.json"))["resultList"]["result"]
    if not d: print("missing", pid); continue
    r = d[0]
    ab = re.sub(r"<[^>]+>", " ", r.get("abstractText") or "")
    t = re.sub(r"<[^>]+>", "", r.get("title", ""))
    out.append(dict(pmid=r.get("pmid"), pmcid=r.get("pmcid"), doi=r.get("doi"), title=t, year=r.get("pubYear"), journal=(r.get("journalInfo", {}).get("journal", {}) or {}).get("title"), authors=r.get("authorString", "")[:200], pubtypes=r.get("pubTypeList", {}).get("pubType", []), abstract=ab, isme=bool(re.search(r"chronic fatigue|myalgic|ME/CFS", t + ab, re.I)), islc=bool(re.search(r"covid|pasc", t + ab, re.I))))
for n, r in enumerate(out): r["pid"] = "P%04d" % n
json.dump(out, open(f"{base}/work/pubs_selected.json", "w"))
print(len(out), len(sel), len(extra))
