import json, sys, re
base = sys.argv[1]
P = json.load(open(f"{base}/work/pubs_selected.json"))
S = {x["nct"]: x for x in json.load(open(f"{base}/work/screen.json"))}
CT = {}
for f in ("ct_mecfs_flat", "ct_lc_flat", "ct_mecfsterm_flat", "ct_lcterm_flat"):
    for r in json.load(open(f"{base}/work/{f}.json")): CT[r["nct"]] = r
extra = set(json.load(open(f"{base}/work/extra_ct.json")))
manual = json.load(open(f"{base}/work/manual_links.json"))
for p in P:
    if p["pid"] in manual: p["nct_links"] = sorted(set(p["nct_links"]) | {manual[p["pid"]]})
inc = {n for n, x in S.items() if x["decision"] == "in"} | extra
trials = {n: {"nct": n, "pubs": []} for n in inc}
flag = []
drop = json.load(open(f"{base}/work/pub_drop.json"))
for p in P:
    if p["pid"] in drop: continue
    ln = [n for n in p["nct_links"] if n in CT]
    if not ln:
        trials["PUB:" + p["pid"]] = {"nct": None, "pubs": [p["pid"]]}; continue
    for n in ln[:1]:
        if n not in trials:
            trials[n] = {"nct": n, "pubs": [], "forced": S.get(n, {}).get("reason", "not screened")}
            flag.append((n, p["pid"], S.get(n, {}).get("reason"), p["title"][:80]))
        trials[n]["pubs"].append(p["pid"])
json.dump(trials, open(f"{base}/work/trials_skeleton.json", "w"))
print(len(trials), sum(1 for t in trials.values() if t["nct"]), sum(1 for t in trials.values() if not t["nct"]))
for f in flag: print("FORCED", f)
