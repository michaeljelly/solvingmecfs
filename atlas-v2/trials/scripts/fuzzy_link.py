import json, re, sys
base = sys.argv[1]
P = json.load(open(f"{base}/work/pubs_selected.json"))
S = {x["nct"]: x for x in json.load(open(f"{base}/work/screen.json"))}
CT = {}
for f in ("ct_mecfs_flat", "ct_lc_flat", "ct_mecfsterm_flat", "ct_lcterm_flat"):
    for r in json.load(open(f"{base}/work/{f}.json")): CT[r["nct"]] = r
STOP = set("a an the of in on for with and or to by at as from vs versus patients patient study trial randomized randomised controlled placebo double blind pilot effect effects efficacy safety treatment syndrome chronic fatigue long covid post covid 19 sars cov 2 condition people adults individuals myalgic encephalomyelitis cfs me clinical among who after acute sequelae infection".split())
def toks(t): return {w for w in re.sub(r"[^a-z0-9 ]", " ", (t or "").lower()).split() if w not in STOP and len(w) > 2}
cands = []
EXTRA_IN = set(json.load(open(f'{base}/work/extra_ct.json'))) if __import__('os').path.exists(f'{base}/work/extra_ct.json') else set()
for p in P:
    if p["nct_links"]: continue
    pt = toks(p["title"])
    if not pt: continue
    best = []
    for n, r in CT.items():
        if S.get(n, {}).get("decision") != "in" and n not in EXTRA_IN: continue
        sy = int((r["start"] or "1900")[:4]); py = int(p["year"] or 2100)
        if py < sy: continue
        ct = toks(r["title"]) | toks(r["brief"]) | toks(" ".join(i[1] or "" for i in r["interventions"]))
        j = len(pt & ct) / max(1, len(pt))
        if j >= 0.6 and len(pt & ct) >= 2: best.append((round(j, 2), n))
    best.sort(reverse=True)
    if best: cands.append((p["pid"], p["year"], p["title"][:100], [(b, n, CT[n]["brief"][:70], (CT[n]["start"] or "")[:4]) for b, n in best[:3]]))
for c in cands:
    print(c[0], c[1], c[2]); [print("    ", x) for x in c[3]]
