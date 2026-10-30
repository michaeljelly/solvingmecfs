#!/usr/bin/env python3
"""Check that every id in data/graph.json and data/treatments.json points at something that exists, and flag atlas findings
that duplicate each other with different grades. usage: python3 -I check_graph.py <site_src_dir> <atlas_research_dir>"""
import glob, json, os, sys
SRC, RES = sys.argv[1], sys.argv[2]
G = json.load(open(os.path.join(SRC, "data", "graph.json"), encoding="utf8"))
T = json.load(open(os.path.join(SRC, "data", "treatments.json"), encoding="utf8"))["treatments"]
F, H = {}, {h["id"] for h in json.load(open(os.path.join(RES, "hypotheses.json"), encoding="utf8"))["hypotheses"]}
for p in glob.glob(os.path.join(RES, "*.json")):
    for f in json.load(open(p, encoding="utf8")).get("findings", []):
        F[f["id"]] = f
steps = {s["n"] for s in G["steps"]}
qs, es = {q["id"] for q in G["questions"]}, {e["id"] for e in G["experiments"]}
err = []
def need(kind, ids, pool, where):
    for i in ids:
        if i not in pool: err.append(f"{where}: unknown {kind} {i!r}")
for s in G["steps"]:
    w = f"step {s['n']}"
    for k in ("evidence", "weak", "tension"): need("finding", s.get(k, []), F, w)
    need("hypothesis", s.get("hypotheses", []), H, w)
    if s["found"] not in G["grades"]: err.append(f"{w}: unknown grade {s['found']}")
for e in G["edges"]:
    need("step", [x for x in (e["from"], e["to"]) if x != "infection"], steps, "edge")
for c in G["constraints"]: need("finding", c["findings"], F, "constraint " + c["fact"][:30])
for q in G["questions"]:
    w = "question " + q["id"]
    need("finding", q["findings"], F, w); need("step", q["steps"], steps, w)
    if q.get("experiment"): need("experiment", [q["experiment"]], es, w)
for q in G["small_questions"]: need("finding", [q["finding"]], F, "small question")
for e in G["experiments"]:
    w = "experiment " + e["id"]
    need("step", e["steps"], steps, w); need("question", e["questions"], qs, w); need("hypothesis", e["hypotheses"], H, w)
for t in T:
    w = "treatment " + t["id"]
    need("step", t.get("links") or [], steps, w); need("hypothesis", t.get("hypotheses") or [], H, w)
    if not t.get("model_rationale"): err.append(f"{w}: no model_rationale")
    if t["verdict"] == "try" and t["evidence"] in ("rct-negative", "none"): err.append(f"{w}: 'try' with evidence {t['evidence']}")
# findings that claim the same thing in two domains should carry the same grade
DUPES = [("neuro-cbf-tilt", "vascular-cbf-tilt"), ("neuro-hypocapnia", "vascular-hypocapnia"), ("gut-low-ipa", "omics-ipa"),
         ("muscle-exertion-muscle-proteins-blood", "omics-proteome-muscle-low"), ("immune-rituximab", "trials-rituximab"),
         ("immune-immunoadsorption", "trials-immunoadsorption"), ("immune-bc007", "trials-bc007"), ("immune-daratumumab", "trials-daratumumab"),
         ("infection-paxlovid-lc-null", "trials-paxlovid"), ("immune-cfrna", "omics-cfrna"), ("endocrine-renin-aldosterone-paradox", "vascular-raas-paradox"),
         ("muscle-exertion-oxygen-extraction", "vascular-oxygen-extraction")]
for a, b in DUPES:
    if a in F and b in F and F[a]["grade"] != F[b]["grade"]:
        err.append(f"grade mismatch: {a} {F[a]['grade']} vs {b} {F[b]['grade']}")
for f in F.values():
    blob = json.dumps(f)
    for bad in ("I found", "from memory"):
        if bad in blob: err.append(f"finding {f['id']}: contains {bad!r}")
print("\n".join(err) if err else "graph ok")
sys.exit(1 if err else 0)
