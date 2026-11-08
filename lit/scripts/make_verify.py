#!/usr/bin/env python3
"""Build verification batches: one item per (paper, atlas finding) link that integrate.py would put in a sweep list,
plus relation=new candidates. Each item carries the full abstract and the full finding text.
usage: python3 -I make_verify.py [size]  -> lit/verify/v-NN.jsonl"""
import json, os, glob, sys
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.dirname(H)
size = int(sys.argv[1]) if len(sys.argv) > 1 else 100
L = {str(r["id"]): r for r in json.load(open(os.path.join(H, "library.json"), encoding="utf8"))}
F = {}
for p in glob.glob(os.path.join(R, "atlas", "build", "research", "*.json")):
    if os.path.basename(p) in ("overview.json", "hypotheses.json"): continue
    for f in json.load(open(p, encoding="utf8"))["findings"]: F[f["id"]] = f
done = set()
for f in glob.glob(os.path.join(H, "verify", "out-*.json")):
    done |= {x.get("key") for x in json.load(open(f, encoding="utf8")) if isinstance(x, dict)}
items = []
for f in sorted(glob.glob(os.path.join(H, "screen", "*.json"))):
    for x in json.load(open(f, encoding="utf8")):
        if not isinstance(x, dict) or str(x.get("id")) not in L: continue
        if x.get("relevant") not in ("core", "related") or (x.get("importance") or 1) < 2: continue
        rel = x.get("relation"); r = L[str(x["id"])]
        paper = {"title": r["title"], "year": r.get("year"), "journal": r.get("journal"), "pop": r.get("pop"), "abstract": r.get("abstract") or ""}
        if rel in ("replicates", "contradicts", "extends"):
            for fid in (x.get("atlas") or [])[:2]:
                if fid in F:
                    key = f"{x['id']}|{fid}"
                    if key in done: continue
                    fd = F[fid]
                    items.append({"key": key, "paper_id": str(x["id"]), "proposed_relation": rel, "screen_summary": x.get("finding"), **paper,
                                  "finding_id": fid, "finding_grade": fd.get("grade"), "finding_claim": fd["finding"], "finding_detail": fd.get("detail", "")[:2500],
                                  "finding_studies": [s.get("cite") for s in (fd.get("studies") or []) + (fd.get("contradicting") or [])]})
        elif rel == "new":
            key = f"{x['id']}|new"
            if key not in done:
                items.append({"key": key, "paper_id": str(x["id"]), "proposed_relation": "new", "screen_summary": x.get("finding"), **paper})
os.makedirs(os.path.join(H, "verify"), exist_ok=True)
for i in range(0, len(items), size):
    with open(os.path.join(H, "verify", f"v-{i // size + 1:02d}.jsonl"), "w", encoding="utf8") as o:
        for it in items[i:i + size]: o.write(json.dumps(it, ensure_ascii=False) + "\n")
print(len(items), "items in", -(-len(items) // size), "batches")
