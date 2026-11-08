#!/usr/bin/env python3
"""Fold screened papers into the atlas research JSONs and write review lists.
- importance >= 2 papers that replicate/extend/contradict an atlas finding go into that finding's "sweep" list
  (kept apart from the curated "studies"/"contradicting", which are only changed by hand after review).
- writes lit/grade_review.md (findings with new evidence, for manual re-grading) and lit/new_candidates.md (new findings).
usage: python3 -I integrate.py [--apply]"""
import json, os, re, glob, sys
from collections import defaultdict
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.dirname(H)
RES = os.path.join(R, "atlas", "build", "research")
APPLY = "--apply" in sys.argv
L = {str(r["id"]): r for r in json.load(open(os.path.join(H, "library.json"), encoding="utf8"))}
S = []
for f in sorted(glob.glob(os.path.join(H, "screen", "*.json"))):
    for x in json.load(open(f, encoding="utf8")):
        if isinstance(x, dict) and str(x.get("id")) in L: S.append(x)
def url(r): return ("https://doi.org/" + r["doi"]) if r.get("doi") else (f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/" if r.get("pmid") else f"https://europepmc.org/article/{r['src']}/{r['id']}")
def cite(r):
    a = (r.get("authors") or "").split(",")[0].strip()
    a = re.sub(r"\s+[A-Z]{1,6}$", "", a)
    a = re.sub(r"^([A-Z][a-z]+) [A-Z]{1,6} (van|von|de|der|den)$", r"\2 \1", a)
    a = re.sub(r"\s+\S*\(\S*$", "", a)
    if a == "Campen": a = "van Campen"
    return f"{a} {r.get('year')}, {r.get('journal') or ('preprint' if r.get('src') == 'PPR' else '')}".strip(", ")
# collapse preprint/journal duplicates: same first author surname + same first 50 normalised title chars
seen, U = {}, []
for x in sorted(S, key=lambda x: L[str(x["id"])].get("src") == "PPR"):
    r = L[str(x["id"])]
    k = (re.sub(r"\W", "", (r.get("authors") or "").split(" ")[0].lower()), re.sub(r"[^a-z]", "", r["title"].lower())[:50])
    if k in seen: continue
    seen[k] = 1; U.append(x)
files = {}
for p in glob.glob(os.path.join(RES, "*.json")):
    if os.path.basename(p) in ("overview.json", "hypotheses.json"): continue
    files[p] = json.load(open(p, encoding="utf8"))
F = {f["id"]: (p, f) for p, d in files.items() for f in d["findings"]}
add = defaultdict(lambda: {"rep": [], "con": [], "ext": []})
newc = defaultdict(list)
# second-pass verdicts (VERIFY_BRIEF.md): once any exist, only verified links are used
V = {}
for vf in glob.glob(os.path.join(H, "verify", "out-*.json")):
    for v in json.load(open(vf, encoding="utf8")):
        if isinstance(v, dict) and v.get("key"): V[v["key"]] = v
VERIFIED = bool(V)
RK = {"replicates": "rep", "contradicts": "con", "extends": "ext"}
def link(fid, rel, imp, r, x, v=None):
    if fid not in F or rel not in RK or imp < 2: return
    y = dict(x)
    if v: y.update({"finding": v.get("result") or x.get("finding"), "verified2": True})
    add[fid][RK[rel]].append((imp, r, y))
for x in U:
    r = L[str(x["id"])]
    if x.get("relevant") not in ("core", "related"): continue
    imp = x.get("importance") or 1
    rel = x.get("relation")
    if not VERIFIED:
        if rel == "new" and imp >= 2: newc[(r.get("domains") or ["other"])[0]].append((imp, r, x))
        if rel in RK and imp >= 2:
            for fid in (x.get("atlas") or [])[:2]: link(fid, rel, imp, r, x)
        continue
    keys = [f"{x['id']}|{fid}" for fid in (x.get("atlas") or [])[:2]] if rel in RK else [f"{x['id']}|new"] if rel == "new" else []
    for k in keys:
        v = V.get(k)
        if not v or v.get("verdict") == "drop": continue
        vi = v.get("importance") or imp
        if v.get("verdict") == "new_ok":
            if vi >= 2: newc[(r.get("domains") or ["other"])[0]].append((vi, r, dict(x, finding=v.get("result") or x.get("finding"))))
            continue
        for fid in (v.get("atlas") or [])[:2]: link(fid, v.get("relation"), vi, r, x, v)
for d in add.values():  # one entry per paper per finding
    for kk in d:
        seenp, keep = set(), []
        for t in d[kk]:
            if t[1]["id"] in seenp: continue
            seenp.add(t[1]["id"]); keep.append(t)
        d[kk] = keep
n_add = 0
for fid0, (p0, f0) in F.items(): f0.pop("sweep", None)
for fid, d in add.items():
    p, f = F[fid]
    f["sweep"] = []  # rebuilt on every run; curated "studies"/"contradicting" are never touched
    have = {(s.get("url") or "").lower() for s in (f.get("studies") or []) + (f.get("contradicting") or [])}
    for key, tgt in (("rep", "sweep"), ("ext", "sweep"), ("con", "sweep")):
        for imp, r, x in sorted(d[key], key=lambda t: -t[0]):
            u = url(r)
            if u.lower() in have: continue
            f.setdefault(tgt, []).append({"cite": cite(r), "n": x.get("n") or "", "design": x.get("design") or "", "result": x.get("finding") or "", "url": u,
                                          "verified": True, "sweep": True, "checked": bool(x.get("verified2")), "relation": {"rep": "replicates", "ext": "extends", "con": "contradicts"}[key], "importance": imp})
            have.add(u.lower()); n_add += 1
if APPLY:
    for p, d in files.items():
        json.dump(d, open(p, "w", encoding="utf8"), ensure_ascii=False, indent=1)
# review lists
out = ["# Findings with new evidence from the literature sweep", "", "Re-grade by hand. rep/ext/con = new replicating, extending, contradicting studies (importance >= 2).", ""]
for fid, d in sorted(add.items(), key=lambda kv: -(len(kv[1]["con"]) * 3 + len(kv[1]["rep"]) * 2 + len(kv[1]["ext"]))):
    p, f = F[fid]
    out.append(f"## {fid} [{f['grade']}] rep {len(d['rep'])} ext {len(d['ext'])} con {len(d['con'])}")
    out.append(f"Claim: {f['finding']}")
    for key in ("con", "rep", "ext"):
        for imp, r, x in sorted(d[key], key=lambda t: -t[0])[:6]:
            out.append(f"- {key} i{imp} {cite(r)} ({x.get('n') or '?'}): {x.get('finding')} {url(r)}")
    out.append("")
open(os.path.join(H, "grade_review.md"), "w", encoding="utf8").write("\n".join(out))
o2 = ["# Candidate new findings (relation=new, importance >= 2)", ""]
for dom, v in sorted(newc.items()):
    o2.append(f"## {dom} ({len(v)})")
    for imp, r, x in sorted(v, key=lambda t: -t[0]):
        o2.append(f"- i{imp} {cite(r)} [{r['pop']}] ({x.get('n') or '?'}, {x.get('design')}): {x.get('finding')} {url(r)}")
    o2.append("")
o2 += ["", "# Verified results with no matching atlas finding (importance >= 2)", ""]
seen_h = set()
for v in sorted(V.values(), key=lambda v: -(v.get("importance") or 1)):
    rs = (v.get("reason") or "").lower()
    if v.get("verdict") != "drop" or (v.get("importance") or 1) < 2 or not ("no atlas finding" in rs or "could be" in rs): continue
    pid = v["key"].split("|")[0]
    if pid in seen_h or pid not in L: continue
    seen_h.add(pid); r = L[pid]
    o2.append(f"- i{v.get('importance')} {cite(r)} [{r['pop']}]: {v.get('result')} {url(r)}")
open(os.path.join(H, "new_candidates.md"), "w", encoding="utf8").write("\n".join(o2))
print("screened", len(S), "unique", len(U), "findings touched", len(add), "entries added", n_add, "new candidates", sum(len(v) for v in newc.values()), "applied" if APPLY else "(dry run)")
