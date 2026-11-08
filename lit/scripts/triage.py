#!/usr/bin/env python3
"""Merge raw Europe PMC pulls into lit/library.json: one row per work (preprint+journal versions merged),
with a type (primary / review / meta / trial / protocol / case / commentary / other), domain tags and whether we already cite it.
usage: python3 -I triage.py"""
import json, os, re, glob, html
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); R = os.path.dirname(H)
def norm(s): return re.sub(r"[^a-z0-9]", "", (s or "").lower())[:120]
rows = {}
for f in sorted(glob.glob(os.path.join(H, "raw", "epmc_*.json"))):
    pop = "ME/CFS" if f.endswith("_me.json") else "long COVID"
    for x in json.load(open(f, encoding="utf8")):
        t = html.unescape(re.sub(r"<[^>]+>", "", x.get("title") or "")).strip()
        k = norm(t)
        if not k: continue
        pt = [p.lower() for p in (x.get("pubTypeList") or {}).get("pubType", [])]
        r = {"id": x.get("id"), "src": x.get("source"), "pmid": x.get("pmid"), "pmcid": x.get("pmcid"), "doi": (x.get("doi") or "").lower() or None,
             "title": t, "authors": (x.get("authorString") or "")[:300], "journal": x.get("journalTitle"), "year": x.get("pubYear"), "date": x.get("firstPublicationDate"),
             "abstract": html.unescape(re.sub(r"<[^>]+>", " ", x.get("abstractText") or "")).strip(), "pubtypes": pt, "lang": x.get("language"),
             "oa": x.get("isOpenAccess") == "Y", "cites": x.get("citedByCount") or 0, "pop": pop,
             "mesh": [m.get("descriptorName") for m in (x.get("meshHeadingList") or {}).get("meshHeading", [])][:25],
             "kw": (x.get("keywordList") or {}).get("keyword", [])[:20]}
        o = rows.get(k)
        if o:  # keep the journal version, remember the preprint
            if o["pop"] != pop: o["pop"] = "both"
            if o["src"] == "PPR" and r["src"] != "PPR": r["preprint_doi"] = o.get("doi"); r["pop"] = o["pop"]; rows[k] = r
            elif r["src"] == "PPR" and o["src"] != "PPR": o["preprint_doi"] = r.get("doi")
            continue
        rows[k] = r
REV = r"\breview\b|systematic-review|meta-analysis"
def kind(r):
    pt, t, a = " ".join(r["pubtypes"]), r["title"].lower(), r["abstract"].lower()
    if re.search(r"erratum|correction|retraction|expression of concern", pt) or re.match(r"^(correction|erratum|retraction)", t): return "other"
    if re.search(r"meta-analysis|systematic", pt) or re.search(r"meta-analys|systematic review|scoping review|umbrella review", t): return "meta"
    if re.search(r"protocol", t) and not re.search(r"result", t): return "protocol"
    if re.search(r"randomi[sz]ed controlled trial|clinical trial", pt) or re.search(r"randomi[sz]ed|placebo-controlled|pilot trial|open-label", t): return "trial"
    if re.search(r"case report", pt) or re.search(r"case report|a case of|case series", t): return "case"
    if re.search(REV, pt) or re.search(r"\breview\b|narrative|perspective|hypothesis|viewpoint|overview|update on|state of the art", t): return "review"
    if re.search(r"editorial|comment|letter|news|interview|biography|historical", pt) and len(a) < 900: return "commentary"
    if len(a) < 200: return "no-abstract"
    return "primary"
DOM = [("vascular", r"vascul|endotheli|microcirc|capillar|blood flow|blood volume|orthostat|POTS|tachycardia|coagul|clot|platelet|complement|hemodynam|haemodynam|cardiac|stroke volume|preload|tilt"),
       ("muscle-exertion", r"muscle|exerti|exercise|CPET|mitochond|bioenerget|lactate|fatigabil|oxygen uptake|VO2|post-exertional"),
       ("neuro", r"brain|neuro|cognit|MRI|PET\b|autonomic|small[- ]fib|microglia|cerebr|sleep|pain|vagus|vagal|dysautonom|EEG"),
       ("immune", r"immun|T[- ]cell|B[- ]cell|\bNK\b|natural killer|cytokine|inflamm|autoantibod|antibod|interferon|mast cell|monocyte|autoimmun"),
       ("infection", r"viral|virus|EBV|Epstein|herpes|HHV|enterovir|persisten|reservoir|reactivat|antiviral|SARS-CoV-2"),
       ("omics", r"metabolom|proteom|lipidom|transcriptom|omic|biomarker|metaboli|single[- ]cell|methylation|epigen|machine learning"),
       ("endocrine", r"hormon|cortisol|HPA|endocrin|thyroid|vasopressin|pituitar|adrenal|estrogen|oestrogen|menopaus|sex differen"),
       ("gut", r"microbio|gut\b|intestin|bowel|gastro|dysbios"), ("genetics", r"genetic|genom|GWAS|heritab|polygenic|variant|twin"),
       ("trials", r"\btrial\b|randomi[sz]ed|placebo|treatment|therap|drug|intervention|rehabilit"),
       ("epi", r"epidemiol|prevalence|incidence|cohort|registry|natural history|diagnos|quality of life|health service|paediatric|pediatric|adolescen|children|mortality|cost")]
cited = ""
for f in glob.glob(os.path.join(R, "atlas", "build", "research", "*.json")) + [os.path.join(R, "site-src", "data", n) for n in ("trials.json", "treatments.json")] + glob.glob(os.path.join(R, "kb", "*.json")):
    if not os.path.exists(f): continue
    t = open(f, encoding="utf8").read()
    if "/research/" in f.replace(os.sep, "/"):  # curated evidence only; the automatic "sweep" lists do not count as cited
        try:
            d = json.loads(t)
            for x in d.get("findings", []): x.pop("sweep", None)
            t = json.dumps(d, ensure_ascii=False)
        except Exception: pass
    cited += t.lower()
PMIDS = set(re.findall(r"\b\d{6,9}\b", cited)); DOIS = set(re.findall(r"10\.\d{4,9}/[^\s\"'<>)\]]+", cited))
out = []
for r in rows.values():
    r["kind"] = kind(r)
    txt = r["title"] + " " + r["abstract"][:3000]
    r["domains"] = [k for k, p in DOM if re.search(p, txt, re.I)][:4]
    r["cited"] = bool((r["doi"] and any(d.startswith(r["doi"]) for d in DOIS if d[:12] == r["doi"][:12])) or (r.get("preprint_doi") and r["preprint_doi"] in DOIS) or (r["pmid"] and r["pmid"] in PMIDS))
    out.append(r)
out.sort(key=lambda r: (r.get("date") or ""), reverse=True)
json.dump(out, open(os.path.join(H, "library.json"), "w", encoding="utf8"), ensure_ascii=False)
from collections import Counter
print(len(out), Counter(r["kind"] for r in out), Counter(r["pop"] for r in out))
p = [r for r in out if r["kind"] in ("primary", "trial", "meta") and not r["cited"]]
print("to screen (uncited primary/trial/meta):", len(p), Counter((r["year"] or "")[:3] + "x" for r in p), Counter(r["pop"] for r in p))
