#!/usr/bin/env python3
"""Pull the ME/CFS literature (all years) and the PEM/ME-relevant long COVID literature from Europe PMC
(covers PubMed, PMC and preprints from medRxiv/bioRxiv/Research Square/SSRN). Writes lit/raw/epmc_<set>.json.
usage: python3 -I fetch_lit.py [me] [lc]"""
import json, os, subprocess, sys, time, urllib.parse
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RAW = os.path.join(H, "raw")
ME = ('(TITLE:"chronic fatigue syndrome" OR TITLE:"myalgic encephalomyelitis" OR TITLE:"ME/CFS" OR TITLE:"CFS/ME" OR TITLE:"myalgic encephalopathy" '
      'OR TITLE:"systemic exertion intolerance" OR ABSTRACT:"myalgic encephalomyelitis" OR ABSTRACT:"ME/CFS")')
LC = ('(TITLE:"long covid" OR TITLE:"post-covid" OR TITLE:"post covid" OR TITLE:"post-acute sequelae" OR TITLE:"PASC" OR TITLE:"post-COVID-19 condition" '
      'OR TITLE:"long-COVID" OR TITLE:"post-acute COVID") AND (ABSTRACT:"post-exertional" OR ABSTRACT:"ME/CFS" OR ABSTRACT:"myalgic" OR TITLE:fatigue '
      'OR TITLE:autonomic OR TITLE:POTS OR TITLE:orthostatic OR TITLE:exercise OR TITLE:exertion OR TITLE:cognitive OR TITLE:"brain fog" OR TITLE:dysautonomia)')
KEEP = ("id", "source", "pmid", "pmcid", "doi", "title", "authorString", "journalTitle", "pubYear", "firstPublicationDate", "abstractText", "pubTypeList",
        "language", "isOpenAccess", "citedByCount", "keywordList", "meshHeadingList", "fullTextUrlList", "commentCorrectionList")
def pull(name, q):
    out, cm = [], "*"
    while True:
        u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 100, "cursorMark": cm, "resultType": "core"})
        for i in range(7):
            r = subprocess.run(["curl", "-sS", "-m", "180", u], capture_output=True, text=True)
            try:
                d = json.loads(r.stdout)
                if "resultList" in d: break
            except Exception: pass
            time.sleep(2 ** i)
        res = d["resultList"]["result"]
        out += [{k: x.get(k) for k in KEEP if x.get(k) is not None} for x in res]
        print(name, len(out), "/", d["hitCount"], flush=True)
        if not res or d.get("nextCursorMark") in (None, cm): break
        cm = d["nextCursorMark"]
    json.dump(out, open(os.path.join(RAW, f"epmc_{name}.json"), "w"), ensure_ascii=False)
os.makedirs(RAW, exist_ok=True)
for s in (sys.argv[1:] or ["me", "lc"]): pull(s, ME if s == "me" else LC)
