#!/usr/bin/env python3
"""Fetch ME/CFS and long COVID research projects from public sources into register/raw/.
Sources: ME/CFS Research Register (MRR, Hamburg), NIH RePORTER, UKRI Gateway to Research.
usage: python3 fetch_register.py"""
import json, os, subprocess, sys, time
H = os.path.dirname(os.path.abspath(__file__)); RAW = os.path.join(H, "raw")
sys.path.insert(0, H)
import devalue

def curl(args):
    for i in range(4):
        r = subprocess.run(["curl", "-sL", "-m", "60"] + args, capture_output=True, text=True)
        if r.returncode == 0 and r.stdout: return r.stdout
        time.sleep(2 ** i)
    raise RuntimeError("fetch failed: " + " ".join(args))

def mrr():
    out = {}
    for t in ["projects", "organisations", "persons", "researchnetworks", "workinggroups"]:
        txt = curl([f"https://mrr.mecfs-research.org/en/mrr/{t}/__data.json"])
        p = json.loads(txt); props = devalue.decode([n for n in p["nodes"] if n and n.get("type") == "data"][-1]["data"])["props"]
        out[t] = next(v for v in props.values() if isinstance(v, list) and len(v) > 3)
    json.dump(out, open(os.path.join(RAW, "mrr.json"), "w"), ensure_ascii=False, default=str)
    print("mrr", {k: len(v) for k, v in out.items()})

NIH_Q = {"mecfs": '"myalgic encephalomyelitis" "chronic fatigue syndrome" "ME/CFS"',
         "lc": '"long COVID" "post-acute sequelae" "post-COVID condition" "post-COVID-19 condition" PASC "long-COVID"'}
def nih():
    res = []
    for pop, q in NIH_Q.items():
        off = 0
        while True:
            body = {"criteria": {"advanced_text_search": {"operator": "or", "search_field": "projecttitle,terms,abstracttext", "search_text": q},
                                 "fiscal_years": list(range(2016, 2027))},
                    "include_fields": ["ApplId", "FiscalYear", "ProjectNum", "CoreProjectNum", "ProjectTitle", "AbstractText", "Organization", "PrincipalInvestigators",
                                       "AgencyIcAdmin", "AwardAmount", "ProjectStartDate", "ProjectEndDate", "IsActive", "ActivityCode", "ProjectDetailUrl", "SubprojectId", "Terms"],
                    "offset": off, "limit": 500}
            d = json.loads(curl(["-X", "POST", "https://api.reporter.nih.gov/v2/projects/search", "-H", "Content-Type: application/json", "-d", json.dumps(body)]))
            for r in d["results"]: r["_query"] = pop
            res += d["results"]; off += 500
            if off >= d["meta"]["total"] or off >= 14500: break
            time.sleep(1)
    json.dump(res, open(os.path.join(RAW, "nih.json"), "w"))
    print("nih", len(res))

def gtr():
    import re as _re
    ME = r"myalgic|chronic fatigue syndrome|ME/CFS|CFS/ME|long[- ]covid|post-acute sequelae|post-covid(-19)? (syndrome|condition)|PASC"
    seen, res = set(), []
    for q in ['"myalgic encephalomyelitis"', '"chronic fatigue syndrome"', '"ME/CFS"', '"long covid"', '"post-acute"', '"post covid syndrome"', '"post-covid condition"']:
        for f in ("pro.t", "pro.a"):
            p = 1
            while True:
                d = json.loads(curl(["-H", "Accept: application/vnd.rcuk.gtr.json-v7", f"https://gtr.ukri.org/gtr/api/projects?q={q.replace(' ', '%20').replace(chr(34), '%22')}&f={f}&s=100&p={p}"]))
                for x in d.get("project", []):
                    if x["id"] in seen or not _re.search(ME, (x.get("title") or "") + " " + (x.get("abstractText") or ""), _re.I): continue
                    seen.add(x["id"])
                    for l in x.get("links", {}).get("link", []):
                        if l["rel"] in ("LEAD_ORG", "PI_PER", "FUND"):
                            try: l["_data"] = json.loads(curl(["-H", "Accept: application/vnd.rcuk.gtr.json-v7", l["href"].replace("http://", "https://")]))
                            except Exception: pass
                    x["_query"] = q; res.append(x)
                if p >= d.get("totalPages", 1): break
                p += 1
    json.dump(res, open(os.path.join(RAW, "gtr.json"), "w"))
    print("gtr", len(res))

def cordis():
    import urllib.parse
    q = "contenttype='project' AND ('myalgic' OR 'chronic fatigue syndrome' OR 'ME/CFS' OR 'long covid' OR 'post-covid' OR 'post-acute sequelae' OR 'post covid-19 condition')"
    res, pg = [], 1
    while True:
        d = json.loads(curl([f"https://cordis.europa.eu/search/en?q={urllib.parse.quote(q)}&format=json&p={pg}&num=10"]))
        hits = (d.get("hits") or {}).get("hit") or []
        if isinstance(hits, dict): hits = [hits]
        res += [h["project"] for h in hits if "project" in h]
        if not hits or pg * 10 >= int(d["result"]["header"]["totalHits"]): break
        pg += 1
    json.dump(res, open(os.path.join(RAW, "cordis.json"), "w"), ensure_ascii=False)
    print("cordis", len(res))

def polybio():
    import re, html as H
    res = []
    for pg in range(1, 30):
        x = curl([f"https://polybio.org/projects/feed/?paged={pg}"])
        items = re.findall(r"(?s)<item>(.*?)</item>", x)
        if not items: break
        for it in items:
            g = lambda t: H.unescape(re.sub(r"<!\[CDATA\[|\]\]>", "", (re.search(rf"(?s)<{t}>(.*?)</{t}>", it) or [None, ""])[1])).strip()
            res.append({"title": g("title"), "url": g("link"), "date": g("pubDate"), "creator": g("dc:creator"), "description": re.sub(r"<[^>]+>", " ", g("description")),
                        "content": re.sub(r"<[^>]+>", " ", g("content:encoded"))[:4000]})
        if len(items) < 10: break
    json.dump(res, open(os.path.join(RAW, "polybio.json"), "w"), ensure_ascii=False)
    print("polybio", len(res))

def cihr():
    """CIHR Funding Decisions Database (Solr behind webapps.cihr-irsc.gc.ca/decisions). No phrase search on acallText,
    so query broad AND terms and filter on title/abstract/keywords in build_register."""
    import urllib.parse
    qs = ["myalgic", "encephalomyelitis", "(chronic AND fatigue AND syndrome)", "(long AND covid)", "(post AND covid AND condition)",
          "PASC", "post-covid", "post-exertional", "ME\\/CFS", "(post-acute AND sequelae)"]
    fl = "id,projecttitle,abstract,keyworddelim,pinamesdelim,orgnameinp2,cihramount2,competitiondate,approvedterm2,competitioncd,programname2,prcname2,instname2"
    res = {}
    for q in qs:
        start = 0
        while True:
            u = "https://webapps.cihr-irsc.gc.ca/decisions/sq?" + urllib.parse.urlencode({"q": "acallText:" + q, "wt": "json", "rows": 500, "start": start, "fl": fl})
            d = json.loads(curl([u]))["response"]
            for x in d["docs"]: res[x["id"]] = x
            start += 500
            if start >= d["numFound"]: break
        print("cihr", q, d["numFound"], len(res))
    json.dump(list(res.values()), open(os.path.join(RAW, "cihr.json"), "w"), ensure_ascii=False)

if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    for f in (sys.argv[1:] or ["mrr", "nih", "gtr", "cordis", "polybio", "cihr"]): globals()[f]()
