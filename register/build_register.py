#!/usr/bin/env python3
"""Merge register/raw/{mrr,nih,gtr}.json into site-src/data/projects.json (one row per project).
usage: python3 -I build_register.py <site_src_dir>"""
import json, os, re, sys, datetime
from collections import defaultdict
H = os.path.dirname(os.path.abspath(__file__)); RAW = os.path.join(H, "raw"); SRC = sys.argv[1]
TODAY = datetime.date.today().isoformat()
ME = r"myalgic|chronic fatigue syndrome|\bME/CFS\b|\bCFS/ME\b|\bSEID\b|post-?exertional malaise"
LC = r"long[- ]covid|post-acute sequelae|\bPASC\b|post-covid|post covid|long-term effects of covid|post-acute covid|post-covid-19 condition"
TOPICS = [  # atlas domain id, label, pattern
    ("vascular", "Blood flow and vessels", r"vascul|endotheli|microcirc|capillar|blood flow|blood volume|orthostat|POTS|tachycardia|coagul|clot|platelet|complement|hemodynam|haemodynam|cardiac|preload"),
    ("muscle-exertion", "Muscle, exertion and energy", r"muscle|exerti|exercise|CPET|mitochond|bioenerget|metabolic rate|lactate|fatigabil"),
    ("neuro", "Brain and nerves", r"brain|neuro|cognit|MRI|PET\b|autonomic|small[- ]fib|microglia|cerebr|sleep|pain|vagus|vagal|dysautonom"),
    ("immune", "Immune system", r"immun|T[- ]cell|B[- ]cell|\bNK\b|natural killer|cytokine|inflamm|autoantibod|antibod|interferon|mast cell|monocyte|autoimmun"),
    ("infection", "Infection and viruses", r"viral|virus|EBV|Epstein|herpes|HHV|enterovir|persisten|reservoir|reactivat|antiviral|SARS-CoV-2 antigen|retrovir"),
    ("omics", "Metabolism and blood omics", r"metabolom|proteom|lipidom|transcriptom|omic|biomarker|metaboli|single[- ]cell|methylation|epigen"),
    ("endocrine", "Hormones", r"hormon|cortisol|HPA|endocrin|thyroid|vasopressin|pituitar|adrenal|estrogen|oestrogen|menopaus|sex differen"),
    ("gut", "Gut and microbiome", r"microbio|gut\b|intestin|bowel|gastro|dysbios"),
    ("genetics", "Genetics", r"genetic|genom|GWAS|heritab|polygenic|variant|sequencing"),
    ("trials", "Treatment trials", r"\btrial\b|randomi[sz]ed|placebo|treatment|therap|drug|repurpos|intervention|rehabilit"),
    ("epi", "Epidemiology and care", r"epidemiol|prevalence|incidence|cohort|registry|biobank|natural history|diagnos|quality of life|health service|care\b|clinic|paediatric|pediatric|adolescen|children"),
]
def topics(title, more=""):
    """Topics named in the title (plus curated research areas); else the two topics the abstract mentions most."""
    t = [k for k, _, p in TOPICS if re.search(p, title, re.I)]
    if t: return t
    sc = sorted(((len(re.findall(p, more, re.I)), k) for k, _, p in TOPICS), reverse=True)
    return [k for c, k in sc[:2] if c >= 2] or ["other"]
def pop(text):
    m, l = bool(re.search(ME, text, re.I)), bool(re.search(LC, text, re.I))
    return "both" if m and l else ("ME/CFS" if m else ("long COVID" if l else None))
def d10(s): return (str(s or "")[:10]) or None
def status(end, given=None):
    if given: return "active" if given.lower().startswith(("ongoing", "active", "recruit")) else "completed"
    return "active" if end and end >= TODAY else ("completed" if end else "unknown")

rows = []
# --- ME/CFS Research Register (Hamburg)
p = os.path.join(RAW, "mrr.json")
if os.path.exists(p):
    M = json.load(open(p, encoding="utf8"))
    for x in M["projects"]:
        text = (x.get("project_title") or "") + " " + (x.get("project_description") or "")
        areas = [a.get("display_value") for a in x.get("research_area_list") or [] if isinstance(a, dict)]
        rows.append({"id": "mrr-" + x["_id"], "source": "MRR", "source_label": "ME/CFS Research Register",
            "url": x.get("url") or x.get("url_other") or "https://mrr.mecfs-research.org/en/mrr/projects",
            "title": x.get("project_title"), "population": pop(text) or "ME/CFS", "focus": "primary",
            "pi": [d.get("display_value") for d in x.get("principal_investigator") or [] if isinstance(d, dict)],
            "org": "; ".join(d.get("display_value") for d in x.get("host_organisations") or [] if isinstance(d, dict)),
            "country": "; ".join(x.get("country") or []), "funder": "; ".join(d.get("display_value", "") for d in (x.get("funding_grant") or []) + (x.get("sponsor") or []) if isinstance(d, dict)),
            "start": d10(x.get("date_start")), "end": d10(x.get("date_end")), "status": status(d10(x.get("date_end")), x.get("project_status")),
            "amount_usd": None, "kind": "; ".join(x.get("research_type_list") or []), "n": x.get("patients_enrolled"),
            "areas": areas, "topics": topics((x.get("project_title") or "") + " " + " ".join(areas), x.get("project_description") or ""), "trial_url": x.get("url_clinical_trial")})
# --- NIH RePORTER: group fiscal years by core project number
p = os.path.join(RAW, "nih.json")
if os.path.exists(p):
    by = defaultdict(list)
    for r in json.load(open(p, encoding="utf8")):
        if r.get("subproject_id"): key = r["core_project_num"] + "-" + str(r["subproject_id"])
        else: key = r["core_project_num"]
        by[key].append(r)
    for key, rs in by.items():
        rs.sort(key=lambda r: r.get("fiscal_year") or 0)
        last = rs[-1]
        title, ab = last.get("project_title") or "", (last.get("abstract_text") or "")[:4000]
        pt = pop(title)
        pa = pop(ab)
        if not (pt or pa): continue  # matched only on indexing terms
        seen_fy, amt = set(), 0
        for r in rs:
            if r["fiscal_year"] not in seen_fy or True: amt += r.get("award_amount") or 0
            seen_fy.add(r["fiscal_year"])
        org = last.get("organization") or {}
        rows.append({"id": "nih-" + key, "source": "NIH", "source_label": "NIH RePORTER",
            "url": f"https://reporter.nih.gov/project-details/{last['appl_id']}",
            "title": title.strip(), "population": pt or pa, "focus": "primary" if pt else "related",
            "pi": [q.get("full_name", "").replace("  ", " ").strip().title() for q in last.get("principal_investigators") or []],
            "org": (org.get("org_name") or "").title(), "country": (org.get("org_country") or "").title().replace("United States", "USA"),
            "funder": "NIH " + ((last.get("agency_ic_admin") or {}).get("abbreviation") or ""),
            "start": d10(last.get("project_start_date")), "end": d10(last.get("project_end_date")), "status": status(d10(last.get("project_end_date"))),
            "amount_usd": amt or None, "kind": last.get("activity_code"), "n": None, "areas": [],
            "topics": topics(title, ab), "years": sorted({r["fiscal_year"] for r in rs})})
# --- UKRI Gateway to Research
p = os.path.join(RAW, "gtr.json")
if os.path.exists(p):
    for x in json.load(open(p, encoding="utf8")):
        title, ab = x.get("title") or "", x.get("abstractText") or ""
        pt, pa = pop(title), pop(ab[:4000])
        if not (pt or pa): continue
        org, pis, amt, st, en = "", [], None, None, None
        for l in x.get("links", {}).get("link", []):
            dd = l.get("_data") or {}
            if l["rel"] == "LEAD_ORG": org = dd.get("name", "")
            if l["rel"] == "PI_PER": pis.append(" ".join(filter(None, [dd.get("firstName"), dd.get("surname")])))
            if l["rel"] == "FUND":
                v = (dd.get("valuePounds") or {}).get("amount") if isinstance(dd.get("valuePounds"), dict) else dd.get("valuePounds")
                amt = round(float(v) * 1.27) if v else None
                st = datetime.datetime.utcfromtimestamp(l["start"] / 1000).date().isoformat() if l.get("start") else None
                en = datetime.datetime.utcfromtimestamp(l["end"] / 1000).date().isoformat() if l.get("end") else None
        rows.append({"id": "gtr-" + x["id"], "source": "UKRI", "source_label": "UKRI Gateway to Research",
            "url": f"https://gtr.ukri.org/projects?ref={x.get('grantReference', '')}", "title": title.strip(), "population": pt or pa,
            "focus": "primary" if pt else "related", "pi": pis, "org": org, "country": "UK", "funder": x.get("leadFunder") or "UKRI",
            "start": st, "end": en, "status": status(en), "amount_usd": amt, "kind": x.get("grantCategory"), "n": None, "areas": [],
            "topics": topics(title, ab)})
# --- EU CORDIS
CC = {"AT": "Austria", "BE": "Belgium", "CH": "Switzerland", "CZ": "Czechia", "DE": "Germany", "DK": "Denmark", "EE": "Estonia", "EL": "Greece", "ES": "Spain", "FI": "Finland",
      "FR": "France", "HU": "Hungary", "IE": "Ireland", "IL": "Israel", "IT": "Italy", "LV": "Latvia", "LT": "Lithuania", "NL": "Netherlands", "NO": "Norway", "PL": "Poland",
      "PT": "Portugal", "SE": "Sweden", "SI": "Slovenia", "SK": "Slovakia", "UK": "UK", "TR": "Turkey", "RO": "Romania", "HR": "Croatia", "BG": "Bulgaria", "CY": "Cyprus", "LU": "Luxembourg", "IS": "Iceland"}
p = os.path.join(RAW, "cordis.json")
if os.path.exists(p):
    for x in json.load(open(p, encoding="utf8")):
        title, ab = x.get("title") or "", (x.get("objective") or "")[:4000]
        pt, pa = pop(title), pop(ab)
        if not (pt or pa): continue
        o = (x.get("relations") or {}).get("associations", {}).get("organization") or []
        o = o if isinstance(o, list) else [o]
        co = next((y for y in o if (y.get("@attributes") or {}).get("type") == "coordinator"), o[0] if o else {})
        amt = x.get("ecMaxContribution")
        rows.append({"id": "cordis-" + str(x["id"]), "source": "CORDIS", "source_label": "EU CORDIS", "url": f"https://cordis.europa.eu/project/id/{x['id']}",
            "title": (x.get("acronym") + ": " if x.get("acronym") else "") + title.strip(), "population": pt or pa, "focus": "primary" if pt else "related",
            "pi": [], "org": (co.get("legalName") or "").title(), "country": CC.get((co.get("address") or {}).get("country"), (co.get("address") or {}).get("country") or ""),
            "funder": "European Commission", "start": d10(x.get("startDate")), "end": d10(x.get("endDate")), "status": status(d10(x.get("endDate"))),
            "amount_usd": round(float(amt) * 1.08) if amt else None, "kind": "EU grant", "n": None, "areas": [], "topics": topics(title, ab)})
# --- PolyBio Research Foundation project archive
p = os.path.join(RAW, "polybio.json")
if os.path.exists(p):
    for x in json.load(open(p, encoding="utf8")):
        text = x["title"] + " " + x["description"] + " " + x["content"]
        pp = pop(text)
        if not pp: continue  # Lyme- or Alzheimer's-only projects
        dt = datetime.datetime.strptime(x["date"][5:16], "%d %b %Y").date().isoformat() if x.get("date") else None
        m = re.search(r"(?:Principal Investigator|Lead Investigator|PI)s?:?\s*([^\n.;]{3,80})", x["content"])
        rows.append({"id": "polybio-" + x["url"].rstrip("/").rsplit("/", 1)[-1], "source": "PolyBio", "source_label": "PolyBio Research Foundation", "url": x["url"],
            "title": x["title"], "population": pp, "focus": "primary" if pop(x["title"]) else "related", "pi": [m.group(1).strip()] if m else [], "org": "",
            "country": "", "funder": "PolyBio Research Foundation", "start": dt, "end": None, "status": "active" if dt and dt >= "2024" else "unknown",
            "amount_usd": None, "kind": "Charity-funded", "n": None, "areas": [], "topics": topics(x["title"], text)})
# --- registered trials (site-src/data/trials.json), so the register covers every country with a trial
p = os.path.join(SRC, "data", "trials.json")
if os.path.exists(p):
    for t in json.load(open(p, encoding="utf8"))["trials"]:
        if t.get("result") == "withdrawn": continue
        end = d10(t.get("completion") or t.get("primary_completion"))
        st = "active" if t.get("result") == "pending" else "completed"
        rows.append({"id": "trial-" + t["id"], "source": "Trials", "source_label": "Trial registries", "url": t.get("url") or (f"https://clinicaltrials.gov/study/{t['id']}" if t["id"].startswith("NCT") else ""),
            "title": t.get("title") or t.get("intervention"), "population": t["population"], "focus": "primary", "pi": [], "org": t.get("sponsor") or "",
            "country": "; ".join(t["country"]) if isinstance(t.get("country"), list) else (t.get("country") or ""), "funder": "", "start": d10(t.get("start")), "end": end, "status": st, "amount_usd": None,
            "kind": "Trial: " + (t.get("design") or ""), "n": t.get("n_actual") or t.get("n_planned"), "areas": [], "topics": sorted(set(["trials"] + [k for k in topics(t.get("title") or "") if k != "other"]))})
CN = {"United States": "USA", "United States Of America": "USA", "US": "USA", "United Kingdom": "UK", "England": "UK", "Scotland": "UK", "Wales": "UK", "Great Britain": "UK", "Korea, Republic of": "South Korea", "Hong Kong": "China (Hong Kong)"}
for r in rows:
    r["country"] = "; ".join(dict.fromkeys(CN.get(c.strip(), c.strip()) for c in (r["country"] or "").split(";") if c.strip()))
rows.sort(key=lambda r: (r["status"] != "active", r["focus"] != "primary", -(int((r["start"] or "0")[:4] or 0))))
meta = {"generated": TODAY, "sources": {"MRR": "https://mrr.mecfs-research.org", "NIH": "https://reporter.nih.gov", "UKRI": "https://gtr.ukri.org", "CORDIS": "https://cordis.europa.eu", "Trials": "/trials", "PolyBio": "https://polybio.org/projects/"},
        "topics": {k: l for k, l, _ in TOPICS} | {"other": "Other"}}
json.dump({"meta": meta, "projects": rows}, open(os.path.join(SRC, "data", "projects.json"), "w", encoding="utf8"), ensure_ascii=False, indent=0)
from collections import Counter
print("projects", len(rows), Counter(r["source"] for r in rows), Counter(r["population"] for r in rows), Counter(r["status"] for r in rows))
