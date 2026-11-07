"""Refresh unfinished trials from their registries: status, sponsor, collaborators, investigators, sites.
Writes site-src/data/trials-live.json, keyed by trial id. Raw responses cached in raw/live/.
Usage: python3 -I enrich_live.py <mecfs-research dir>   (add --refresh to re-download)"""
import json, os, re, sys, time, html, urllib.request, urllib.parse, datetime
R = sys.argv[1].rstrip("/"); REF = "--refresh" in sys.argv
RAW = R + "/atlas-v2/trials/raw/live/"; os.makedirs(RAW, exist_ok=True)
T = json.load(open(R + "/site-src/data/trials.json", encoding="utf8"))["trials"]
CN = {"United States": "USA", "United Kingdom": "UK", "England": "UK", "Scotland": "UK", "Wales": "UK", "Korea, Republic of": "South Korea"}

def fetch(url, path):
    if REF or not os.path.exists(path):
        for a in range(4):
            try:
                b = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "solvingmecfs-trials/1.0"}), timeout=60).read()
                open(path, "wb").write(b); time.sleep(0.3); break
            except Exception as ex:
                if a == 3: print("fail", url, ex); return None
                time.sleep(2 ** a)
    return open(path, "rb").read()

def ctgov(nct):
    b = fetch(f"https://clinicaltrials.gov/api/v2/studies/{nct}", RAW + nct + ".json")
    if not b: return None
    p = json.loads(b)["protocolSection"]; sm = p["statusModule"]; sp = p.get("sponsorCollaboratorsModule", {})
    locs = p.get("contactsLocationsModule", {}).get("locations", [])
    offs = p.get("contactsLocationsModule", {}).get("overallOfficials", [])
    rp = sp.get("responsibleParty", {})
    pis = [o.get("name") for o in offs if o.get("name")] or ([rp["investigatorFullName"]] if rp.get("investigatorFullName") else [])
    return {"registry_status": sm.get("overallStatus"), "last_update": sm.get("lastUpdatePostDateStruct", {}).get("date"),
            "primary_completion": sm.get("primaryCompletionDateStruct", {}).get("date"), "completion": sm.get("completionDateStruct", {}).get("date"),
            "sponsor": sp.get("leadSponsor", {}).get("name"), "collaborators": [c["name"] for c in sp.get("collaborators", [])],
            "investigators": list(dict.fromkeys(re.sub(r"^(Prof\.?|Dr\.?|Professor) +|\s+(M\.?D\.?|Ph\.?D\.?|PT|MAS|MBBS)$", "", x.split(",")[0]).strip() for x in pis))[:3],
            "investigator_affiliations": list(dict.fromkeys(o.get("affiliation") for o in offs if o.get("affiliation")))[:3],
            "sites": [{"facility": l.get("facility"), "city": l.get("city"), "country": CN.get(l.get("country"), l.get("country")), "status": l.get("status")} for l in locs],
            "interventions": [i.get("name") for i in p.get("armsInterventionsModule", {}).get("interventions", []) if i.get("type") != "OTHER" or not re.match(r"(?i)placebo|sham|control|usual|standard", i.get("name", ""))],
            "enrolment": p.get("designModule", {}).get("enrollmentInfo", {}).get("count")}

def isrctn(i):
    b = fetch("https://www.isrctn.com/api/query/format/default?q=" + i, RAW + i + ".xml")
    if not b: return None
    x = b.decode("utf8", "ignore")
    g = lambda tag, s=x: [html.unescape(re.sub(r"<[^>]+>", " ", m)).strip() for m in re.findall(rf"<{tag}\b[^>]*>(.*?)</{tag}>", s, re.S)]
    pis = []
    for c in re.findall(r"<contact\b[^>]*>(.*?)</contact>", x, re.S):
        if "Principal investigator" in c or "Scientific" in c:
            n = " ".join(filter(None, g("forename", c)[:1] + g("surname", c)[:1]))
            if n and n not in ("Trial Lead",) and not re.search(r"(?i)trial|study|team", n): pis.append(n)
    cen = re.findall(r"<trialCentre\b[^>]*>(.*?)</trialCentre>", x, re.S)
    end = (g("overallEndDate") or [""])[0][:10]
    return {"registry_status": "Completed (end date passed)" if end and end < str(datetime.date.today()) else "Ongoing",
            "completion": end, "sponsor": (g("organisation", "".join(re.findall(r"<sponsor\b[^>]*>(.*?)</sponsor>", x, re.S))) or [None])[0], "collaborators": [],
            "funders": [f for f in g("name", "".join(re.findall(r"<funder\b[^>]*>(.*?)</funder>", x, re.S)))],
            "investigators": list(dict.fromkeys(pis))[:3],
            "sites": [{"facility": (g("name", c) or [""])[0], "city": (g("city", c) or [""])[0], "country": CN.get((g("country", c) or [""])[0], (g("country", c) or [""])[0])} for c in cen],
            "enrolment": None}

out = {}
for t in T:
    if t.get("result") != "pending": continue
    ids = t.get("all_ids") or [t["id"]]
    n = next((i for i in ids if i.startswith("NCT")), None); s = next((i for i in ids if i.startswith("ISRCTN")), None)
    rec = ctgov(n) if n else isrctn(s) if s else None
    if rec: rec["source_id"] = n or s; out[t["id"]] = rec
json.dump({"meta": {"generated": str(datetime.date.today()), "note": "Live registry fields for trials with result=pending; read by register/build_register.py and make_research.py"}, "trials": out},
          open(R + "/site-src/data/trials-live.json", "w", encoding="utf8"), ensure_ascii=False, indent=1)
print("live", len(out), "of", sum(1 for t in T if t.get("result") == "pending"))
