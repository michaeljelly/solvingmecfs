import re, json, glob, sys, html
base = sys.argv[1]
rd = f"{base}/raw/registries"
def g(t, tag):
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", t, re.S)
    return html.unescape(re.sub(r"\s+", " ", m.group(1))).strip() if m else ""
def ga(t, tag):
    return [html.unescape(re.sub(r"\s+", " ", x)).strip() for x in re.findall(rf"<{tag}[^>]*>(.*?)</{tag}>", t, re.S)]
out = {}
for f in sorted(glob.glob(f"{rd}/isrctn_*.xml")):
    t = open(f, encoding="utf-8").read()
    for tr in t.split("<fullTrial>")[1:]:
        iid = "ISRCTN" + g(tr, "isrctn")
        if iid in out: continue
        cond = tr.split("<conditions>")[1].split("</conditions>")[0] if "<conditions>" in tr else ""
        ivs = tr.split("<interventions>")[1].split("</interventions>")[0] if "<interventions>" in tr else ""
        po = tr.split("<primaryOutcomes>")[1].split("</primaryOutcomes>")[0] if "<primaryOutcomes>" in tr else ""
        ext = tr.split("<externalRefs>")[1].split("</externalRefs>")[0] if "<externalRefs>" in tr else ""
        outputs = tr.split("<outputs>")[1].split("</outputs>")[0] if "<outputs>" in tr else ""
        out[iid] = {
            "id": iid, "title": g(tr, "title"), "scientific_title": g(tr, "scientificTitle"), "acronym": g(tr, "acronym"),
            "design": g(tr, "primaryStudyDesign") + " / " + g(tr, "secondaryStudyDesign"),
            "conditions": g(cond, "description") or g(cond, "condition"),
            "intervention": g(ivs, "description")[:400], "intervention_type": g(ivs, "interventionType"),
            "primary_outcome": g(po, "variable")[:300] or g(po, "method")[:300],
            "start": g(tr, "recruitmentStart")[:10], "end": g(tr, "overallEndDate")[:10],
            "target_n": g(tr, "targetEnrolment"), "final_n": g(tr, "totalFinalEnrolment"),
            "countries": ga(tr, "country")[:3], "sponsor": g(tr.split("<sponsor>")[1] if "<sponsor>" in tr else "", "organisation"),
            "ext_refs": re.findall(r"NCT\d{8}|\d{4}-\d{6}-\d{2}", ext + tr[:5000]),
            "outputs": re.findall(r"(?:pubmed|doi)[^<\"]{0,80}", outputs)[:5],
            "overall_status_hint": g(tr, "overallStatusOverride") or "",
        }
json.dump(list(out.values()), open(f"{base}/work/isrctn_flat.json", "w"), indent=0)
ct = {}
for f in sorted(glob.glob(f"{rd}/ctis_*.json")):
    for r in json.load(open(f)).get("data", []):
        ct.setdefault(r["ctNumber"], r)
json.dump(list(ct.values()), open(f"{base}/work/ctis_flat.json", "w"), indent=0)
print(len(out), len(ct))
