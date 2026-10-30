import json, glob, sys, re
# usage: flatten_ctgov.py <rawdir> <prefix> <out.json>
rawdir, prefix, out = sys.argv[1:4]
rows = []
for f in sorted(glob.glob(f"{rawdir}/{prefix}_*.json")):
    for s in json.load(open(f))["studies"]:
        ps = s["protocolSection"]; idm = ps["identificationModule"]; st = ps["statusModule"]
        dm = ps.get("designModule", {}); di = dm.get("designInfo", {})
        aim = ps.get("armsInterventionsModule", {}); om = ps.get("outcomesModule", {})
        el = ps.get("eligibilityModule", {}); cl = ps.get("contactsLocationsModule", {})
        sec = [x.get("id") for x in idm.get("secondaryIdInfos", [])]
        countries = sorted({l.get("country") for l in cl.get("locations", []) if l.get("country")})
        refs = ps.get("referencesModule", {}).get("references", [])
        r = dict(
            nct=idm["nctId"], title=idm.get("officialTitle") or idm.get("briefTitle"), brief=idm.get("briefTitle"),
            acronym=idm.get("acronym"), secondary_ids=sec, org_id=idm.get("orgStudyIdInfo", {}).get("id"),
            conditions=ps.get("conditionsModule", {}).get("conditions", []),
            keywords=ps.get("conditionsModule", {}).get("keywords", []),
            interventions=[(i.get("type"), i.get("name"), (i.get("description") or "")[:200]) for i in aim.get("interventions", [])],
            arms=[(a.get("type"), a.get("label")) for a in aim.get("armGroups", [])],
            phases=dm.get("phases", []), allocation=di.get("allocation"), model=di.get("interventionModel"),
            masking=di.get("maskingInfo", {}).get("masking"), purpose=di.get("primaryPurpose"),
            enrollment=dm.get("enrollmentInfo", {}).get("count"), enroll_type=dm.get("enrollmentInfo", {}).get("type"),
            start=st.get("startDateStruct", {}).get("date"), completion=(st.get("primaryCompletionDateStruct") or st.get("completionDateStruct") or {}).get("date"),
            completion_full=st.get("completionDateStruct", {}).get("date"),
            status=st.get("overallStatus"), why_stopped=st.get("whyStopped"),
            sponsor=ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}).get("name"),
            countries=countries,
            primary_outcomes=[(o.get("measure"), o.get("timeFrame")) for o in om.get("primaryOutcomes", [])],
            criteria=el.get("eligibilityCriteria", "")[:4000],
            summary=ps.get("descriptionModule", {}).get("briefSummary", "")[:1500],
            refs=[(x.get("pmid"), x.get("type"), x.get("citation", "")[:300]) for x in refs],
            has_results=s.get("hasResults", False),
        )
        if s.get("hasResults"):
            res = []
            for o in s["resultsSection"].get("outcomeMeasuresModule", {}).get("outcomeMeasures", []):
                if o.get("type") != "PRIMARY": continue
                an = [(a.get("pValue"), a.get("statisticalMethod"), a.get("paramType"), a.get("paramValue"), a.get("ciLowerLimit"), a.get("ciUpperLimit"), (a.get("estimateComment") or a.get("statisticalComment") or "")[:150]) for a in o.get("analyses", [])]
                groups = {g["id"]: g["title"] for g in o.get("groups", [])}
                meas = []
                for c in o.get("classes", [])[:3]:
                    for cat in c.get("categories", [])[:2]:
                        for m in cat.get("measurements", []):
                            meas.append((c.get("title", ""), groups.get(m["groupId"], m["groupId"])[:40], m.get("value"), m.get("spread") or m.get("lowerLimit")))
                res.append(dict(title=o.get("title"), unit=o.get("unitOfMeasure"), param=o.get("paramType"), analyses=an, meas=meas[:12]))
            r["results_primary"] = res
            pf = s["resultsSection"].get("participantFlowModule", {})
            r["results_n"] = None
            try:
                bl = s["resultsSection"]["baselineCharacteristicsModule"]["denoms"][0]["counts"]
                r["results_n"] = sum(int(x["value"]) for x in bl if x["groupId"] != "BG%03d" % (len(bl) - 1)) if len(bl) > 1 else int(bl[0]["value"])
                tot = [x for x in bl]
                r["results_n"] = int(tot[-1]["value"]) if len(tot) > 1 else int(tot[0]["value"])
            except Exception:
                pass
        rows.append(r)
json.dump(rows, open(out, "w"))
print(len(rows))
