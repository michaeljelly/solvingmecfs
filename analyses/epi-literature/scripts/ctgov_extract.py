"""Extract posted results (primary+secondary outcomes, statistical analyses) from CT.gov v2 JSON."""
import json, pandas as pd, re
D = "/tmp/work_epi/ctgov"; OUTD = "/mnt/project-files/mecfs-research/deep2/epi_literature"
rows = []; trials = []
seen = set()
for cond in ["mecfs", "longcovid", "pots"]:
    for s in json.load(open(f"{D}/{cond}.json")):
        ps = s["protocolSection"]; idm = ps["identificationModule"]; nct = idm["nctId"]
        st = ps.get("statusModule", {}); des = ps.get("designModule", {})
        arms = ps.get("armsInterventionsModule", {})
        ivs = "; ".join(f"{i.get('type','')}:{i.get('name','')}" for i in arms.get("interventions", []))
        enr = des.get("enrollmentInfo", {}).get("count")
        alloc = des.get("designInfo", {}).get("allocation", ""); mask = des.get("designInfo", {}).get("maskingInfo", {}).get("masking", "")
        refs = ps.get("referencesModule", {}).get("references", [])
        pmids = ";".join(r.get("pmid", "") for r in refs if r.get("type") == "RESULT" and r.get("pmid"))
        trials.append(dict(cond_query=cond, nct=nct, title=idm.get("briefTitle"), status=st.get("overallStatus"), phase=";".join(des.get("phases", [])),
                           enrollment=enr, allocation=alloc, masking=mask, interventions=ivs, has_results=s.get("hasResults", False),
                           conditions="; ".join(ps.get("conditionsModule", {}).get("conditions", [])), result_pmids=pmids,
                           start=st.get("startDateStruct", {}).get("date"), completion=st.get("completionDateStruct", {}).get("date")))
        if not s.get("hasResults") or (nct, cond) in seen: continue
        seen.add((nct, cond))
        rs = s["resultsSection"]
        for om in rs.get("outcomeMeasuresModule", {}).get("outcomeMeasures", []):
            groups = {g["id"]: g.get("title", "") for g in om.get("groups", [])}
            # values
            vals = []
            for cl in om.get("classes", []):
                for cat in cl.get("categories", []):
                    for m in cat.get("measurements", []):
                        vals.append(f"{groups.get(m.get('groupId'),'')[:40]}={m.get('value')}({m.get('spread') or m.get('lowerLimit','')})")
            base = dict(cond_query=cond, nct=nct, title=idm.get("briefTitle"), interventions=ivs, enrollment=enr, alloc=alloc, masking=mask,
                        otype=om.get("type"), outcome=om.get("title"), timeframe=om.get("timeFrame"), param=om.get("paramType"),
                        unit=om.get("unitOfMeasure"), values=" | ".join(vals)[:600])
            ans = om.get("analyses", [])
            if not ans:
                rows.append({**base, "p": None, "method": None, "est": None, "est_type": None, "ci": None, "groups_cmp": None})
            for a in ans:
                rows.append({**base, "p": a.get("pValue"), "method": a.get("statisticalMethod"), "est": a.get("paramValue"),
                             "est_type": a.get("paramType"), "ci": f"{a.get('ciLowerLimit','')}..{a.get('ciUpperLimit','')}",
                             "groups_cmp": ",".join(groups.get(g, g)[:30] for g in a.get("groupIds", []))})
T = pd.DataFrame(trials).drop_duplicates("nct"); R = pd.DataFrame(rows)
T.to_csv(f"{OUTD}/ctgov_all_interventional_trials.csv", index=False)
R.to_csv(f"{OUTD}/ctgov_posted_results_outcomes.csv", index=False)
print(T.shape, R.shape, T.has_results.sum())
def pnum(p):
    if p is None or (isinstance(p, float)): return None
    m = re.search(r"([\d.]+(?:[eE]-?\d+)?)", str(p))
    return float(m.group(1)) if m else None
R["pnum"] = R.p.map(pnum)
prim = R[R.otype == "PRIMARY"]
summ = prim.groupby(["nct"]).agg(cond=("cond_query", "first"), title=("title", "first"), iv=("interventions", "first"), n=("enrollment", "first"),
      alloc=("alloc", "first"), n_primary=("outcome", "nunique"), min_p=("pnum", "min"), any_analysis=("p", lambda x: x.notna().any())).reset_index()
summ.to_csv(f"{OUTD}/ctgov_primary_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70); pd.set_option("display.max_rows", 300)
print(summ.sort_values(["cond", "min_p"]).to_string())
