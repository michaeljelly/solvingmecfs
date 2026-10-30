"""For curated ME/CFS, long COVID and POTS treatment trials with posted results: compute standardized effects
(Hedges g between first two arms) for primary and key secondary outcomes from posted means/SD."""
import json, pandas as pd, numpy as np, re
D = "/tmp/work_epi/ctgov"; OUTD = "/mnt/project-files/mecfs-research/deep2/epi_literature"
keep = """NCT01071044 NCT03613129 NCT00498485 NCT03674541 NCT00598585 NCT00375973 NCT03891667 NCT02055898 NCT00977171
NCT04871815 NCT05199233 NCT05074888 NCT05453201 NCT05965726 NCT05965739 NCT05999435 NCT05126563 NCT05618574 NCT05618587 NCT05445427
NCT04652518 NCT05576662 NCT05592418 NCT04809974 NCT04876417 NCT05092516 NCT05121766 NCT05472090 NCT05633407 NCT05918978 NCT05874037
NCT05877508 NCT05946551 NCT06161688 NCT06214455 NCT06305780 NCT06524739 NCT05047952 NCT05200858 NCT05608629 NCT05595369 NCT05965752
NCT05205460 NCT06136871 NCT05096884 NCT05104749 NCT04663945 NCT05597722 NCT06476496 NCT06208696
NCT03484273 NCT03182725 NCT03919773 NCT00685919 NCT00770484 NCT06133075 NCT01547117 NCT01563107 NCT01978535 NCT03070730 NCT06268288
NCT00962949 NCT01064739 NCT01764711 NCT03602482 NCT01522235 NCT00555880 NCT01518946 NCT02154243 NCT03808740 NCT02500732""".split()
studies = {}
for c in ["mecfs", "longcovid", "pots"]:
    for s in json.load(open(f"{D}/{c}.json")):
        studies[s["protocolSection"]["identificationModule"]["nctId"]] = s
rows = []
for nct in keep:
    s = studies.get(nct)
    if not s or not s.get("hasResults"): continue
    ps = s["protocolSection"]; title = ps["identificationModule"]["briefTitle"]
    des = ps.get("designModule", {}); n = des.get("enrollmentInfo", {}).get("count")
    for om in s["resultsSection"].get("outcomeMeasuresModule", {}).get("outcomeMeasures", []):
        if om.get("type") not in ("PRIMARY", "SECONDARY"): continue
        groups = [g["id"] for g in om.get("groups", [])]; gt = {g["id"]: g.get("title", "") for g in om.get("groups", [])}
        dens = {}
        for d in om.get("denoms", []):
            for c_ in d.get("counts", []): dens[c_["groupId"]] = float(c_.get("value", "nan") or "nan")
        meas = {}
        cls = om.get("classes", [])
        if len(cls) >= 1:
            cats = cls[0].get("categories", [])
            if cats:
                for m in cats[0].get("measurements", []):
                    meas[m["groupId"]] = m
        ptxt = ";".join(str(a.get("pValue")) for a in om.get("analyses", []))
        g = None; diff = None
        disp = om.get("dispersionType", "")
        if len(groups) >= 2 and groups[0] in meas and groups[1] in meas and om.get("paramType") in ("MEAN", "LEAST_SQUARES_MEAN", "MEDIAN"):
            try:
                m1, m2 = float(meas[groups[0]]["value"]), float(meas[groups[1]]["value"])
                diff = m1 - m2
                if "Standard Deviation" in disp or disp == "STANDARD_DEVIATION":
                    s1, s2 = float(meas[groups[0]]["spread"]), float(meas[groups[1]]["spread"])
                    n1, n2 = dens.get(groups[0], np.nan), dens.get(groups[1], np.nan)
                    sp = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
                    g = diff / sp * (1 - 3 / (4 * (n1 + n2) - 9))
            except Exception: pass
        rows.append(dict(nct=nct, title=title[:80], n=n, otype=om["type"], outcome=om.get("title", "")[:150], timeframe=(om.get("timeFrame") or "")[:60],
                         param=om.get("paramType"), disp=disp, unit=om.get("unitOfMeasure"),
                         arm1=gt.get(groups[0], "")[:35] if groups else "", arm2=gt.get(groups[1], "")[:35] if len(groups) > 1 else "",
                         v1=meas.get(groups[0], {}).get("value") if groups else None, v2=meas.get(groups[1], {}).get("value") if len(groups) > 1 else None,
                         n1=dens.get(groups[0]) if groups else None, n2=dens.get(groups[1]) if len(groups) > 1 else None,
                         diff_arm1_minus_arm2=diff, hedges_g=g, posted_p=ptxt, n_groups=len(groups)))
R = pd.DataFrame(rows); R.to_csv(f"{OUTD}/ctgov_curated_effects.csv", index=False)
pd.set_option("display.width", 260); pd.set_option("display.max_colwidth", 60); pd.set_option("display.max_rows", 500)
P = R[R.otype == "PRIMARY"]
print(P[["nct", "title", "n", "outcome", "arm1", "arm2", "v1", "v2", "n1", "n2", "hedges_g", "posted_p", "n_groups"]].to_string())
