#!/usr/bin/env python3
"""Write pages/trials.html and <site>/trials-data.json from data/trials.json.
usage: python3 -I make_trials.py <site_src_dir> <site_dir>"""
import html, json, os, sys
from collections import Counter, defaultdict
SRC, SITE = sys.argv[1], sys.argv[2]
D = json.load(open(os.path.join(SRC, "data", "trials.json"), encoding="utf8"))
T = D["trials"]
e = lambda s: html.escape(str(s or ""), quote=True)
EFF = ("positive", "mixed", "null", "negative/harm")
def year(t):
    y = (t.get("start") or "")[:4]
    if not y:
        ys = [p.get("year") for p in t.get("publications") or [] if p.get("year")]
        y = min(ys) if ys else ""
    return int(y) if y.isdigit() else None
def n(t): return t.get("n_actual") or t.get("n_planned") or 0
def rig(t):
    d = t.get("design") or ""
    if d in ("single-arm", "open-label", "n-of-1"): return 0
    if t.get("blinded") and d in ("RCT-placebo", "crossover", "RCT-active"): return 3 if n(t) >= 100 else 2
    if d.startswith("RCT") or d == "crossover": return 1
    return None
COL = {"positive": "var(--rep)", "mixed": "var(--early)", "null": "var(--faint)", "negative/harm": "var(--ref)",
       "pending": "var(--rule)", "unreported-after-completion": "var(--cont)"}
LAB = {"positive": "Positive", "mixed": "Mixed", "null": "No benefit", "negative/harm": "Harm", "pending": "Not finished",
       "unreported-after-completion": "Finished, no result found"}

def bar(counts, keys, total=None):
    total = total or sum(counts.get(k, 0) for k in keys) or 1
    segs = "".join(f'<span style="width:{100*counts.get(k,0)/total:.2f}%;background:{COL[k]}" title="{LAB[k]}: {counts.get(k,0)}"></span>' for k in keys if counts.get(k, 0))
    return f'<span class="sb">{segs}</span>'
def legend(keys): return '<p class="lg">' + "".join(f'<span><i style="background:{COL[k]}"></i>{LAB[k]}</span>' for k in keys) + "</p>"

# Figure 1: result by rigour (verified efficacy results)
RL = ["Uncontrolled (single arm, open label)", "Randomised, not blinded", "Blinded, placebo or sham", "Blinded, 100 or more people"]
rc = defaultdict(Counter)
for t in T:
    if t.get("verified") and t.get("result") in EFF:
        r = rig(t)
        if r is None: continue
        rc[r][t["result"]] += 1
        if r == 3: rc[2][t["result"]] += 1
f1 = "".join(f'<div class="row"><span class="rl">{RL[i]}</span>{bar(rc[i], EFF)}<span class="rv">{100*rc[i]["positive"]//max(1,sum(rc[i].values()))}% positive · n={sum(rc[i].values())}</span></div>' for i in range(4))

# Figure 2: trials per year by population
yc = defaultdict(Counter)
for t in T:
    y = year(t)
    if y and y >= 1995: yc[y]["LC" if t["population"] == "long COVID" else "ME"] += 1
ys = list(range(min(yc), max(yc) + 1)); mx = max(sum(yc[y].values()) for y in ys)
cols = "".join(f'<div class="yc" title="{y}: {yc[y]["ME"]} ME/CFS, {yc[y]["LC"]} long COVID"><span style="height:{100*yc[y]["LC"]/mx:.1f}%;background:var(--cool)"></span><span style="height:{100*yc[y]["ME"]/mx:.1f}%;background:var(--warm)"></span><em>{str(y)[2:] if y % 5 == 0 else ""}</em></div>' for y in ys)

# Figure 3: by class
ALLK = list(EFF) + ["pending", "unreported-after-completion"]
cc = defaultdict(Counter)
for t in T: cc[t.get("intervention_class") or "other"][t.get("result")] += 1
order = sorted(cc, key=lambda k: -sum(cc[k].values()))
mxc = max(sum(cc[k].values()) for k in order)
f3 = "".join(f'<div class="row"><span class="rl">{e(k[0].upper()+k[1:])}</span><span class="sbw"><span style="width:{100*sum(cc[k][x] for x in ALLK)/mxc:.1f}%;display:block">{bar(cc[k], ALLK)}</span></span><span class="rv">{sum(cc[k].values())}</span></div>' for k in order)

# Tables
def link(t):
    i = t["id"]
    if i.startswith("NCT"): return f"https://clinicaltrials.gov/study/{i}"
    if t.get("doi"): return "https://doi.org/" + t["doi"]
    if t.get("pmid"): return f"https://pubmed.ncbi.nlm.nih.gov/{t['pmid']}/"
    return t.get("url") or ""
def nm(t): return t.get("acronym") or t["intervention"]
big = sorted([t for t in T if rig(t) == 3 and t.get("verified") and t.get("result") in EFF], key=lambda t: -n(t))
bigrows = "".join(f'<tr><td><a href="{e(link(t))}" rel="noopener">{e(nm(t))}</a></td><td>{e(t["population"])}</td><td class="num">{n(t)}</td><td><span class="dot" style="background:{COL[t["result"]]}"></span>{LAB[t["result"]]}</td><td>{e(t.get("result_summary"))}</td></tr>' for t in big[:25])
unrep = sorted([t for t in T if t.get("result") == "unreported-after-completion"], key=lambda t: -n(t))
unrows = "".join(f'<tr><td><a href="{e(link(t))}" rel="noopener">{e(nm(t))}</a></td><td>{e(t["population"])}</td><td class="num">{n(t)}</td><td>{e((t.get("completion") or "")[:7])}</td><td>{e(t.get("sponsor"))}</td></tr>' for t in unrep[:20])
pend = sorted([t for t in T if t.get("result") == "pending" and (t.get("design") or "").startswith("RCT")], key=lambda t: -n(t))
pendrows = "".join(f'<tr><td><a href="{e(link(t))}" rel="noopener">{e(nm(t))}</a></td><td>{e(t["population"])}</td><td class="num">{n(t)}</td><td>{e((t.get("completion") or t.get("primary_completion") or "")[:7])}</td><td>{e(t.get("sponsor"))}</td></tr>' for t in pend[:15])

# ---- Running trials: grade, context, and what looks promising (data/trials-watch.json, data/trials-live.json)
import re, datetime
TODAY = datetime.date.today().isoformat()
W = json.load(open(os.path.join(SRC, "data", "trials-watch.json"), encoding="utf8"))
LIVE = {}
_lp = os.path.join(SRC, "data", "trials-live.json")
if os.path.exists(_lp): LIVE = json.load(open(_lp, encoding="utf8")).get("trials", {})
TX = {x["id"]: x for x in json.load(open(os.path.join(SRC, "data", "treatments.json"), encoding="utf8"))["treatments"]}
STEP = D["meta"].get("model_steps_key", {})
FEAS = r"feasib|recruit|retention|acceptab|adherence|safety|adverse|tolerab|usability|satisfaction|fidelity|maximum tolerated"
OBJ = r"VO2|oxygen|CPET|exercise test|walk|6MWT|ESWT|steps|actigraph|upright|heart rate|blood volume|Cogstate|reaction time|cognit\w* (test|battery|index)|Symbol Digit|DSST|MMSE|Addenbrooke|NeuroTrax|Montreal|GSH|lactate"
GL = {"A": ("Could settle it", "var(--rep)"), "B": ("Informative", "var(--sev)"), "C": ("Early signal only", "var(--early)"), "D": ("Cannot show benefit", "var(--faint)")}
def grade(t):
    d, po, nn = t.get("design") or "", t.get("primary_outcome") or "", n(t)
    why = []
    if t.get("result") == "feasibility-only" or (po and re.match(r"^\W*(" + FEAS + ")", po, re.I)) or (po and re.search(FEAS, po, re.I) and not re.search(r"fatigue|function|quality|symptom|cognit", po, re.I)):
        return "D", "main outcome is feasibility or safety"
    if d in ("single-arm", "open-label", "n-of-1", "unclear", "") or d is None: return "D", "no control group"
    ctl = 3 if (t.get("blinded") and d in ("RCT-placebo", "crossover")) else (2 if t.get("blinded") else 1)
    why.append({3: "placebo or sham, blinded", 2: "randomised, assessor-blinded", 1: "randomised, not blinded"}[ctl])
    why.append(f"{nn} people" if nn else "size not stated")
    if po and re.search(OBJ, po, re.I): why.append("objective main outcome")
    elif po: why.append("self-rated main outcome")
    if ctl == 3 and nn >= 150: g = "A"
    elif (ctl == 3 and nn >= 60) or (ctl == 2 and nn >= 150): g = "B"
    else: g = "C"
    return g, ", ".join(why)
def due(t):
    lv = LIVE.get(t["id"]) or {}
    return (lv.get("completion") or t.get("primary_completion") or t.get("completion") or "")[:7]
def state(t):
    lv = LIVE.get(t["id"]) or {}
    rs = (lv.get("registry_status") or t.get("status") or "").replace("_", " ").lower()
    dd = due(t)
    if t.get("result") == "unreported-after-completion": return "Finished, no result found"
    if "complet" in rs or (dd and dd < TODAY[:7] and "recruit" not in rs): return "Finished, awaiting result"
    if "not yet" in rs: return "Not yet recruiting"
    if "active" in rs: return "Running, recruitment closed"
    return "Recruiting" if "recruit" in rs else "Running"
TXRE = [("pyridostigmine", r"pyridostigmine"), ("low-dose-naltrexone", r"naltrexone"), ("glp1", r"GLP-1|semaglutide|tirzepatide"), ("metformin", r"metformin"),
        ("colchicine", r"colchicine"), ("fluvoxamine", r"fluvoxamine"), ("low-dose-aripiprazole", r"aripiprazole"), ("ivig", r"immunoglobulin|IVIG"), ("ivabradine", r"ivabradine"),
        ("daratumumab", r"daratumumab"), ("rituximab", r"rituximab"), ("iv-saline", r"saline"), ("hbot-clinical", r"hyperbaric"), ("rapamycin", r"rapamycin|sirolimus"),
        ("tvns", r"vagus"), ("stellate-ganglion-block", r"stellate"), ("immunoadsorption", r"immunoadsorption|apheresis"), ("fmt", r"faecal|fecal|FMT"), ("probiotics", r"probiotic|microbiome therapy"),
        ("mito-supplements", r"NADH/NAD|nicotinamide|creatine|taurine"), ("covid-antivirals", r"nirmatrelvir|monoclonal|Paxlovid"), ("herpes-antivirals", r"valacyclovir|valganciclovir"),
        ("cell-therapy", r"cell therapy|MSC|cord blood"), ("steroids", r"corticosteroid|prednis|hydrocortisone"), ("stimulants", r"methylphenidate|modafinil|stimulant"), ("amantadine", r"amantadine"),
        ("cognitive-rehab-tdcs", r"cognitive rehabilitation|transcranial"), ("graded-exercise", r"graded exercise"), ("long-covid-rehab", r"exercise/rehabilitation"), ("cbt-supportive", r"CBT"),
        ("red-light", r"photobiomodulation"), ("lau7b-deupirfenidone", r"pirfenidone|fenretinide"), ("coq10-nadh", r"coenzyme Q10|CoQ10"), ("sleep-medicines", r"melatonin"),
        ("mast-cell", r"famotidine|loratadine|antihistamine"), ("pacing-heart-rate", r"pacing")]
def txof(t):
    w = (W["watch"].get(t["id"]) or {})
    if "tx" in w: return w["tx"]
    s = (t.get("intervention") or "") + " " + (t.get("title") or "")
    return next((k for k, p in TXRE if re.search(p, s, re.I) and k in TX), None)
VL = {"try": "Reasonable to try", "doctor": "Ask a doctor", "specialist": "Specialist only", "trial-only": "Only in a trial", "avoid": "Avoid"}
def chips(steps): return "".join(f'<a class="stp" href="/model#step-{s}" title="{e(STEP.get(str(s), ""))}">{s}</a>' for s in steps or [])
def gb(g): return f'<span class="gr" style="--g:{GL[g][1]}" title="{GL[g][0]}">{g}</span>'
RUN = [t for t in T if t.get("result") == "pending" or (t.get("result") == "unreported-after-completion" and (t.get("completion") or "") >= "2025-01")]
for t in RUN: t["_g"], t["_gw"] = grade(t)
gcount = Counter(t["_g"] for t in RUN)

# What looks promising
prom = ""
for tier in W["promising"]:
    items = ""
    for k in tier["ids"]:
        x = TX.get(k)
        if not x: continue
        items += (f'<li><a href="/treatments#{e(k)}"><b>{e(x["name"])}</b></a> <span class="vd v-{e(x["verdict"])}">{e(VL.get(x["verdict"], x["verdict"]))}</span> {chips(x.get("links"))}'
                  f'<span class="el">{e(x["evidence_line"])}</span>' + (f'<span class="pn">Next readout: {e(x["pending"])}</span>' if x.get("pending") else "") + "</li>")
    prom += f'<details class="tier"{" open" if tier is W["promising"][0] else ""}><summary><b>{e(tier["tier"])}</b> <span class="small">{len(tier["ids"])}</span><span class="tb">{e(tier["blurb"])}</span></summary><ul class="pl">{items}</ul></details>'

# Readouts to watch
wl = [t for t in T if t["id"] in W["watch"]]
for t in wl:
    if "_g" not in t: t["_g"], t["_gw"] = grade(t)
wl.sort(key=lambda t: (due(t) or "9999"))
cards = ""
for t in wl:
    w = W["watch"][t["id"]]; x = TX.get(w.get("tx") or "")
    st = state(t); late = st.startswith("Finished")
    prior = f'<p><b>Evidence so far.</b> {e(x["evidence_line"])} <a href="/treatments#{e(x["id"])}">Treatment page</a></p>' if x else ""
    cards += (f'<article class="wc{" late" if late else ""}"><header>{gb(t["_g"])}<div><h3><a href="{e(link(t))}" rel="noopener">{e(w.get("name") or t.get("acronym") or t["intervention"])}</a></h3>'
              f'<p class="wm">{e(t["population"])} · {e(t["_gw"])} · {e(t["id"])}</p></div>'
              f'<p class="wd"><span class="{"ov" if late else ""}">{e(st)}</span><br>{"due " + e(due(t)) if due(t) else ""}</p></header>'
              f'<p>{e(w["why"])} <span class="stps">Model steps {chips(t.get("model_steps") or (x or {}).get("links"))}</span></p>{prior}'
              f'<div class="io"><p><b class="ip">If it works</b>{e(w["pos"])}</p><p><b class="in">If it fails</b>{e(w["neg"])}</p></div></article>')

# Every running trial, graded
def rrow(t):
    x = TX.get(txof(t) or "")
    return (f'<tr><td>{gb(t["_g"])}</td><td><a href="{e(link(t))}" rel="noopener">{e(t.get("acronym") or t["intervention"])}</a>' + (f'<br><span class="note">{e(t["intervention"])}</span>' if t.get("acronym") else "") +
            f'</td><td>{e(t["population"])}</td><td class="num">{n(t) or ""}</td><td>{chips(t.get("model_steps"))}</td><td>{e(due(t))}<br><span class="note">{e(state(t))}</span></td>'
            f'<td>{f"""<a href="/treatments#{e(x["id"])}">{e(VL.get(x["verdict"]))}</a>""" if x else ""}</td></tr>')
RUN.sort(key=lambda t: ("ABCD".index(t["_g"]), -n(t)))
TH = '<thead><tr><th>Grade</th><th>Trial</th><th>Population</th><th>n</th><th>Steps</th><th>Due</th><th>Our verdict on the treatment</th></tr></thead>'
runA = "".join(rrow(t) for t in RUN if t["_g"] == "A"); runB = "".join(rrow(t) for t in RUN if t["_g"] == "B"); runCD = "".join(rrow(t) for t in RUN if t["_g"] in "CD")
gbar = "".join(f'<span style="width:{100*gcount[g]/max(1,len(RUN)):.1f}%;background:{GL[g][1]}" title="{g}: {gcount[g]}"></span>' for g in "ABCD")
glegend = "".join(f'<span><i style="background:{GL[g][1]}"></i>{g} · {GL[g][0]} ({gcount[g]})</span>' for g in "ABCD")
TW_CSS = """
.p-trials .tier{border-top:1px solid var(--rule);padding:.6rem 0}
.p-trials .tier summary{cursor:pointer;font-family:var(--f-sans);list-style:none}
.p-trials .tier summary::-webkit-details-marker{display:none}
.p-trials .tier summary::before{content:"+";display:inline-block;width:1rem;color:var(--muted)}
.p-trials .tier[open] summary::before{content:"−"}
.p-trials .tb{display:block;color:var(--muted);font-size:.85rem;margin:.2rem 0 0 1rem}
.p-trials .pl{list-style:none;padding:0;margin:.6rem 0 0 1rem}
.p-trials .pl li{margin:0 0 .8rem;font-size:.92rem;line-height:1.45}
.p-trials .pl .el{display:block;color:var(--ink);font-size:.86rem}
.p-trials .pl .pn{display:block;color:var(--muted);font-size:.8rem;font-family:var(--f-sans)}
.p-trials .vd{font-family:var(--f-sans);font-size:.72rem;padding:.05rem .4rem;border:1px solid var(--rule);border-radius:2px;color:var(--muted)}
.p-trials .vd.v-try{border-color:var(--rep);color:var(--rep)} .p-trials .vd.v-avoid{border-color:var(--ref);color:var(--ref)}
.p-trials .stp{display:inline-block;min-width:1.15rem;text-align:center;font-family:var(--f-sans);font-size:.7rem;line-height:1.15rem;border-radius:50%;background:var(--sunk);color:var(--ink);text-decoration:none;margin-left:.15rem}
.p-trials .gr{display:inline-block;width:1.5rem;height:1.5rem;line-height:1.5rem;text-align:center;font-family:var(--f-sans);font-weight:700;font-size:.8rem;color:var(--bg);background:var(--g);border-radius:2px}
.p-trials .wc{border-top:1.5px solid var(--ink);padding:.8rem 0 1rem}
.p-trials .wc.late{border-top-color:var(--warm)}
.p-trials .wc header{display:grid;grid-template-columns:1.9rem 1fr auto;gap:.7rem;align-items:start}
.p-trials .wc h3{margin:0;font-size:1.05rem} .p-trials .wc h3 a{text-decoration:none}
.p-trials .wm{margin:.1rem 0 0;color:var(--muted);font-family:var(--f-sans);font-size:.78rem}
.p-trials .wd{margin:0;text-align:right;font-family:var(--f-sans);font-size:.78rem;color:var(--muted)} .p-trials .wd .ov{color:var(--warm);font-weight:600}
.p-trials .wc>p{font-size:.93rem;margin:.6rem 0 0}
.p-trials .stps{font-family:var(--f-sans);font-size:.75rem;color:var(--muted);white-space:nowrap}
.p-trials .io{display:grid;grid-template-columns:1fr 1fr;gap:1rem;margin-top:.6rem;font-size:.88rem}
.p-trials .io p{margin:0} .p-trials .io b{display:block;font-family:var(--f-sans);font-size:.75rem;text-transform:uppercase;letter-spacing:.04em}
.p-trials .ip{color:var(--rep)} .p-trials .in{color:var(--ref)}
.p-trials .gscale{display:grid;grid-template-columns:repeat(4,1fr);gap:.8rem;font-family:var(--f-sans);font-size:.8rem;margin:.8rem 0}
.p-trials .gscale div{border-top:3px solid var(--g);padding-top:.35rem} .p-trials .gscale b{display:block}
@media (max-width:640px){.p-trials .io{grid-template-columns:1fr} .p-trials .gscale{grid-template-columns:1fr 1fr} .p-trials .wc header{grid-template-columns:1.9rem 1fr} .p-trials .wd{grid-column:2;text-align:left}}
"""

# Data for the browser
pub = [{"i": t["id"], "u": link(t), "t": t.get("acronym") or "", "x": t["intervention"], "c": t.get("intervention_class") or "other",
        "p": t["population"], "d": t.get("design") or "", "b": 1 if t.get("blinded") else 0, "n": n(t), "y": year(t) or "",
        "r": t.get("result") or "", "s": t.get("result_summary") or "", "v": 1 if t.get("verified") else 0, "st": t.get("model_steps") or []} for t in T]
json.dump(pub, open(os.path.join(SITE, "trials-data.json"), "w", encoding="utf8"), ensure_ascii=False, separators=(",", ":"))

cnt = Counter(t["population"] for t in T); res = Counter(t.get("result") for t in T)
nver = sum(1 for t in T if t.get("verified") and t.get("result") in EFF)
out = f'''<!--meta {{"title":"Treatment trials in ME/CFS and long COVID","description":"{len(T)} treatment trials in ME/CFS and long COVID, with design, size and result, and the trials that finished without reporting.","nav":"trials"}}-->
<style>
.p-trials .fig{{font-family:var(--f-sans);font-size:.8rem;margin:1rem 0}}
.p-trials .row{{display:grid;grid-template-columns:13rem 1fr 9rem;gap:.8rem;align-items:center;padding:.3rem 0}}
.p-trials .rl{{color:var(--ink)}} .p-trials .rv{{color:var(--muted);font-variant-numeric:tabular-nums}}
.p-trials .sb{{display:flex;height:.9rem;width:100%;background:var(--sunk)}} .p-trials .sb span{{display:block;height:100%}}
.p-trials .sbw{{display:block}}
.p-trials .lg{{display:flex;flex-wrap:wrap;gap:.3rem 1rem;color:var(--muted);margin:.2rem 0 0}} .p-trials .lg i{{display:inline-block;width:.7rem;height:.7rem;margin-right:.3rem;vertical-align:-.05rem}}
.p-trials .years{{display:flex;align-items:flex-end;gap:2px;height:9rem;border-bottom:1px solid var(--rule);padding-top:.5rem}}
.p-trials .yc{{flex:1;display:flex;flex-direction:column-reverse;height:100%;position:relative}} .p-trials .yc span{{display:block}}
.p-trials .yc em{{position:absolute;bottom:-1.2rem;left:50%;transform:translateX(-50%);font-style:normal;color:var(--muted);font-size:.7rem}}
.p-trials .dot{{display:inline-block;width:.55rem;height:.55rem;border-radius:50%;margin-right:.35rem}}
.p-trials table td{{font-size:.86rem;vertical-align:top}}
.p-trials .flt{{display:flex;flex-wrap:wrap;gap:.5rem;margin:1rem 0;font-family:var(--f-sans);font-size:.85rem;position:sticky;top:0;background:var(--bg);padding:.5rem 0;z-index:2}}
.p-trials .flt input,.p-trials .flt select{{font:inherit;padding:.35rem .5rem;border:1px solid var(--rule);background:var(--bg);color:var(--ink)}}
.p-trials .flt input{{flex:1 1 12rem}}
.p-trials #more{{font:inherit;font-family:var(--f-sans);padding:.4rem .8rem;border:1px solid var(--rule);background:none;color:var(--link);cursor:pointer}}
{TW_CSS}
@media (max-width:640px){{.p-trials .row{{grid-template-columns:1fr;gap:.2rem}} .p-trials .rv{{font-size:.75rem}}}}
</style>
<div class="wrap">
<h1>Treatment trials</h1>
<p class="deck">No treatment has yet shown a large, replicated, placebo-controlled benefit in ME/CFS. This page ranks the best bets, lists the results due soonest, and grades all {len(T)} trials in ME/CFS ({cnt["ME/CFS"]}), long COVID ({cnt["long COVID"]}) or both ({cnt["both"]}). Long COVID trials are included when they target fatigue, crashes after exertion, cognition or autonomic symptoms.</p>
<p class="byline">Built {e(D["meta"].get("generated","")[:10])} from ClinicalTrials.gov, ISRCTN, EU CTIS and Europe PMC · <a href="/trials-data.json">Download as JSON</a></p>

<nav class="small toc">On this page: <a href="#promising">What looks promising</a> · <a href="#watch">Readouts to watch</a> · <a href="#running">Every running trial, graded</a> · <a href="#design">Design and results</a> · <a href="#all">All trials</a></nav>

<h2 id="promising">What looks promising</h2>
<p>The best bets on current evidence, grouped by how strong that evidence is. Each links to its <a href="/treatments">treatment page</a>; the numbered circles are the <a href="/model">working model</a> steps it acts on.</p>
{prom}

<h2 id="watch">Readouts to watch</h2>
<p>The running trials whose results would most change what we think or do, soonest first. Each says which step of the model it tests and what we would conclude either way. Orange means the trial has finished and no result is public yet.</p>
{cards}

<h2 id="running">Every running trial, graded</h2>
<p>{len(RUN)} trials are running or finished in the last two years without a result. We grade each on what its design can show, not on whether the idea is good.</p>
<div class="gscale">{"".join(f'<div style="--g:{GL[g][1]}"><b>{g} · {GL[g][0]}</b>{d}</div>' for g, d in zip("ABCD", ["Blinded, placebo or sham, 150 or more people", "Blinded and 60 or more people, or assessor-blinded and 150 or more", "Randomised but small or unblinded", "No control group, or only feasibility or safety measured"]))}</div>
<figure class="fig"><span class="sb" style="height:1.1rem">{gbar}</span><p class="lg">{glegend}</p></figure>
<h3>Grade A: the trials that could settle a question</h3>
<div class="tbl"><table>{TH}<tbody>{runA}</tbody></table></div>
<details><summary class="small">Show the {gcount["B"]} grade B trials</summary><div class="tbl"><table>{TH}<tbody>{runB}</tbody></table></div></details>
<details><summary class="small">Show the {sum(1 for t in RUN if t["_g"] in "CD")} grade C and D trials</summary><div class="tbl"><table>{TH}<tbody>{runCD}</tbody></table></div></details>
<p class="note">Grades are computed from the registry: design, blinding, size and main outcome. A trial can be well designed and still test a weak idea; the cards above say which ideas matter. Steps are assigned from the intervention and may be wrong for unusual trials.</p>

<h2 id="design">Design and results</h2>
<p>The better a trial is designed, the less often it finds a benefit. Most uncontrolled trials report one; most large blinded trials do not. Of {res["unreported-after-completion"]} trials that finished more than a year ago, we found no published result.</p>

<figure class="fig l-page"><div>{f1}</div>{legend(EFF)}
<figcaption><b>Figure 1. Result by trial design.</b> {nver} trials whose efficacy result we read in the abstract or posted results. "Positive" means the authors reported a benefit on the primary outcome; it says nothing about bias. Traditional Chinese medicine and acupuncture trials, mostly unblinded, account for many of the uncontrolled positives.</figcaption></figure>

<figure class="fig l-page"><div class="years">{cols}</div><p class="lg" style="margin-top:1.6rem"><span><i style="background:var(--warm)"></i>ME/CFS</span><span><i style="background:var(--cool)"></i>Long COVID or both</span></p>
<figcaption><b>Figure 2. Trials by start year</b> (or first publication, where unregistered).</figcaption></figure>

<figure class="fig l-page"><div>{f3}</div>{legend(ALLK)}
<figcaption><b>Figure 3. Trials by type of treatment and result.</b> Bar length is the number of trials.</figcaption></figure>

<h2 id="large">Large blinded trials</h2>
<p>Placebo- or sham-controlled trials of 100 or more people, with a result.</p>
<div class="tbl"><table><thead><tr><th>Trial</th><th>Population</th><th>n</th><th>Result</th><th>Summary</th></tr></thead><tbody>{bigrows}</tbody></table></div>

<h2 id="unreported">Finished, no result found</h2>
<p>The largest of {res["unreported-after-completion"]} trials that completed more than 12 months ago with no posted or published result we could find. Some will have reported somewhere we missed; tell us.</p>
<div class="tbl"><table><thead><tr><th>Trial</th><th>Population</th><th>n</th><th>Completed</th><th>Sponsor</th></tr></thead><tbody>{unrows}</tbody></table></div>


<h2 id="all">All trials</h2>
<div class="flt"><input id="q" type="search" placeholder="Search drug, therapy or trial id" aria-label="Search trials">
<select id="fp" aria-label="Population"><option value="">All populations</option><option>ME/CFS</option><option>long COVID</option><option>both</option></select>
<select id="fr" aria-label="Result"><option value="">All results</option>{"".join(f'<option value="{k}">{LAB[k]}</option>' for k in ALLK)}</select>
<select id="fb" aria-label="Design"><option value="">All designs</option><option value="b">Blinded only</option></select></div>
<p class="note" id="cnt"></p>
<div class="tbl"><table><thead><tr><th>Treatment</th><th>Population</th><th>Design</th><th>n</th><th>Year</th><th>Result</th></tr></thead><tbody id="tb"></tbody></table></div>
<p><button id="more" type="button">Show more</button></p>
<p class="note">Results marked with * were classified from the registry or title without reading an abstract. Classes and model steps are assigned by keyword. Corrections to <a href="mailto:hello@solvingmecfs.org">hello@solvingmecfs.org</a>.</p>
</div>
<script>
(function(){{var D=[],lim=60,L={json.dumps(LAB)},C={json.dumps(COL)};
function esc(s){{return String(s).replace(/[&<>"]/g,function(c){{return{{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}}[c]}})}}
function draw(){{var q=document.getElementById('q').value.toLowerCase(),p=document.getElementById('fp').value,r=document.getElementById('fr').value,b=document.getElementById('fb').value;
var f=D.filter(function(t){{return(!p||t.p===p)&&(!r||t.r===r)&&(!b||t.b)&&(!q||(t.x+' '+t.t+' '+t.i+' '+t.s+' '+t.c).toLowerCase().indexOf(q)>=0)}});
document.getElementById('cnt').textContent=f.length+' trials';
document.getElementById('tb').innerHTML=f.slice(0,lim).map(function(t){{return '<tr><td><a href="'+esc(t.u)+'" rel="noopener">'+esc(t.t||t.x)+'</a>'+(t.t?'<br><span class="note">'+esc(t.x)+'</span>':'')+'</td><td>'+esc(t.p)+'</td><td>'+esc(t.d)+(t.b?', blinded':'')+'</td><td class="num">'+(t.n||'')+'</td><td class="num">'+t.y+'</td><td>'+(L[t.r]?'<span class="dot" style="background:'+C[t.r]+'"></span>'+L[t.r]+(t.v?'':'*'):esc(t.r))+(t.s?'<br><span class="note">'+esc(t.s)+'</span>':'')+'</td></tr>'}}).join('');
document.getElementById('more').style.display=f.length>lim?'':'none'}}
['q','fp','fr','fb'].forEach(function(id){{document.getElementById(id).addEventListener('input',function(){{lim=60;draw()}})}});
document.getElementById('more').addEventListener('click',function(){{lim+=100;draw()}});
fetch('/trials-data.json').then(function(r){{return r.json()}}).then(function(d){{D=d.sort(function(a,b){{return(b.n||0)-(a.n||0)}});draw()}});}})();
</script>
'''
open(os.path.join(SRC, "pages", "trials.html"), "w", encoding="utf8").write(out)
print("built trials", len(T), dict(rc[0]), dict(rc[2]), dict(rc[3]))
