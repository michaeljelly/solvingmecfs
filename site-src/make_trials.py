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
@media (max-width:640px){{.p-trials .row{{grid-template-columns:1fr;gap:.2rem}} .p-trials .rv{{font-size:.75rem}}}}
</style>
<div class="wrap">
<h1>Treatment trials</h1>
<p class="deck">{len(T)} treatment trials in ME/CFS ({cnt["ME/CFS"]}), long COVID ({cnt["long COVID"]}) or both ({cnt["both"]}), with design, size and result. Long COVID trials are included when they target fatigue, crashes after exertion, cognition or autonomic symptoms.</p>
<p class="byline">Built {e(D["meta"].get("generated","")[:10])} from ClinicalTrials.gov, ISRCTN, EU CTIS and Europe PMC · <a href="/trials-data.json">Download as JSON</a></p>

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

<h2 id="pending">Largest randomised trials still running</h2>
<div class="tbl"><table><thead><tr><th>Trial</th><th>Population</th><th>n</th><th>Due</th><th>Sponsor</th></tr></thead><tbody>{pendrows}</tbody></table></div>

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
