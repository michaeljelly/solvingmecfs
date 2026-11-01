#!/usr/bin/env python3
"""Write pages/research.html, <site>/research-data.json and <site>/orgs-data.json
from data/projects.json (register/build_register.py) and ../kb/registry.json (data thread).
usage: python3 -I make_research.py <site_src_dir> <site_dir>"""
import html, json, os, re, sys
from collections import Counter, defaultdict
SRC, SITE = sys.argv[1], sys.argv[2]
e = lambda s: html.escape(str(s or ""), quote=True)
D = json.load(open(os.path.join(SRC, "data", "projects.json"), encoding="utf8"))
P, TOP = D["projects"], D["meta"]["topics"]
RP = os.path.join(SRC, "..", "kb", "registry.json")
ORG = [o for o in json.load(open(RP, encoding="utf8"))["entities"] if o["type"] != "funding_call"] if os.path.exists(RP) else []
TYPE = {"research_group": "Research group", "research_network": "Research network", "biobank": "Biobank", "cohort": "Cohort study",
        "clinic": "Clinic", "data_portal": "Data portal", "registry": "Patient registry", "patient_org": "Patient organisation",
        "funder": "Funder", "government_plan": "Government plan", "initiative": "Initiative", "conference": "Conference",
        "journal": "Journal", "industry_sponsor": "Company", "trial_network": "Trial network"}
PC = {"ME/CFS": "var(--warm)", "both": "var(--early)", "long COVID": "var(--cool)"}
act = [p for p in P if p["status"] == "active"]
# Figure 1: active projects by topic and population
tc = defaultdict(Counter)
for p in act:
    for t in p["topics"]: tc[t][p["population"]] += 1
order = sorted(tc, key=lambda k: -sum(tc[k].values()))
mx = max(sum(tc[k].values()) for k in order)
def sbar(c, scale):
    return "".join(f'<span style="width:{100 * c[k] / scale:.2f}%;background:{PC[k]}" title="{c[k]} {k}"></span>' for k in PC if c[k])
f1 = "".join(f'<a class="row" href="#all" data-topic="{k}"><span class="rl">{e(TOP.get(k, k))}</span><span class="sb">{sbar(tc[k], mx)}</span><span class="rv">{sum(tc[k].values())}</span></a>' for k in order)
# Figure 2: active projects by country
cc = defaultdict(Counter)
for p in act:
    for c in (p["country"] or "Unknown").split("; ")[:1]: cc[c or "Unknown"][p["population"]] += 1
corder = [k for k in sorted(cc, key=lambda k: -sum(cc[k].values())) if k != "Unknown"][:12]
mxc = max(sum(cc[k].values()) for k in corder)
f2 = "".join(f'<div class="row"><span class="rl">{e(k)}</span><span class="sb">{sbar(cc[k], mxc)}</span><span class="rv">{sum(cc[k].values())}</span></div>' for k in corder)
legend = '<p class="lg">' + "".join(f'<span><i style="background:{v}"></i>{k}</span>' for k, v in PC.items()) + "</p>"
# Gaps: how many active projects touch what the model says matters most
def cnt(rx): return sum(1 for p in act if re.search(rx, p["title"], re.I))
g_sev = cnt(r"severe|bedbound|housebound|homebound|very ill")
g_pem = cnt(r"post-?exertional|\bPEM\b|exertion|exercise|crash")
g_early = cnt(r"early|onset|incident|prospective|after infection|acute")
srcs = Counter(p["source"] for p in P)
pub = [{"t": p["title"], "u": p["url"], "p": p["population"], "f": 1 if p["focus"] == "primary" else 0, "o": p["org"], "c": p["country"],
        "pi": ", ".join(p["pi"][:3]), "fu": p["funder"], "s": (p["start"] or "")[:4], "en": (p["end"] or "")[:4], "st": p["status"],
        "tp": p["topics"], "src": p["source"], "a": p.get("amount_usd") or 0} for p in P]
json.dump(pub, open(os.path.join(SITE, "research-data.json"), "w", encoding="utf8"), ensure_ascii=False, separators=(",", ":"))
opub = [{"n": o["name"], "u": o["url"] if o.get("link_status") in ("ok", "ok-no-mention", "blocked") or o.get("link_status", "").startswith("unreachable") else "",
         "ty": o["type"], "c": o.get("country") or "", "d": o.get("description") or "", "st": o.get("status") or ""} for o in ORG if o.get("link_status") not in ("hijacked", "dead")]
json.dump(opub, open(os.path.join(SITE, "orgs-data.json"), "w", encoding="utf8"), ensure_ascii=False, separators=(",", ":"))
otc = Counter(o["ty"] for o in opub)
out = f'''<!--meta {{"title":"Who is working on what","description":"{len(P)} research projects on ME/CFS and long COVID, and {len(opub)} research groups, networks, biobanks, registries and organisations, searchable by topic, country and status.","nav":"research"}}-->
<style>
.p-research .fig{{font-family:var(--f-sans);font-size:.8rem;margin:1rem 0}}
.p-research .row{{display:grid;grid-template-columns:13rem 1fr 3rem;gap:.8rem;align-items:center;padding:.28rem 0;color:var(--ink);text-decoration:none}}
.p-research a.row:hover .rl{{text-decoration:underline}}
.p-research .rv{{color:var(--muted);font-variant-numeric:tabular-nums;text-align:right}}
.p-research .sb{{display:flex;height:.8rem}} .p-research .sb span{{display:block;height:100%}} .p-research .sb span+span{{border-left:1.5px solid var(--bg)}}
.p-research .lg{{display:flex;flex-wrap:wrap;gap:.3rem 1rem;color:var(--muted);margin:.4rem 0 0}} .p-research .lg i{{display:inline-block;width:.7rem;height:.7rem;margin-right:.3rem;vertical-align:-.05rem}}
.p-research .stats{{display:grid;grid-template-columns:repeat(3,1fr);gap:1rem;margin:1.2rem 0;font-family:var(--f-sans)}}
.p-research .stats div{{border-top:1.5px solid var(--ink);padding-top:.4rem}} .p-research .stats b{{display:block;font-size:1.6rem;font-weight:600;font-variant-numeric:lining-nums}} .p-research .stats span{{font-size:.8rem;color:var(--muted)}}
.p-research .flt{{display:flex;flex-wrap:wrap;gap:.5rem;margin:1rem 0 .4rem;font-family:var(--f-sans);font-size:.85rem;position:sticky;top:0;background:var(--bg);padding:.5rem 0;z-index:2}}
.p-research .flt input,.p-research .flt select{{font:inherit;padding:.35rem .5rem;border:1px solid var(--rule);background:var(--bg);color:var(--ink);max-width:100%}}
.p-research .flt input{{flex:1 1 12rem}}
.p-research .pl{{list-style:none;padding:0;margin:0}}
.p-research .pl li{{border-top:1px solid var(--rule);padding:.6rem 0;display:grid;grid-template-columns:1fr auto;gap:.15rem 1rem}}
.p-research .pl .tt{{font-family:var(--f-sans);font-size:.9rem;font-weight:600;line-height:1.35}} .p-research .pl .tt a{{color:var(--ink);text-decoration:none}} .p-research .pl .tt a:hover{{color:var(--link)}}
.p-research .pl .mt{{grid-column:1/-1;font-family:var(--f-sans);font-size:.78rem;color:var(--muted);line-height:1.4}}
.p-research .pl .pp{{font-family:var(--f-sans);font-size:.74rem;white-space:nowrap;font-weight:600}}
.p-research .pl .tg{{display:inline-block;margin-right:.5rem;color:var(--ink-2)}}
.p-research .pl .ds{{grid-column:1/-1;font-size:.86rem;line-height:1.4;margin:.1rem 0 0;color:var(--ink-2)}}
.p-research .more{{font:inherit;font-family:var(--f-sans);padding:.4rem .8rem;border:1px solid var(--rule);background:none;color:var(--link);cursor:pointer;margin-top:.8rem}}
@media (max-width:640px){{.p-research .row{{grid-template-columns:1fr 2.5rem;gap:.1rem .6rem}} .p-research .row .rl{{grid-column:1/-1}} .p-research .stats{{grid-template-columns:1fr 1fr 1fr;gap:.6rem}} .p-research .stats b{{font-size:1.2rem}}}}
</style>
<div class="wrap">
<h1>Who is working on what</h1>
<p class="deck">{len(P):,} research projects on ME/CFS and long COVID, and {len(opub)} research groups, networks, biobanks, registries and organisations.</p>
<p class="byline">Built {e(D["meta"]["generated"])} from the <a href="https://mrr.mecfs-research.org">ME/CFS Research Register</a> ({srcs["MRR"]} European projects), <a href="https://reporter.nih.gov">NIH RePORTER</a> ({srcs["NIH"]} grants since 2016){f', <a href="https://gtr.ukri.org">UKRI Gateway to Research</a> ({srcs["UKRI"]})' if srcs["UKRI"] else ""}, <a href="https://cordis.europa.eu">EU CORDIS</a> ({srcs["CORDIS"]}), <a href="https://polybio.org/projects/">PolyBio</a> ({srcs["PolyBio"]}), the <a href="/trials">trial registries</a> ({srcs["Trials"]} trials) and our own survey of organisations · <a href="/research-data.json">Projects as JSON</a> · <a href="/orgs-data.json">Organisations as JSON</a></p>

<div class="stats"><div><b>{len(act)}</b><span>active projects</span></div><div><b>{g_sev}</b><span>name severe or housebound patients in the title</span></div><div><b>{g_pem}</b><span>name exertion, exercise or PEM in the title</span></div></div>
<p>About a quarter of patients are housebound or bedbound, yet {g_sev} of {len(act)} active projects say in their title that they study severe patients. Crashes after exertion define the illness, and {g_pem} name exertion or PEM. The <a href="/questions">open questions</a> page sets out which experiments would most change what we know.</p>

<h2 id="topics">Active projects by topic</h2>
<figure class="fig"><div>{f1}</div>{legend}
<figcaption><b>Figure 1.</b> Active projects, by topic and population. A project can count under several topics; topics are assigned from the title and abstract by keyword. Tap a topic to list its projects.</figcaption></figure>
<h2 id="countries">Where the work is</h2>
<figure class="fig"><div>{f2}</div>{legend}
<figcaption><b>Figure 2.</b> Active projects by lead country, top 12. Trials are counted for every country. Other research is complete only where a public grants database exists (NIH, UKRI, the EU, and the six countries of the European register), so the United States and Germany are over-represented, and national funders in Australia, Canada and Asia are not yet included.</figcaption></figure>

<h2 id="all">All projects</h2>
<div class="flt"><input id="q" type="search" placeholder="Search title, lab, person or institution" aria-label="Search projects">
<select id="fp" aria-label="Population"><option value="">ME/CFS and long COVID</option><option>ME/CFS</option><option>long COVID</option><option>both</option></select>
<select id="ft" aria-label="Topic"><option value="">All topics</option>{"".join(f'<option value="{k}">{e(v)}</option>' for k, v in TOP.items())}</select>
<select id="fs" aria-label="Status"><option value="active">Active</option><option value="">Active and finished</option></select>
<select id="ff" aria-label="Focus"><option value="1">Main focus</option><option value="">Including related work</option></select></div>
<p class="note" id="cnt"></p>
<ul class="pl" id="pl"></ul>
<button class="more" id="more" type="button">Show more</button>

<h2 id="orgs">Groups, networks and resources</h2>
<p>Research groups, networks, biobanks, cohorts, registries, data portals, patient organisations, government plans and conferences, each with a checked link. Funders and their open calls are on the <a href="/funding">funding</a> page.</p>
<div class="flt"><input id="oq" type="search" placeholder="Search name, country or description" aria-label="Search organisations">
<select id="ot" aria-label="Type"><option value="">All types</option>{"".join(f'<option value="{k}">{e(TYPE.get(k, k))} ({n})</option>' for k, n in otc.most_common())}</select></div>
<p class="note" id="ocnt"></p>
<ul class="pl" id="ol"></ul>
<button class="more" id="omore" type="button">Show more</button>
<p class="note">Missing a project, or something out of date? Write to <a href="mailto:hello@solvingmecfs.org">hello@solvingmecfs.org</a>.</p>
</div>
<script>
(function(){{var P=[],O=[],lim=40,olim=40,TOP={json.dumps(TOP)},TY={json.dumps(TYPE)},PC={json.dumps(PC)};
function $(i){{return document.getElementById(i)}}
function esc(s){{return String(s==null?'':s).replace(/[&<>"]/g,function(c){{return{{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}}[c]}})}}
function money(a){{return a>=1e6?'$'+(a/1e6).toFixed(1)+'M':a>=1e3?'$'+Math.round(a/1e3)+'k':''}}
function draw(){{var q=$('q').value.toLowerCase(),p=$('fp').value,t=$('ft').value,s=$('fs').value,f=$('ff').value;
var r=P.filter(function(x){{return(!p||x.p===p)&&(!t||x.tp.indexOf(t)>=0)&&(!s||x.st===s)&&(!f||x.f)&&(!q||(x.t+' '+x.o+' '+x.pi+' '+x.c+' '+x.fu).toLowerCase().indexOf(q)>=0)}});
$('cnt').textContent=r.length+' projects';
$('pl').innerHTML=r.slice(0,lim).map(function(x){{return '<li><span class="tt"><a href="'+esc(x.u)+'" rel="noopener">'+esc(x.t)+'</a></span><span class="pp" style="color:'+PC[x.p]+'">'+esc(x.p)+'</span><span class="mt">'+esc([x.pi,x.o,x.c].filter(Boolean).join(' · '))+'<br>'+esc([x.fu,(x.s||'?')+'–'+(x.en||''),money(x.a),x.st==='active'?'':'finished'].filter(Boolean).join(' · '))+' · '+x.tp.map(function(k){{return '<span class="tg">'+esc(TOP[k]||k)+'</span>'}}).join('')+'</span></li>'}}).join('');
$('more').style.display=r.length>lim?'':'none'}}
function odraw(){{var q=$('oq').value.toLowerCase(),t=$('ot').value;
var r=O.filter(function(x){{return(!t||x.ty===t)&&(!q||(x.n+' '+x.c+' '+x.d).toLowerCase().indexOf(q)>=0)}});
$('ocnt').textContent=r.length+' entries';
$('ol').innerHTML=r.slice(0,olim).map(function(x){{return '<li><span class="tt">'+(x.u?'<a href="'+esc(x.u)+'" rel="noopener">'+esc(x.n)+'</a>':esc(x.n))+'</span><span class="pp">'+esc(TY[x.ty]||x.ty)+'</span><span class="mt">'+esc([x.c,x.st&&x.st!=='active'?x.st:''].filter(Boolean).join(' · '))+'</span>'+(x.d?'<p class="ds">'+esc(x.d)+'</p>':'')+'</li>'}}).join('');
$('omore').style.display=r.length>olim?'':'none'}}
['q','fp','ft','fs','ff'].forEach(function(i){{$(i).addEventListener('input',function(){{lim=40;draw()}})}});
['oq','ot'].forEach(function(i){{$(i).addEventListener('input',function(){{olim=40;odraw()}})}});
$('more').addEventListener('click',function(){{lim+=60;draw()}});$('omore').addEventListener('click',function(){{olim+=60;odraw()}});
document.querySelectorAll('a.row[data-topic]').forEach(function(a){{a.addEventListener('click',function(){{$('ft').value=a.dataset.topic;$('fs').value='active';lim=40;draw()}})}});
var hm=location.hash.match(/^#t=([a-z-]+)/);if(hm){{$('ft').value=hm[1];}}
fetch('/research-data.json').then(function(r){{return r.json()}}).then(function(d){{P=d;draw();if(hm)$('all').scrollIntoView()}});
fetch('/orgs-data.json').then(function(r){{return r.json()}}).then(function(d){{O=d.sort(function(a,b){{return a.n.localeCompare(b.n)}});odraw()}});}})();
</script>
'''
json.dump({k: sum(tc[k].values()) for k in tc}, open(os.path.join(SRC, "data", "project-counts.json"), "w"))
open(os.path.join(SRC, "pages", "research.html"), "w", encoding="utf8").write(out)
print("built research", len(P), "projects,", len(opub), "orgs; active", len(act), "severe", g_sev, "pem", g_pem)
