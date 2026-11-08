#!/usr/bin/env python3
"""Write pages/library.html and <site>/library-data.json from lit/library.json and lit/screen/*.json.
usage: python3 -I make_library.py <site_src_dir> <site_dir> <lit_dir>"""
import datetime, glob, html, json, os, sys
from collections import Counter
SRC, SITE, LIT = sys.argv[1:4]
e = lambda s: html.escape(str(s or ""), quote=True)
L = json.load(open(os.path.join(LIT, "library.json"), encoding="utf8"))
S = {}
for f in sorted(glob.glob(os.path.join(LIT, "screen", "*.json"))):
    try:
        for x in json.load(open(f, encoding="utf8")):
            if isinstance(x, dict) and x.get("id"): S[str(x["id"])] = x
    except Exception as ex:
        print("skip", f, ex)
F = {}
for p in glob.glob(os.path.join(os.path.dirname(LIT), "atlas", "build", "research", "*.json")):
    if os.path.basename(p) in ("overview.json", "hypotheses.json"): continue
    for x in json.load(open(p, encoding="utf8")).get("findings", []): F[x["id"]] = x
DL = {"vascular": "Blood flow and vessels", "muscle-exertion": "Muscle and exertion", "neuro": "Brain and nerves", "immune": "Immune system", "infection": "Infection",
      "omics": "Omics and biomarkers", "endocrine": "Hormones", "gut": "Gut", "genetics": "Genetics", "trials": "Treatment", "epi": "Epidemiology and care"}
KL = {"primary": "Primary study", "trial": "Trial", "meta": "Systematic review or meta-analysis", "review": "Review", "commentary": "Commentary or letter",
      "case": "Case report", "protocol": "Protocol", "no-abstract": "No abstract", "other": "Erratum or other"}
RL = {"replicates": "Replicates", "contradicts": "Contradicts", "extends": "Extends", "new": "New finding", "none": "Background"}
def url(r):
    if r.get("doi"): return "https://doi.org/" + r["doi"]
    if r.get("pmid"): return f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/"
    return f"https://europepmc.org/article/{r.get('src')}/{r.get('id')}"
def au(a):
    a = (a or "").split(",")
    return (a[0].strip() + (" et al." if len(a) > 1 else "")) if a and a[0].strip() else ""
pub, rows = [], []
for r in L:
    s = S.get(str(r["id"])) or {}
    rel = s.get("relation") if s.get("relevant") in ("core", "related") else None
    d = {"i": r["id"], "u": url(r), "t": r["title"], "a": au(r.get("authors")), "j": r.get("journal") or ("Preprint" if r.get("src") == "PPR" else ""), "y": r.get("year") or "",
         "dt": r.get("date") or "", "p": r["pop"], "k": r["kind"], "pp": 1 if r.get("src") == "PPR" else 0, "dm": r.get("domains") or [], "c": 1 if r.get("cited") else 0,
         "sc": 1 if s else 0, "rv": s.get("relevant") or "", "r": rel or "", "f": s.get("finding") or "", "at": [x for x in (s.get("atlas") or []) if x in F][:3],
         "im": s.get("importance") or 0, "n": s.get("n") or "", "ds": s.get("design") or "", "nt": s.get("note") or ""}
    pub.append(d)
json.dump(pub, open(os.path.join(SITE, "library-data.json"), "w", encoding="utf8"), ensure_ascii=False, separators=(",", ":"))

tot = len(pub); scr = sum(d["sc"] for d in pub); cit = sum(d["c"] for d in pub); pre = sum(d["pp"] for d in pub)
core = sum(1 for d in pub if d["k"] in ("primary", "trial", "meta"))
since = (datetime.date.today() - datetime.timedelta(days=90)).isoformat()
def item(d, show_rel=True):
    at = "".join(f' <a class="fa" href="/atlas#f-{e(x)}" title="{e(F[x]["finding"])}">{e(F[x].get("subdomain") or x)}</a>' for x in d["at"])
    rel = f'<span class="rl r-{e(d["r"])}">{RL.get(d["r"], "")}</span> ' if (show_rel and d["r"] and d["r"] != "none") else ""
    return (f'<li>{rel}<a href="{e(d["u"])}" rel="noopener">{e(d["t"])}</a>'
            f'<span class="m">{e(d["a"])} · {e(d["j"])} {e(d["y"])}{" · " + e(d["n"]) if d["n"] else ""}{" · preprint" if d["pp"] and d["j"] != "Preprint" else ""}</span>'
            + (f'<span class="fd">{e(d["f"])}{at}</span>' if d["f"] else "") + "</li>")
key = [d for d in pub if d["im"] == 3 and d["r"] in ("replicates", "contradicts", "new", "extends")]
key.sort(key=lambda d: (d["dt"]), reverse=True)
CUR = {(s.get("url") or "").lower() for f in F.values() for s in (f.get("studies") or []) + (f.get("contradicting") or [])}
TXT = "".join(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", n), encoding="utf8").read().lower() for n in ("treatments.json", "trials.json") if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", n)))
key = [d for d in key if not (d["r"] == "new" and (d["u"].lower() in CUR or d["u"].lower().split("doi.org/")[-1] in TXT or f'{(d["a"] or "x").split()[0]} {d["y"]}'.lower() in TXT))]
grp = {k: [d for d in key if d["r"] == k] for k in ("contradicts", "new", "replicates", "extends")}
keyhtml = "".join(f'<details class="kg"{" open" if k in ("contradicts", "new") else ""}><summary><b>{ {"contradicts": "Results that contradict an atlas finding", "new": "Findings the atlas does not cover yet", "replicates": "Independent replications", "extends": "Results that extend a finding"}[k]}</b> <span class="small">{len(v)}</span></summary><ul class="ul">{"".join(item(d, False) for d in v[:80])}</ul></details>' for k, v in grp.items() if v)
newest = [d for d in pub if d["dt"] >= since and d["k"] in ("primary", "trial", "meta") and d["rv"] != "not"]
newest.sort(key=lambda d: d["dt"], reverse=True)
kc = Counter(d["k"] for d in pub)
out = f'''<!--meta {{"title":"Library: every ME/CFS paper","description":"{tot:,} ME/CFS and long COVID papers and preprints, screened against the evidence atlas: which replicate, contradict or add to what is known.","nav":"library"}}-->
<style>
.p-library .tiles{{display:grid;grid-template-columns:repeat(4,1fr);gap:1rem;margin:1.2rem 0 1.6rem;font-family:var(--f-sans)}}
.p-library .tiles div{{border-top:1.5px solid var(--ink);padding-top:.4rem}} .p-library .tiles b{{display:block;font-size:1.4rem;font-variant-numeric:lining-nums}} .p-library .tiles span{{font-size:.8rem;color:var(--muted)}}
@media (max-width:640px){{.p-library .tiles{{grid-template-columns:1fr 1fr}}}}
.p-library .ul{{list-style:none;padding:0;margin:.4rem 0 1rem}}
.p-library .ul li{{border-top:1px solid var(--rule);padding:.55rem 0;font-size:.92rem;line-height:1.4}}
.p-library .ul li>a{{text-decoration:none;color:var(--ink);font-weight:600}} .p-library .ul li>a:hover{{text-decoration:underline}}
.p-library .m{{display:block;font-family:var(--f-sans);font-size:.76rem;color:var(--muted);margin-top:.1rem}}
.p-library .fd{{display:block;font-size:.88rem;margin-top:.25rem}}
.p-library .fa{{font-family:var(--f-sans);font-size:.72rem;padding:.05rem .35rem;border:1px solid var(--rule);border-radius:2px;text-decoration:none;white-space:nowrap}}
.p-library .rl{{font-family:var(--f-sans);font-size:.7rem;font-weight:600;text-transform:uppercase;letter-spacing:.03em;margin-right:.2rem}}
.p-library .r-contradicts{{color:var(--ref)}} .p-library .r-replicates{{color:var(--rep)}} .p-library .r-new{{color:var(--warm)}} .p-library .r-extends{{color:var(--sev)}}
.p-library .kg{{border-top:1.5px solid var(--ink);padding:.6rem 0}} .p-library .kg summary{{cursor:pointer;font-family:var(--f-sans)}}
.p-library .flt{{display:flex;flex-wrap:wrap;gap:.5rem;margin:1rem 0 .4rem;font-family:var(--f-sans);font-size:.85rem;position:sticky;top:0;background:var(--bg);padding:.5rem 0;z-index:2}}
.p-library .flt input,.p-library .flt select{{font:inherit;padding:.35rem .5rem;border:1px solid var(--rule);background:var(--bg);color:var(--ink)}} .p-library .flt input{{flex:1 1 14rem}}
.p-library #more{{font:inherit;font-family:var(--f-sans);padding:.4rem .8rem;border:1px solid var(--rule);background:none;color:var(--link);cursor:pointer}}
</style>
<div class="wrap">
<h1>Library</h1>
<p class="deck">Every ME/CFS paper and preprint we can find, plus long COVID papers on exertion, fatigue, autonomic and cognitive symptoms. Each primary study is read and marked against the <a href="/atlas">evidence atlas</a>: does it replicate, contradict or add to what is known?</p>
<p class="byline">Updated {{{{built}}}} from Europe PMC (PubMed, PMC, medRxiv, bioRxiv, Research Square) · <a href="/library-data.json">Download as JSON</a></p>
<div class="tiles"><div><b>{tot:,}</b><span>papers and preprints</span></div><div><b>{core:,}</b><span>primary studies, trials and meta-analyses</span></div><div><b>{cit + scr:,}</b><span>read and placed against the atlas</span></div><div><b>{pre:,}</b><span>preprints</span></div></div>

<h2 id="key">Papers that change the picture</h2>
<p>Large or well-designed studies that contradict, replicate or add to an atlas finding. Each links to the finding it bears on.</p>
{keyhtml or '<p class="note">Screening in progress.</p>'}

<h2 id="new">Newest studies</h2>
<p>Primary studies, trials and meta-analyses published in the last three months.</p>
<ul class="ul">{"".join(item(d) for d in newest[:40])}</ul>

<h2 id="all">Search everything</h2>
<div class="flt"><input id="q" type="search" placeholder="Search title, finding, author or journal" aria-label="Search">
<select id="fp" aria-label="Population"><option value="">All populations</option><option>ME/CFS</option><option>long COVID</option><option>both</option></select>
<select id="fk" aria-label="Type"><option value="">All types</option>{"".join(f'<option value="{k}">{e(v)} ({kc[k]})</option>' for k, v in KL.items() if kc[k])}</select>
<select id="fd" aria-label="Topic"><option value="">All topics</option>{"".join(f'<option value="{k}">{e(v)}</option>' for k, v in DL.items())}</select>
<select id="fr" aria-label="Bearing on the atlas"><option value="">Any bearing</option>{"".join(f'<option value="{k}">{e(v)}</option>' for k, v in RL.items())}</select>
<select id="fs" aria-label="Sort"><option value="d">Newest first</option><option value="i">Most important first</option></select></div>
<p class="note" id="cnt"></p>
<ul class="ul" id="list"></ul>
<p><button id="more" type="button">Show more</button></p>
<p class="note">Screening reads the title and abstract only. "Important" means a large or well-designed study that should move an atlas finding. Missing a paper? <a href="mailto:hello@solvingmecfs.org">Tell us</a>.</p>
</div>
<script>
(function(){{var D=[],lim=50,R={json.dumps(RL)},A={json.dumps({k: (F[k].get("subdomain") or k) for k in F})};
function esc(s){{return String(s).replace(/[&<>"]/g,function(c){{return{{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}}[c]}})}}
function li(d){{return '<li>'+(d.r&&d.r!=='none'?'<span class="rl r-'+d.r+'">'+R[d.r]+'</span> ':'')+'<a href="'+esc(d.u)+'" rel="noopener">'+esc(d.t)+'</a><span class="m">'+esc(d.a)+' · '+esc(d.j)+' '+esc(d.y)+(d.n?' · '+esc(d.n):'')+(d.pp&&d.j!=='Preprint'?' · preprint':'')+'</span>'+(d.f?'<span class="fd">'+esc(d.f)+d.at.map(function(x){{return ' <a class="fa" href="/atlas#f-'+esc(x)+'">'+esc(A[x]||x)+'</a>'}}).join('')+'</span>':'')+'</li>'}}
function g(id){{return document.getElementById(id).value}}
function draw(){{var q=g('q').toLowerCase(),p=g('fp'),k=g('fk'),dm=g('fd'),r=g('fr'),s=g('fs');
var f=D.filter(function(d){{return(!p||d.p===p)&&(!k||d.k===k)&&(!dm||d.dm.indexOf(dm)>=0)&&(!r||d.r===r)&&(!q||(d.t+' '+d.f+' '+d.a+' '+d.j).toLowerCase().indexOf(q)>=0)}});
f.sort(s==='i'?function(a,b){{return(b.im-a.im)||(b.dt<a.dt?-1:1)}}:function(a,b){{return a.dt<b.dt?1:-1}});
document.getElementById('cnt').textContent=f.length.toLocaleString()+' papers';
document.getElementById('list').innerHTML=f.slice(0,lim).map(li).join('');
document.getElementById('more').style.display=f.length>lim?'':'none'}}
['q','fp','fk','fd','fr','fs'].forEach(function(id){{document.getElementById(id).addEventListener('input',function(){{lim=50;draw()}})}});
document.getElementById('more').addEventListener('click',function(){{lim+=100;draw()}});
fetch('/library-data.json').then(function(r){{return r.json()}}).then(function(d){{D=d;draw()}});}})();
</script>
'''
open(os.path.join(SRC, "pages", "library.html"), "w", encoding="utf8").write(out)
print("library", tot, "screened", len(S), "key", len(key), "newest", len(newest))
