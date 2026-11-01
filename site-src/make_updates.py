#!/usr/bin/env python3
"""Write pages/updates.html and <site>/feed.xml: what changed recently across projects, trials and the atlas.
usage: python3 -I make_updates.py <site_src_dir> <site_dir> <atlas_research_dir>"""
import datetime, glob, html, json, os, re, sys
SRC, SITE, RES = sys.argv[1:4]
e = lambda s: html.escape(str(s or ""), quote=True)
today = datetime.date.today()
P = json.load(open(os.path.join(SRC, "data", "projects.json"), encoding="utf8"))["projects"]
T = json.load(open(os.path.join(SRC, "data", "trials.json"), encoding="utf8"))["trials"]
F = [f | {"_d": json.load(open(p, encoding="utf8"))["title"]} for p in sorted(glob.glob(os.path.join(RES, "*.json"))) if os.path.basename(p) not in ("hypotheses.json", "overview.json")
     for f in json.load(open(p, encoding="utf8")).get("findings", [])]
since = (today - datetime.timedelta(days=183)).isoformat()
started = sorted([p for p in P if p.get("start") and since <= p["start"] <= today.isoformat() and p["focus"] == "primary"], key=lambda p: p["start"], reverse=True)
soon = (today + datetime.timedelta(days=183)).isoformat()
due = sorted([t for t in T if t.get("result") == "pending" and (t.get("primary_completion") or t.get("completion") or "") and today.isoformat()[:7] <= (t.get("primary_completion") or t.get("completion"))[:10] <= soon],
             key=lambda t: (t.get("primary_completion") or t.get("completion")))
yr = str(today.year)
recent = [f for f in F if any(re.search(rf"\b{yr}\b", s.get("cite") or "") for s in (f.get("studies") or []) + (f.get("contradicting") or []))]
newpubs = sorted([(p, t) for t in T for p in t.get("publications") or [] if str(p.get("year")) == yr and not p.get("protocol")], key=lambda x: x[1].get("n_actual") or 0, reverse=True)
GR = {"replicated": "go", "several": "go2", "early": "maybe", "contested": "contest", "refuted": "no"}
GL = {"replicated": "Replicated", "several": "Several studies", "early": "Early", "contested": "Contested", "refuted": "Refuted"}
def tlink(t): return t.get("url") or (f"https://clinicaltrials.gov/study/{t['id']}" if t["id"].startswith("NCT") else "")
li_p = "".join(f'<li><a href="{e(p["url"])}" rel="noopener">{e(p["title"])}</a><span class="m">{e(p["start"])} · {e(p["population"])} · {e(", ".join(x for x in [p.get("org"), p.get("country")] if x))}</span></li>' for p in started[:40])
li_t = "".join(f'<li><a href="{e(tlink(t))}" rel="noopener">{e(t.get("acronym") or t.get("intervention"))}</a> <span class="m">{e(t["population"])} · n={e(t.get("n_actual") or t.get("n_planned"))} · {e(t.get("design"))} · due {e((t.get("primary_completion") or t.get("completion"))[:7])}</span></li>' for t in due[:40])
li_f = "".join(f'<li><a href="/atlas#f-{e(f["id"])}"><span class="pill {GR.get(f["grade"], "contest")}">{GL.get(f["grade"], f["grade"])}</span> {e(f["finding"])}</a><span class="m">{e(f["_d"])}</span></li>' for f in recent[:60])
li_r = "".join(f'<li><a href="{"https://doi.org/" + e(p["doi"]) if p.get("doi") else "https://pubmed.ncbi.nlm.nih.gov/" + e(p.get("pmid")) + "/"}" rel="noopener">{e(p.get("title"))}</a><span class="m">{e(t.get("acronym") or t.get("intervention"))} · {e(t["population"])} · result: {e(t.get("result"))}</span></li>' for p, t in newpubs[:30])
out = f'''<!--meta {{"title":"What's new","description":"New projects, trial results due, and {yr} evidence added to the atlas.","nav":"updates"}}-->
<style>
.p-updates .ul{{list-style:none;padding:0;margin:.6rem 0 1.6rem}}
.p-updates .ul li{{border-top:1px solid var(--rule);padding:.5rem 0;font-family:var(--f-sans);font-size:.88rem;line-height:1.4}}
.p-updates .ul li a{{color:var(--ink);text-decoration:none}} .p-updates .ul li a:hover{{color:var(--link)}}
.p-updates .m{{display:block;color:var(--muted);font-size:.76rem;margin-top:.1rem}}
.p-updates .jump{{font-family:var(--f-sans);font-size:.85rem}}
</style>
<div class="wrap">
<h1>What's new</h1>
<p class="deck">New projects, trials due to report, and this year's studies, gathered from across the site.</p>
<p class="byline">Rebuilt {today.isoformat()} · <a href="/feed.xml">Subscribe (Atom feed)</a></p>
<p class="jump"><a href="#started">Projects started</a> ({len(started)}) · <a href="#due">Trials due to report</a> ({len(due)}) · <a href="#results">Trial results published in {yr}</a> ({len(newpubs)}) · <a href="#evidence">Atlas findings with {yr} studies</a> ({len(recent)})</p>
<h2 id="started">Projects started in the last six months</h2><ul class="ul">{li_p or "<li>None found.</li>"}</ul>
<h2 id="due">Trials due to finish in the next six months</h2><p>Results to watch for. Dates are the registered primary completion and often slip.</p><ul class="ul">{li_t or "<li>None found.</li>"}</ul>
<h2 id="results">Trial results published in {yr}</h2><ul class="ul">{li_r or "<li>None found.</li>"}</ul>
<h2 id="evidence">Atlas findings with {yr} studies</h2><ul class="ul">{li_f or "<li>None found.</li>"}</ul>
<p class="small">Open funding calls and deadlines are on the <a href="/funding">funding</a> page.</p>
</div>
'''
open(os.path.join(SRC, "pages", "updates.html"), "w", encoding="utf8").write(out)
# Atom feed: new projects and new trial publications
items = [(p["start"], f"New project: {p['title']}", p["url"], f"{p['population']} · {p.get('org') or ''} {p.get('country') or ''}", "p-" + p["id"]) for p in started[:50]]
items += [(f"{yr}-01-01", f"Trial result: {p.get('title')}", ("https://doi.org/" + p["doi"]) if p.get("doi") else f"https://pubmed.ncbi.nlm.nih.gov/{p.get('pmid')}/", f"{t.get('intervention')} · {t['population']} · {t.get('result')}", "t-" + str(p.get("pmid") or p.get("doi"))) for p, t in newpubs[:50]]
items.sort(key=lambda x: x[0], reverse=True)
x = lambda s: html.escape(str(s or ""), quote=False)
feed = '<?xml version="1.0" encoding="utf-8"?>\n<feed xmlns="http://www.w3.org/2005/Atom"><title>Solving ME/CFS: what\'s new</title><link href="https://solvingmecfs.org/updates"/><link rel="self" href="https://solvingmecfs.org/feed.xml"/>' \
    f'<id>https://solvingmecfs.org/feed.xml</id><updated>{today.isoformat()}T00:00:00Z</updated><author><name>Solving ME/CFS</name></author>' + \
    "".join(f'<entry><title>{x(t)}</title><link href="{x(u)}"/><id>tag:solvingmecfs.org,2026:{x(i)}</id><updated>{d[:10]}T00:00:00Z</updated><summary>{x(s)}</summary></entry>' for d, t, u, s, i in items) + "</feed>\n"
open(os.path.join(SITE, "feed.xml"), "w", encoding="utf8").write(feed)
print("updates", len(started), "started,", len(due), "due,", len(newpubs), "results,", len(recent), "findings")
