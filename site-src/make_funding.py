#!/usr/bin/env python3
"""Write pages/funding.html from data/funding.json. usage: python3 -I make_funding.py <site_src_dir>"""
import html, json, os, sys, datetime
SRC = sys.argv[1]
D = json.load(open(os.path.join(SRC, "data", "funding.json"), encoding="utf8"))
e = lambda s: html.escape(str(s or ""), quote=True)
C = [c for c in D["calls"] if c["id"] not in ("amar",)]
TODAY = datetime.date.today()
FIT = {"crash_window_sampling": "Crash-window sampling", "circulatory_phenotyping": "Circulatory tests", "serum_factor_tests": "Serum factor",
       "pituitary_testing": "Pituitary", "patient_tracking_software": "Tracking software", "trials": "Trials", "data_reanalysis": "Data re-analysis"}
def group(c):
    s = c["status"].lower()
    if s.startswith("open") or s.startswith("rolling"): return 0
    if s.startswith("recurring") or s.startswith("upcoming"): return 1
    return 2
GL = ["Open now or rolling", "Recurring or upcoming", "Closed, or no open call"]
def dl(c):
    d = (c.get("deadline") or {}).get("date")
    if d:
        try:
            dd = datetime.date.fromisoformat(d)
            return dd.strftime("%-d %b %Y"), dd
        except ValueError: pass
    return "", None
def ver(c):
    lv = (c.get("verification") or {}).get("level", "unverified")
    return {"verified": "Checked on funder's site", "mixed": "Partly checked", "snippet": "From search results", "unverified": "Unverified"}.get(lv, lv)
upcoming = sorted([(dl(c)[1], c) for c in C if dl(c)[1] and dl(c)[1] >= TODAY], key=lambda x: x[0])[:8]
rows_up = "".join(f'<tr><td class="num">{d.strftime("%-d %b %Y")}</td><td><a href="{e(c["url"])}" rel="noopener">{e(c["call_name"])}</a><br><span class="fn">{e(c["funder"])}</span></td><td>{e(c["amount"].get("text"))}</td></tr>' for d, c in upcoming)
def fits(c): return ", ".join(FIT.get(f, f) for f in c.get("fit") or [])
sections = []
for gi, gl in enumerate(GL):
    items = sorted([c for c in C if group(c) == gi], key=lambda c: (dl(c)[1] or datetime.date(2100, 1, 1)))
    trs = "".join(f'''<tr><td><a href="{e(c["url"])}" rel="noopener">{e(c["call_name"])}</a><br><span class="fn">{e(c["funder"])} · {e(c.get("mecfs_specificity"))}</span></td>
<td>{e(c["amount"].get("text"))}</td><td>{e((c.get("deadline") or {}).get("text"))}</td><td>{e(fits(c))}</td><td class="fn">{e(ver(c))}</td></tr>''' for c in items)
    sections.append(f'<h3>{gl} <span class="cnt">{len(items)}</span></h3><div class="tbl"><table class="calls"><thead><tr><th>Call</th><th>Amount</th><th>Deadline</th><th>Fits</th><th>Source</th></tr></thead><tbody>{trs}</tbody></table></div>')
fitmap = {}
for c in C:
    if group(c) < 2:
        for f in c.get("fit") or []: fitmap.setdefault(f, []).append(c)
fitrows = "".join(f'<tr><td>{FIT[f]}</td><td>{", ".join(f"""<a href="{e(c["url"])}" rel="noopener">{e(c["funder"].split(" (")[0])}: {e(c["call_name"])}</a>""" for c in fitmap.get(f, [])[:6]) or "None open"}</td></tr>' for f in FIT)
out = f'''<!--meta {{"title":"Funding for ME/CFS research","description":"Open and recurring funding calls for ME/CFS and long COVID research, which ones to apply to now, and which experiments each call fits.","nav":"funding"}}-->
<style>
.p-funding .fn{{font-family:var(--f-sans);font-size:.78rem;color:var(--muted)}}
.p-funding table.calls td{{vertical-align:top;font-size:.88rem}}
.p-funding h3 .cnt{{font-family:var(--f-sans);font-weight:400;color:var(--muted);font-size:.85rem}}
.p-funding ol.now li{{margin:.7rem 0}}
.p-funding .dd{{font-family:var(--f-sans);font-weight:600;font-variant-numeric:tabular-nums;color:var(--warm)}}
</style>
<div class="wrap">
<h1>Funding</h1>
<p class="deck">{len(C)} calls that can fund ME/CFS research, the ones worth applying to now, and which proposed experiment each one fits.</p>
<p class="byline">Checked {e(D.get("as_of"))}. Deadlines change; confirm on the funder's page before applying.</p>

<h2 id="now">Apply now</h2>
<ol class="now">
<li><b>Anthropic guest investigators</b> <span class="dd">about 5 Nov 2026</span>. Eight to twelve weeks working with Anthropic researchers on public data; principal investigators only. Best for <a href="/data">re-analysing existing data</a> or building an open crash-tracking tool.</li>
<li><b>Solve M.E. Ramsay grants</b> <span class="dd">13 Nov 2026</span>. $50k–150k over one to two years, for early-career researchers, researchers new to ME/CFS and collaborations, with no country limit. Best for a pilot of <a href="/#critical-path">crash-window sampling</a>, the serum-factor test or pituitary testing.</li>
<li><b>MRC Proof of Concept, stage one</b> <span class="dd">18 Nov 2026</span>. For turning a candidate marker into a diagnostic, or preparing a trial of a repurposed drug.</li>
<li><b>MRC, under its ME/CFS highlight notice</b> <span class="dd">rolling</span>. The main UK route, with decisions in June and December. MRC's first stated priority is the mechanism of post-exertional malaise.</li>
<li><b>Wellcome Discovery Award</b> <span class="dd">rolling</span>. Large awards for established researchers based in the UK.</li>
<li><b>NIH R21 and R01</b> <span class="dd">5 and 16 Feb 2027</span>. Open to non-US institutions; the R21 (about $275k over two years) suits feasibility pilots.</li>
</ol>

<h2 id="deadlines">Next deadlines</h2>
<div class="tbl"><table><thead><tr><th>Date</th><th>Call</th><th>Amount</th></tr></thead><tbody>{rows_up}</tbody></table></div>

<h2 id="fit">Which calls fit which experiment</h2>
<div class="tbl"><table><thead><tr><th>Experiment</th><th>Open or recurring calls</th></tr></thead><tbody>{fitrows}</tbody></table></div>

<h2 id="gaps">Where there is no funding</h2>
<ul>
<li>No government call specific to ME/CFS is open in the UK, US, Canada, Norway or the Netherlands.</li>
<li>Nothing is dedicated to sampling severely ill patients at home, or to pituitary and endocrine work.</li>
<li>Little supports patient-tracking software, and almost nothing is open to developers outside universities.</li>
<li>The Open Medicine Foundation, PolyBio and Simmaron fund by invitation or partnership only.</li>
</ul>

<h2 id="all">All calls</h2>
{"".join(sections)}
<p class="fn">Know of a call we have missed? Write to <a href="mailto:hello@solvingmecfs.org">hello@solvingmecfs.org</a>.</p>
</div>
'''
open(os.path.join(SRC, "pages", "funding.html"), "w", encoding="utf8").write(out)
print("built funding", len(C))
