#!/usr/bin/env python3
"""Build site/atlas.html (static) and site/atlas-data.json from per-domain research JSON files.

usage: python3 -I build_atlas.py <research_dir> <site_dir> <site_src_dir>
Writes site/atlas-data.json and the page fragment site-src/pages/atlas.html; site_build.py wraps it.
"""
import html, json, os, re, sys, datetime

RES, SITE, SRC = sys.argv[1], sys.argv[2], sys.argv[3]
ORDER = ["vascular", "muscle-exertion", "neuro", "immune", "infection", "omics", "endocrine", "gut", "genetics", "trials"]
GRADES = {
    "replicated": ("go", "Replicated", "Found in several independent cohorts"),
    "several": ("go2", "Several studies", "More than one study, not yet independently settled"),
    "early": ("maybe", "Early", "One or two studies, usually small"),
    "contested": ("contest", "Contested", "Conflicting results or method doubts"),
    "refuted": ("no", "Refuted", "Well-designed work found no effect"),
}
STATUS = {"active": "go", "weakened": "maybe", "largely refuted": "no", "speculative": "contest"}

def e(s):
    return html.escape(str(s or ""), quote=True)

REPO_TREE = "https://github.com/michaeljelly/solvingmecfs/tree/main/analyses/"
def link(cite, url):
    url = (url or "").strip()
    if (cite or "").startswith("Solving ME/CFS re-analysis") and url in ("", "/#data", "/data", "https://solvingmecfs.org/#data"):
        hit = next((cid for pat, cid in REANALYSIS if re.search(pat, cite, re.I)), None)
        url = REPO_TREE + (hit.removeprefix("analysis-") if hit else "")
    if url and not url.startswith(("http", "/")):
        url = "https://doi.org/" + url.removeprefix("doi:").strip() if re.match(r"^(doi:)?10\.", url) else ""
    if url == "/#data": url = "/data"
    return f'<a href="{e(url)}" rel="noopener">{e(cite)}</a>' if url else e(cite)

domains = []
for d in ORDER:
    p = os.path.join(RES, d + ".json")
    if os.path.exists(p):
        domains.append(json.load(open(p, encoding="utf8")))
hyp_path = os.path.join(RES, "hypotheses.json")
hyps = json.load(open(hyp_path, encoding="utf8"))["hypotheses"] if os.path.exists(hyp_path) else []
hyp_names = {h["id"]: h["name"] for h in hyps}

# where each finding and hypothesis sits in the working model, open questions and experiments (site-src/data/graph.json)
GRAPH = json.load(open(os.path.join(SRC, "data", "graph.json"), encoding="utf8"))
ROLE = {"evidence": "evidence for", "weak": "weak evidence for", "tension": "evidence against"}
in_model, in_q, h_steps = {}, {}, {}
for s_ in GRAPH["steps"]:
    for k, lab in ROLE.items():
        for fid in s_.get(k, []):
            in_model.setdefault(fid, []).append((s_["n"], s_["name"], lab))
    for hid in s_.get("hypotheses", []):
        h_steps.setdefault(hid, []).append((s_["n"], s_["name"]))
for q_ in GRAPH["questions"]:
    for fid in q_["findings"]:
        in_q.setdefault(fid, []).append(q_)

def graph_line(fid):
    parts = [f'<a href="/model#step-{n}">step {n}, {e(nm)}</a> ({lab})' for n, nm, lab in in_model.get(fid, [])]
    qs = [f'<a href="/questions#q-{e(q_["id"])}">{e(q_["title"])}</a>' for q_ in in_q.get(fid, [])]
    out = []
    if parts: out.append("Working model: " + "; ".join(parts))
    if qs: out.append("Open question: " + "; ".join(qs))
    return f'<p class="small inmodel">{" · ".join(out)}</p>' if out else ""

# merge domain-level hypothesis notes into the register
notes = {}
for d in domains:
    for n in d.get("hypotheses_notes", []):
        notes.setdefault(n["id"], []).append((d["title"], n))

all_findings, all_studies, glossary, gaps = [], [], {}, []
for d in domains:
    for f in d.get("findings", []):
        f["_domain"] = d["domain"]; f["_dtitle"] = d["title"]
        all_findings.append(f)
        for kind in ("studies", "contradicting"):
            for s in f.get(kind, []) or []:
                all_studies.append({**s, "finding": f["id"], "ftext": f["finding"], "domain": d["title"], "kind": "supports" if kind == "studies" else "contradicts"})
    for g in d.get("glossary", []):
        glossary.setdefault(g["term"].strip(), g["plain"])
    for g in d.get("gaps", []):
        gaps.append((d["title"], g))

# link studies to the data catalogue (site/data/catalogue.json), if present
CAT = []
cp = os.path.join(SITE, "data", "catalogue.json")
if os.path.exists(cp):
    CAT = json.load(open(cp, encoding="utf8"))
REANALYSIS = [(r"biomap", "analysis-biomapai-reanalysis"), (r"proteom", "analysis-proteome"), (r"metabolome", "analysis-aromatic-lactates"),
    (r"microbiome|species meta", "analysis-microbiome-ipa"), (r"genetics_mr|mendelian|\bMR\b", "analysis-mr-proteome"),
    (r"single-cell|methylation", "analysis-singlecell-methylation"), (r"epi|finngen|finregistry", "analysis-epi-literature"),
    (r"\bPEM\b", "analysis-pem-germain2022"),
    (r"decodeme|HLA|MVP", "analysis-decodeme-genetics"), (r"datasets", "analysis-round1-replication")]
by_doi, by_key = {}, {}
for c in CAT:
    if c.get("doi"):
        by_doi.setdefault(c["doi"].lower().removeprefix("doi:").strip(), []).append(c["id"])
    if c["id"].startswith(("github-", "analysis-")): continue
    cit = c.get("citation") or ""
    m, y = re.match(r"([A-Z][A-Za-z'\-]+)[ ,]", cit), re.search(r"\b(19|20)\d\d\b", cit) or re.search(r"(19|20)\d\d", c["id"])
    if m and y and m.group(1) not in ("The", "National", "University", "Long", "Passive", "Zenodo", "GEO", "NCBI", "NIH", "UK", "U"):
        by_key.setdefault((m.group(1).lower(), y.group(0)), []).append(c["id"])
cat_ids = {c["id"] for c in CAT}
cat_doi = {c["id"]: c.get("doi", "").removeprefix("doi:").strip() for c in CAT}

def datasets_for(s):
    cite, url = s.get("cite") or "", (s.get("url") or "").lower()
    if "solving me/cfs" in cite.lower():
        for pat, cid in REANALYSIS:
            if cid in cat_ids and re.search(pat, cite, re.I):
                return [cid]
        return []
    out = []
    for doi, ids in by_doi.items():
        if doi and doi in url: out += ids
    m, y = re.match(r"([A-Z][A-Za-z'\-]+)", cite), re.search(r"\b(19|20)\d\d\b", cite)
    sdoi = re.search(r"10\.\d{4,}/[^\s?#]+", url)
    if m and y:
        for cid in by_key.get((m.group(1).lower(), y.group(0)), []):
            cdoi = (cat_doi.get(cid) or "").lower()
            if sdoi and cdoi and cdoi not in sdoi.group(0) and sdoi.group(0) not in cdoi: continue  # different paper, same author-year
            out.append(cid)
    out = list(dict.fromkeys(out))
    if any(not i.startswith("geo-") for i in out): out = [i for i in out if not i.startswith("geo-")]
    return out[:3]

def dlinks(s):
    ids = datasets_for(s)
    return "".join(f' <a class="dtag" href="/data#{e(i)}" title="Dataset in our catalogue: {e(i)}">{"data" if len(ids) == 1 else e(i)}</a>' for i in ids)

for st in all_studies:
    st["datasets"] = datasets_for(st)

today = datetime.date.today().isoformat()
counts = {k: sum(1 for f in all_findings if f.get("grade") == k) for k in GRADES}

def pill(grade):
    cls, label, _ = GRADES.get(grade, ("contest", grade or "?", ""))
    return f'<span class="pill {cls}">{e(label)}</span>'

def studies_table(rows, title):
    if not rows: return ""
    out = [f'<p class="sublabel">{title}</p><div class="tbl"><table><thead><tr><th>Study</th><th>n</th><th>Design</th><th>Result</th></tr></thead><tbody>']
    for s in rows:
        v = "" if s.get("verified") else ' <span class="unv" title="Not re-checked against the paper this session">unverified</span>'
        out.append(f'<tr><td>{link(s.get("cite"), s.get("url"))}{v}{dlinks(s)}</td><td class="num">{e(s.get("n"))}</td><td>{e(s.get("design"))}</td><td>{e(s.get("result"))}</td></tr>')
    out.append("</tbody></table></div>")
    return "".join(out)

def finding_html(f):
    hy = " ".join(f'<a class="tag" href="#h-{e(h)}">{e(hyp_names.get(h, h))}</a>' for h in f.get("hypotheses", []) or [])
    oq = "".join(f"<li>{e(q)}</li>" for q in f.get("open_questions", []) or [])
    search = " ".join([f.get("finding",""), f.get("plain",""), f.get("subdomain",""), f["_dtitle"]] + [s.get("cite","") for s in f.get("studies",[]) or []]).lower()
    return f'''<details class="finding" id="f-{e(f["id"])}" data-grade="{e(f.get("grade"))}" data-search="{e(search)}">
<summary><span class="ftext">{e(f["finding"])}</span> {pill(f.get("grade"))}</summary>
<div class="fbody">
<p class="plain">{e(f.get("plain"))}</p>
<p>{e(f.get("detail"))}</p>
{studies_table(f.get("studies"), "Evidence")}
{studies_table(f.get("contradicting"), "Against or not replicated")}
{f'<p class="sublabel">Open questions</p><ul>{oq}</ul>' if oq else ""}
{graph_line(f["id"])}
<p class="small">{f'Bears on: {hy}' if hy else ""}{' · ' if hy and f.get("test_available") else ""}{f'Test: {e(f.get("test_available"))}' if f.get("test_available") else ""}</p>
</div></details>'''

body = []
body.append(f'''<div class="wrap">

<h1>ME/CFS evidence atlas</h1>
<p class="deck">Published findings in ME/CFS and related post-infectious illness, including long COVID, organised by body system, graded by replication, and linked to the hypotheses they bear on.</p>
<p class="byline">Rebuilt {today} · {len(all_findings)} findings, {len(all_studies)} study entries, {len(hyps)} hypotheses · <a href="/atlas-data.json">Download as JSON</a></p>

<div class="box key">
<p>For an overview, read the summaries by body system, then the hypotheses. To find a specific result, use the search box, the grade filters or the <a href="#studies">study index</a>.</p>
<p class="grades">{" ".join(f'{pill(k)} <span class="small">{e(v[2])} ({counts[k]})</span>' for k, v in GRADES.items())}</p>
</div>

<h2 id="overview">Overview by body system</h2>
<div class="cards">''')
for d in domains:
    n = len(d.get("findings", []))
    rep = sum(1 for f in d.get("findings", []) if f.get("grade") in ("replicated", "several"))
    body.append(f'<a class="card" href="#d-{e(d["domain"])}"><h3>{e(d["title"])}</h3><p>{e(d.get("overview"))}</p><p class="small">{n} findings · {rep} replicated or multi-study</p></a>')
body.append("</div>")

# hypotheses
if hyps:
    body.append('<h2 id="hypotheses">Hypotheses</h2><p>The main explanations proposed for ME/CFS, with the evidence for and against each. Several can be true at once, and many describe different links in the same chain. Tap one to open it.</p>')
    for h in hyps:
        ef = "".join(f"<li>{e(x)}</li>" for x in h.get("evidence_for", []))
        ea = "".join(f"<li>{e(x)}</li>" for x in h.get("evidence_against", []))
        pr = "".join(f"<li>{e(x)}</li>" for x in h.get("predictions", []))
        dt = "".join(f"<li>{e(x)}</li>" for x in h.get("decisive_tests", []))
        dn = "".join(f'<li><strong>{e(t)}:</strong> {"for: " + e(n.get("for")) if n.get("for") else ""}{"; " if n.get("for") and n.get("against") else ""}{"against: " + e(n.get("against")) if n.get("against") else ""}</li>' for t, n in notes.get(h["id"], []))
        rel = " ".join(f'<a class="tag" href="#h-{e(r)}">{e(hyp_names.get(r, r))}</a>' for r in h.get("related", []))
        linked = [f for f in all_findings if h["id"] in (f.get("hypotheses") or [])]
        lf = "".join(f'<li><a href="#f-{e(f["id"])}">{e(f["finding"])}</a> {pill(f.get("grade"))}</li>' for f in linked[:25])
        refs = "; ".join(link(r.get("cite"), r.get("url")) for r in h.get("key_refs", []))
        body.append(f'''<details class="finding hyp" id="h-{e(h["id"])}"><summary><span class="ftext">{e(h["name"])}<span class="oneline">{e(h.get("one_line"))}</span></span> <span class="pill {STATUS.get(h.get("status"),"contest")}">{e(h.get("status"))}</span></summary><div class="fbody">
<p class="plain">{e(h.get("claim"))}</p>
{f'<p class="small inmodel">In the working model: {", ".join(f"<a href=/model#step-{n}>step {n}, {e(nm)}</a>" for n, nm in h_steps[h["id"]])}</p>' if h["id"] in h_steps else ""}
<p class="small">Proponents: {e(", ".join(h.get("proponents", [])))}</p>
<div class="cols"><div><p class="sublabel">Evidence for</p><ul>{ef}</ul></div><div><p class="sublabel">Evidence against</p><ul>{ea}</ul></div></div>
{f'<p class="sublabel">Predicts</p><ul>{pr}</ul>' if pr else ""}
{f'<p class="sublabel">Tests that would settle it</p><ul>{dt}</ul>' if dt else ""}
{f'<p class="sublabel">Notes from each body system</p><ul>{dn}</ul>' if dn else ""}
{f'<p class="sublabel">Linked findings</p><ul>{lf}</ul>' if lf else ""}
<p class="small">{f"Related: {rel}" if rel else ""}{"<br>" if rel and refs else ""}{f"Key references: {refs}" if refs else ""}</p>
</div></details>''')

# findings by domain
body.append('''<h2 id="findings">Findings by body system</h2>
<div class="filters"><label for="q" class="sr">Search findings</label><input id="q" type="search" placeholder="Search findings, studies, authors" autocomplete="off">
<div class="gbtns" role="group" aria-label="Filter by grade"><button type="button" data-g="" class="on">All</button>''' + "".join(f'<button type="button" data-g="{k}">{v[1]}</button>' for k, v in GRADES.items()) + '''</div><p class="small" id="countline"></p></div>''')
for d in domains:
    body.append(f'<section class="domain" id="d-{e(d["domain"])}"><h3 class="dtitle">{e(d["title"])}</h3><p class="small">{e(d.get("overview"))}</p>')
    subs = {}
    for f in d.get("findings", []):
        subs.setdefault(f.get("subdomain") or "General", []).append(f)
    for s, fs in subs.items():
        body.append(f'<p class="sublabel subd">{e(s)}</p>' + "".join(finding_html(f) for f in fs))
    kr = "".join(f'<li>{link(r.get("cite"), r.get("url"))}: {e(r.get("why"))}</li>' for r in d.get("key_reviews", []))
    if kr: body.append(f'<p class="sublabel">Key reviews</p><ul class="sources">{kr}</ul>')
    body.append("</section>")

# deepest blood studies and what they miss
comp = [c for d in domains for c in d.get("comprehensive_studies", [])]
if comp:
    body.append('<h2 id="blood-coverage">The deepest blood studies, and what they miss</h2><p>A null result only rules out what was measured, in the people sampled, at the time sampled.</p><div class="tbl"><table><thead><tr><th>Study</th><th>People</th><th>Compartments</th><th>Time points</th><th>Severity</th><th>Live-cell tests</th><th>Main limits</th></tr></thead><tbody>')
    for c in comp:
        body.append("<tr>" + "".join(f"<td>{e(c.get(k))}</td>" for k in ("name", "n", "compartments", "timepoints", "severity", "fresh_cell_function", "limits")) + "</tr>")
    body.append("</tbody></table></div>")

# gaps
if gaps:
    body.append('<h2 id="gaps">Gaps nobody has closed</h2><ul>' + "".join(f'<li><strong>{e(t)}:</strong> {e(g)}</li>' for t, g in gaps) + "</ul>")

# study index
body.append('<h2 id="studies">Study index</h2><p class="small">Every study cited above. Tap a column to sort.</p><div class="tbl"><table id="stbl"><thead><tr><th data-k="0">Study</th><th data-k="1">n</th><th data-k="2">Design</th><th data-k="3">Area</th><th data-k="4">Finding</th><th data-k="5">Direction</th></tr></thead><tbody>')
seen = set()
for s in sorted(all_studies, key=lambda s: s.get("cite","")):
    key = (s.get("cite"), s["finding"], s["kind"])
    if key in seen: continue
    seen.add(key)
    body.append(f'<tr><td>{link(s.get("cite"), s.get("url"))}{dlinks(s)}</td><td class="num">{e(s.get("n"))}</td><td>{e(s.get("design"))}</td><td>{e(s["domain"])}</td><td><a href="#f-{e(s["finding"])}">{e(s["ftext"][:90])}</a></td><td>{e(s["kind"])}</td></tr>')
body.append("</tbody></table></div>")

# glossary
if glossary:
    body.append('<h2 id="glossary">Glossary</h2><div class="tbl"><table><tbody>' + "".join(f'<tr><td><strong>{e(t)}</strong></td><td>{e(p)}</td></tr>' for t, p in sorted(glossary.items(), key=lambda x: x[0].lower())) + "</tbody></table></div>")

body.append(f'''<h2 id="about">About this atlas</h2>
<p>Compiled by Michael Jelly, with AI agents assisting, from published papers, preprints and re-analyses of public datasets. Grades are a judgement about replication; tell us where one is wrong. Findings used by the <a href="/model">working model</a> link to the step they support or contradict. The underlying data is in <a href="/atlas-data.json">atlas-data.json</a>.</p>
</div>''')

script = r'''<script>
(function(){
var q=document.getElementById('q'),btns=document.querySelectorAll('.gbtns button'),g='',cl=document.getElementById('countline');
var items=[].slice.call(document.querySelectorAll('.domain details.finding'));
function run(){var t=(q.value||'').toLowerCase().trim(),n=0;items.forEach(function(d){var ok=(!g||d.dataset.grade===g)&&(!t||d.dataset.search.indexOf(t)>-1);d.hidden=!ok;if(ok)n++;if(t&&ok)d.open=true;});
document.querySelectorAll('.domain').forEach(function(s){s.hidden=!s.querySelector('details.finding:not([hidden])');});cl.textContent=n+' of '+items.length+' findings shown';}
q.addEventListener('input',run);btns.forEach(function(b){b.addEventListener('click',function(){btns.forEach(function(x){x.classList.remove('on')});b.classList.add('on');g=b.dataset.g;run();});});
function openHash(){var h=location.hash.slice(1);if(!h)return;var el=document.getElementById(h);if(el&&el.tagName==='DETAILS'){el.open=true;}}
window.addEventListener('hashchange',openHash);openHash();run();
var tb=document.getElementById('stbl');if(tb){var dir={};tb.querySelectorAll('th').forEach(function(th){th.style.cursor='pointer';th.addEventListener('click',function(){var k=+th.dataset.k,rows=[].slice.call(tb.tBodies[0].rows);dir[k]=!dir[k];rows.sort(function(a,b){var x=a.cells[k].textContent,y=b.cells[k].textContent,nx=parseFloat(x),ny=parseFloat(y);var r=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y);return dir[k]?r:-r;});rows.forEach(function(r){tb.tBodies[0].appendChild(r)});});});}
})();
</script>'''

extra_css = '''
.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(15rem,1fr));gap:.4rem 1.4rem;margin:1rem 0}
.cards .card{display:block;border-top:1px solid var(--rule);padding:.6rem 0 .2rem;color:var(--ink);text-decoration:none;min-width:0}
.cards .card:hover h3{text-decoration:underline}
.cards .card h3{margin:0 0 .3rem;font-size:.95rem}.cards .card p{margin:.25rem 0;font-size:.88rem;line-height:1.45}
details.finding{border-top:1px solid var(--rule);margin:0}
details.finding>summary{cursor:pointer;padding:.6rem 0;display:flex;gap:.8rem;justify-content:space-between;align-items:flex-start;list-style:none}
details.finding>summary::-webkit-details-marker{display:none}
details.finding>summary .ftext{min-width:0}
.fbody{padding:0 0 .9rem;min-width:0}
.plain{font-size:1rem}
.sublabel{font-family:var(--f-sans);font-weight:600;font-size:.78rem;color:var(--muted);margin:1rem 0 .2rem}
.subd{margin-top:1.4rem;font-size:.95rem;color:var(--muted)}
.dtitle{font-size:1.2rem;margin-top:2.4rem}
.tag{font-family:var(--f-sans);font-size:.8rem;color:var(--link);margin-right:.5rem}
.unv{font-family:var(--f-sans);font-size:.7rem;color:var(--muted);font-style:italic;white-space:nowrap}
.filters{position:sticky;top:env(safe-area-inset-top,0px);background:var(--bg);padding:.6rem 0;z-index:5;border-bottom:1px solid var(--rule);font-family:var(--f-sans)}
.filters input{width:100%;font:inherit;font-size:.9rem;padding:.5rem .6rem;border:1px solid var(--rule);border-radius:3px;background:var(--bg);color:var(--ink)}
.gbtns{display:flex;flex-wrap:wrap;gap:.35rem;margin-top:.5rem}
.gbtns button{font:inherit;font-size:.78rem;padding:.2rem .55rem;border:1px solid var(--rule);border-radius:3px;background:none;color:var(--ink);cursor:pointer}
.gbtns button.on{background:var(--ink);color:var(--bg);border-color:var(--ink)}
.grades{display:flex;flex-wrap:wrap;gap:.3rem .8rem;align-items:center}
.oneline{display:block;font-size:.9rem;color:var(--muted);margin-top:.15rem;font-weight:400}
.dtag{font-family:var(--f-sans);font-size:.72rem;color:var(--link);border:1px solid currentColor;border-radius:2px;padding:0 .25rem;text-decoration:none;white-space:nowrap}
.hyp>summary .ftext{font-family:var(--f-sans);font-weight:600;font-size:.95rem}
.inmodel{margin:.6rem 0 .2rem;padding-left:.6rem;border-left:2px solid var(--accent)}
.sr{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0)}
'''
desc = f"A living review of ME/CFS research: {len(all_findings)} findings by body system graded by replication, {len(hyps)} hypotheses with evidence for and against, and an index of {len(all_studies)} studies."
frag = "<!--meta " + json.dumps({"title": "ME/CFS evidence atlas", "description": desc, "nav": "atlas"}) + "-->\n<style>" + extra_css + "</style>\n" + "\n".join(body) + script + "\n"
open(os.path.join(SRC, "pages", "atlas.html"), "w", encoding="utf8").write(frag)
data = {"built": today, "domains": [{k: v for k, v in d.items()} for d in domains], "hypotheses": hyps}
for d in data["domains"]:
    for f in d.get("findings", []):
        f.pop("_domain", None); f.pop("_dtitle", None)
json.dump(data, open(os.path.join(SITE, "atlas-data.json"), "w", encoding="utf8"), ensure_ascii=False, indent=1)
print(f"built: {len(domains)} domains, {len(all_findings)} findings, {len(all_studies)} studies, {len(hyps)} hypotheses")
