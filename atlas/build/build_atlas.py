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
_cg = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "site-src", "check_graph.py")
SECOND = set(re.findall(r'\("[\w-]+", "([\w-]+)"\)', open(_cg, encoding="utf8").read().split("DUPES", 1)[1].split("]", 1)[0])) if os.path.exists(_cg) else set()
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

def sweep_html(rows):
    if not rows: return ""
    c = {k: sum(1 for s in rows if s.get("relation") == k) for k in ("replicates", "contradicts", "extends")}
    lab = ", ".join(f"{v} {k[:-1] if k.endswith('s') else k}" for k, v in (("replicates", c["replicates"]), ("contradicts", c["contradicts"]), ("extends", c["extends"])) if v)
    order = {"contradicts": 0, "replicates": 1, "extends": 2}
    rs = sorted(rows, key=lambda s: (order.get(s.get("relation"), 3), -(s.get("importance") or 0)))
    tr = "".join(f'<tr><td>{link(s.get("cite"), s.get("url"))}</td><td><span class="swr sw-{e(s.get("relation"))}">{e(s.get("relation"))}</span></td><td class="num">{e(s.get("n"))}</td><td>{e(s.get("result"))}</td></tr>' for s in rs)
    return (f'<details class="sweep"><summary class="sublabel">More studies from the literature sweep ({len(rows)}: {lab})</summary>'
            f'<p class="small">Found by screening every ME/CFS paper\'s abstract against this finding. Not yet weighed into the grade; read before relying on one.</p>'
            f'<div class="tbl"><table><thead><tr><th>Study</th><th>Bearing</th><th>n</th><th>Result (from the abstract)</th></tr></thead><tbody>{tr}</tbody></table></div></details>')

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
{sweep_html(f.get("sweep"))}
{f'<p class="sublabel">Open questions</p><ul>{oq}</ul>' if oq else ""}
{graph_line(f["id"])}
<p class="small">{f'Bears on: {hy}' if hy else ""}{' · ' if hy and f.get("test_available") else ""}{f'Test: {e(f.get("test_available"))}' if f.get("test_available") else ""}</p>
</div></details>'''

body = []
body.append(f'''<div class="wrap">

<h1>ME/CFS evidence atlas</h1>
<p class="deck">Published findings in ME/CFS and related post-infectious illness, including long COVID, organised by body system, graded by replication, and linked to the hypotheses they bear on.</p>
<p class="byline">Rebuilt {today} · {len(all_findings)} findings, {len(all_studies)} study entries, {len(hyps)} hypotheses · <a href="/atlas-data.json">Download as JSON</a></p>

<h2 id="solid">What is solid, and what has not held up</h2>
<p>Findings seen in several independent cohorts are the firmest ground. Refuted ones were tested well and failed. Everything else is in between and is listed by body system below.</p>
<details class="key-list" open><summary class="sublabel">Replicated ({sum(1 for f in all_findings if f.get("grade") == "replicated" and f["id"] not in SECOND)}, counting each cross-listed finding once)</summary><ul class="kl">{"".join(f'<li><a href="#f-{e(f["id"])}">{e(f["finding"])}</a></li>' for f in all_findings if f.get("grade") == "replicated" and f["id"] not in SECOND)}</ul></details>
<details class="key-list"><summary class="sublabel">Refuted ({sum(1 for f in all_findings if f.get("grade") == "refuted" and f["id"] not in SECOND)})</summary><ul class="kl">{"".join(f'<li><a href="#f-{e(f["id"])}">{e(f["finding"])}</a></li>' for f in all_findings if f.get("grade") == "refuted" and f["id"] not in SECOND)}</ul></details>

<div class="box key">
<p>For an overview, read the summaries by body system, then the hypotheses. To find a specific result, use the search box, the grade filters or the <a href="#studies">study index</a>.</p>
<p class="grades">{" ".join(f'{pill(k)} <span class="small">{e(v[2])} ({counts[k]})</span>' for k, v in GRADES.items())}</p>
</div>

<h2 id="overview">Overview by body system</h2>
<p>Each bar shows how the findings in one body system are graded. Bar length is the number of findings; green segments are the ones that have replicated.</p>
<figure class="gfig">''')
GVAR = {"replicated": "--rep", "several": "--sev", "early": "--early", "contested": "--cont", "refuted": "--ref"}
def gcounts(fs):
    return {k: sum(1 for f in fs if f.get("grade") == k) for k in GRADES}
def gbar(c, scale=None, cls="gbar"):
    tot = sum(c.values()) or 1
    scale = scale or tot
    segs = "".join(f'<span style="width:{100 * c[k] / scale:.2f}%;background:var({GVAR[k]})" title="{c[k]} {GRADES[k][1].lower()}"></span>' for k in GRADES if c.get(k))
    return f'<span class="{cls}" role="img" aria-label="{", ".join(f"{c[k]} {GRADES[k][1].lower()}" for k in GRADES if c.get(k))}">{segs}</span>'
OV = {}
ovp = os.path.join(RES, "overview.json")
if os.path.exists(ovp): OV = json.load(open(ovp, encoding="utf8"))
fby = {f["id"]: f for f in all_findings}
_pc = os.path.join(SRC, "data", "project-counts.json")
PCOUNT = json.load(open(_pc)) if os.path.exists(_pc) else {}
mx = max(len(d.get("findings", [])) for d in domains)
for d in domains:
    c = gcounts(d.get("findings", []))
    body.append(f'<a class="grow" href="#c-{e(d["domain"])}"><span class="gname">{e(d["title"])}</span>{gbar(c, mx)}<span class="gn">{sum(c.values())}</span></a>')
body.append('<figcaption class="glegend">' + " ".join(f'<span><i style="background:var({GVAR[k]})"></i>{e(v[1])}</span>' for k, v in GRADES.items()) + '</figcaption></figure>')
body.append('<div class="sys">')
for d in domains:
    fs = d.get("findings", [])
    c = gcounts(fs)
    od = OV.get("domains", {}).get(d["domain"], {})
    groups = ""
    for gl, items in od.get("groups", []):
        lis = "".join(f'<li><a href="#f-{e(fid)}"><span class="pill {GRADES.get(fby[fid].get("grade"), ("contest",))[0]}" title="{e(GRADES.get(fby[fid].get("grade"), ("", "?"))[1])}"></span>{e(t)}</a></li>' for t, fid in items if fid in fby)
        groups += f'<p class="sublabel">{e(gl)}</p><ul class="kf">{lis}</ul>'
    take = od.get("takeaway") or (d.get("overview") or "").split(". ")[0] + "."
    body.append(f'''<section class="sysc" id="c-{e(d["domain"])}"><h3><a href="#d-{e(d["domain"])}">{e(d["title"])}</a></h3>
<p class="take">{e(take)}</p>{groups}
<details class="more"><summary>Full summary</summary><p>{e(d.get("overview"))}</p></details>
<p class="small foot"><a href="#d-{e(d["domain"])}">All {len(fs)} findings →</a>{f' · <a href="/research#t={e(d["domain"])}">{PCOUNT[d["domain"]]} active projects →</a>' if PCOUNT.get(d["domain"]) else ""}</p></section>''')
body.append("</div>")

# hypotheses
STATUS_ORDER = [("active", "Active", "Supported by current data and still being tested"),
    ("speculative", "Speculative", "Plausible, with little direct data"),
    ("weakened", "Weakened", "Key predictions have failed or not replicated"),
    ("largely refuted", "Largely refuted", "Contradicted by well-designed studies")]
VERD = OV.get("verdicts", {})
def ul(xs, cls=""):
    return f'<ul class="{cls}">' + "".join(f"<li>{e(x)}</li>" for x in xs) + "</ul>"
if hyps:
    body.append('<h2 id="hypotheses">Hypotheses</h2><p>The main explanations proposed for ME/CFS. Several can be true at once, and many describe different links in the same chain. The bar beside each shows the grades of the atlas findings linked to it, which includes findings that count against it. Tap one to open it.</p>')
    body.append('<p class="hjump">' + " · ".join(f'<a href="#hs-{k.replace(" ", "-")}">{lab} ({sum(1 for h in hyps if h.get("status") == k)})</a>' for k, lab, _ in STATUS_ORDER) + "</p>")
    hmax = max([len([f for f in all_findings if h["id"] in (f.get("hypotheses") or [])]) for h in hyps] + [1])
    for sk, slab, sdesc in STATUS_ORDER:
        group = [h for h in hyps if h.get("status") == sk]
        if not group: continue
        def strength(h):
            lf = [f for f in all_findings if h["id"] in (f.get("hypotheses") or [])]
            return (-sum(1 for f in lf if f.get("grade") in ("replicated", "several")), -len(lf))
        group.sort(key=strength)
        body.append(f'<h3 class="hstat" id="hs-{sk.replace(" ", "-")}"><span class="pill {STATUS.get(sk, "contest")}">{e(slab)}</span><span class="small">{e(sdesc)}</span></h3>')
        for h in group:
            linked = [f for f in all_findings if h["id"] in (f.get("hypotheses") or [])]
            gorder = list(GRADES)
            linked.sort(key=lambda f: gorder.index(f.get("grade")) if f.get("grade") in gorder else 9)
            c = gcounts(linked)
            ef, ea = h.get("evidence_for", []), h.get("evidence_against", [])
            bar = (gbar(c, hmax, "gbar sm") + f'<span class="gn">{len(linked)}</span>') if linked else '<span class="small nolink">umbrella model</span>'
            steps = "".join(f'<a class="chip" href="/model#step-{n}">Step {n} · {e(nm)}</a>' for n, nm in h_steps.get(h["id"], []))
            meta = []
            if steps: meta.append(f'<div><span class="ml">Working model</span><span class="chips">{steps}</span></div>')
            if h.get("proponents"): meta.append(f'<div><span class="ml">Proposed by</span><span>{e(", ".join(h["proponents"]))}</span></div>')
            vis, rest = linked[:10], linked[10:]
            def lfli(fs):
                return "".join(f'<li><a href="#f-{e(f["id"])}"><span class="pill {GRADES.get(f.get("grade"), ("contest",))[0]}"></span>{e(f["finding"])}</a><span class="ld">{e(f["_dtitle"])}</span></li>' for f in fs)
            lf = ""
            if linked:
                summ = ", ".join(f"{c[k]} {GRADES[k][1].lower()}" for k in GRADES if c[k])
                lf = (f'<div class="hsec"><p class="bh">Linked findings <span class="cnt">{len(linked)}</span></p><p class="small lsum">{gbar(c, None, "gbar")} {summ}</p><ul class="lf">{lfli(vis)}</ul>'
                      + (f'<details class="more"><summary>Show {len(rest)} more</summary><ul class="lf">{lfli(rest)}</ul></details>' if rest else "") + "</div>")
            dn = "".join(f'<li><strong>{e(t)}.</strong> {"For: " + e(n.get("for")) if n.get("for") else ""}{" " if n.get("for") and n.get("against") else ""}{"Against: " + e(n.get("against")) if n.get("against") else ""}</li>' for t, n in notes.get(h["id"], []))
            rel = "".join(f'<a class="chip" href="#h-{e(r)}">{e(hyp_names.get(r, r))}</a>' for r in h.get("related", []))
            refs = "; ".join(link(r.get("cite"), r.get("url")) for r in h.get("key_refs", []))
            pr, dt = h.get("predictions", []), h.get("decisive_tests", [])
            settle = ""
            if pr or dt:
                settle = '<div class="hsec settle">' + (f'<div><p class="bh">If true, we should see</p>{ul(pr)}</div>' if pr else "") + (f'<div><p class="bh">Tests that would settle it</p>{ul(dt, "tests")}</div>' if dt else "") + "</div>"
            body.append(f'''<details class="hyp" id="h-{e(h["id"])}"><summary><span class="hname">{e(h["name"])}<span class="oneline">{e(h.get("one_line"))}</span></span><span class="hbar">{bar}<span class="tally" title="Evidence items for and against"><b class="tf">+{len(ef)}</b> <b class="ta">−{len(ea)}</b></span></span></summary><div class="hbody">
{f'<p class="verdict"><b>Where it stands.</b> {e(VERD[h["id"]])}</p>' if h["id"] in VERD else ""}
<p class="claim">{e(h.get("claim"))}</p>
{f'<div class="hmeta">{"".join(meta)}</div>' if meta else ""}
<div class="hsec bal"><div class="for"><p class="bh">Evidence for <span class="cnt">{len(ef)}</span></p>{ul(ef)}</div><div class="ag"><p class="bh">Evidence against <span class="cnt">{len(ea)}</span></p>{ul(ea)}</div></div>
{settle}
{lf}
{f'<details class="more"><summary>Notes from each body system</summary><ul class="notes">{dn}</ul></details>' if dn else ""}
{f'<div class="hsec"><p class="bh">Related hypotheses</p><span class="chips">{rel}</span></div>' if rel else ""}
{f'<p class="small refs">Key references: {refs}</p>' if refs else ""}
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

# tests: which abnormalities can be measured, clinically or in research
def short(t, n=120):
    t = t.split(". ")[0]
    return t if len(t) <= n else t[:n].rsplit(" ", 1)[0] + "…"
TT = {"clin": [], "res": []}
for f in all_findings:
    ta = (f.get("test_available") or "").strip()
    if f.get("grade") not in ("replicated", "several", "early") or not ta: continue
    low = ta.lower()
    if low.startswith(("not applicable", "no", "n/a", "family history", "clinical history")) or f["_domain"] == "trials": continue
    k = "res" if low.startswith(("research", "private")) or "research only" in low[:40] or "not standard clinical" in low or "no validated clinical" in low else "clin"
    TT[k].append(f)
go = list(GRADES)
if TT["clin"] or TT["res"]:
    body.append(f'''<h2 id="tests">What can be measured</h2>
<p>Abnormalities graded replicated, several studies or early, with the test that measures each. {len(TT["clin"])} can be ordered through clinical services; {len(TT["res"])} need a research setting. None is diagnostic on its own: each separates groups on average, with overlap between patients and controls.</p>''')
    def ttable(fs):
        rows = "".join(f'<tr><td><a href="#f-{e(f["id"])}"><span class="pill {GRADES[f["grade"]][0]}"></span> {e(short(f["finding"]))}</a></td><td>{e(f["_dtitle"])}</td><td>{e(re.sub(r"^(Yes|Research only)[:;,.]?\\s*", "", f["test_available"]))}</td></tr>' for f in sorted(fs, key=lambda f: (go.index(f["grade"]), f["_dtitle"])))
        return f'<div class="tbl ttab"><table><thead><tr><th>Finding</th><th>Area</th><th>Test</th></tr></thead><tbody>{rows}</tbody></table></div>'
    strong = [f for f in TT["clin"] if f["grade"] in ("replicated", "several")]
    body.append(f'<h3 class="ttl">Replicated or multi-study, available clinically <span class="small">({len(strong)})</span></h3>' + ttable(strong))
    rest_c = [f for f in TT["clin"] if f not in strong]
    body.append(f'<details class="more"><summary>Early findings with a clinical test ({len(rest_c)})</summary>{ttable(rest_c)}</details>')
    body.append(f'<details class="more"><summary>Measurable only in research settings ({len(TT["res"])})</summary>{ttable(TT["res"])}</details>')

# who has been studied: severity, from the kb cohorts table (data thread)
KBC = os.path.join(RES, "..", "..", "..", "kb", "cohorts.json")
if os.path.exists(KBC):
    coh = json.load(open(KBC, encoding="utf8"))
    SK = [("n_mild", "Mild", "--sev"), ("n_moderate", "Moderate", "--early"), ("n_severe", "Severe", "--warm"), ("n_very_severe", "Very severe", "--ink")]
    rows = [c for c in coh if any(c.get(k) for k, _, _ in SK)]
    def basis(c):
        bt = c.get("severity_basis_type")
        s = (c.get("severity_scale") or "").lower()
        if bt == "self_report": return "self-report"
        if bt == "clinician_or_inclusion": return "by inclusion" if "inclusion" in s else "clinician"
        return "by inclusion" if "inclusion" in s else ("self-report" if "self" in s else "clinician")
    rows.sort(key=lambda c: -((c.get("n_severe") or 0) + (c.get("n_very_severe") or 0)))
    def dedup(cs):  # within an overlap group, count only the largest cohort
        best = {}
        for c in cs:
            k = c.get("overlap_group") or c["cohort_id"]
            v = (c.get("n_severe") or 0) + (c.get("n_very_severe") or 0)
            best[k] = max(best.get(k, 0), v)
        return sum(best.values())
    clin = dedup([c for c in rows if basis(c) != "self-report"])
    selfr = dedup([c for c in rows if basis(c) == "self-report"])
    hb = [c for c in coh if c.get("n_housebound") or c.get("n_bedbound")]
    body.append(f'''<h2 id="severity">Who has been studied</h2>
<p>About a quarter of people with ME/CFS are housebound or bedbound, but research samples are drawn mostly from people well enough to travel to a clinic. Of {len(coh)} cohorts behind the atlas, {len(rows)} report how many participants were mild, moderate or severe. Across those, {clin:,} severe or very severe participants were assessed by a clinician or recruited for severity, against {selfr:,} by self-report, most of them in DecodeME, which collected saliva only. Overlapping cohorts are counted once.</p>
<figure class="sevfig">''')
    for c in rows:
        tot = sum(c.get(k) or 0 for k, _, _ in SK) or 1
        segs = "".join(f'<span style="width:{100 * (c.get(k) or 0) / tot:.2f}%;background:var({v})" title="{c.get(k)} {l.lower()}"></span>' for k, l, v in SK if c.get(k))
        sev = (c.get("n_severe") or 0) + (c.get("n_very_severe") or 0)
        body.append(f'<div class="srow"><span class="sname">{e(c["name"].split(" (")[0])}<small>{e(basis(c))} · samples: {e((c.get("sample_types") or "none listed")[:60])}</small></span><span class="gbar">{segs}</span><span class="sn"><b>{sev:,}</b> severe<br><small>of {tot:,} graded</small></span></div>')
    body.append('<figcaption class="glegend">' + " ".join(f'<span><i style="background:var({v})"></i>{l}</span>' for _, l, v in SK) + f'</figcaption><figcaption>Each bar is one cohort, split by reported severity. Cohorts that report only housebound or bedbound counts: {"; ".join(e(c["name"].split(" (")[0]) + " (" + ", ".join(f"{c[k]} {k[2:]}" for k in ("n_housebound", "n_bedbound") if c.get(k)) + ")" for c in hb)}. Source: <a href="https://data.solvingmecfs.org/kb/cohorts.csv">cohorts table</a>, with the quote behind each count.</figcaption></figure>')

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
.kl{font-size:.92rem;line-height:1.45;padding-left:1.1rem;margin:.3rem 0 .8rem}.kl li{margin:.3rem 0}.kl a{color:var(--ink);text-decoration-color:var(--rule)}.key-list summary{cursor:pointer} .sweep{margin:.6rem 0} .sweep summary{cursor:pointer} @media (max-width:640px){.sweep table,.sweep tbody,.sweep tr,.sweep td{display:block;width:auto} .sweep thead{display:none} .sweep tr{padding:.5rem 0;border-bottom:1px solid var(--rule,#ddd)} .sweep td{border:0;padding:.1rem 0}} .swr{font-family:var(--f-sans);font-size:.7rem;text-transform:uppercase;letter-spacing:.03em;font-weight:600} .sw-contradicts{color:var(--ref)} .sw-replicates{color:var(--rep)} .sw-extends{color:var(--sev)}
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
.gfig{margin:1rem 0 1.6rem}
.grow{display:grid;grid-template-columns:13rem 1fr 2rem;gap:.8rem;align-items:center;padding:.28rem 0;color:var(--ink);text-decoration:none;font-family:var(--f-sans);font-size:.85rem}
.grow:hover .gname{text-decoration:underline}
.gbar{display:flex;height:.75rem;min-width:0;border-radius:2px;overflow:hidden;background:var(--sunk)}
.grow .gbar{background:none}
.gbar span{display:block;height:100%;min-width:2px}
.gbar span+span{border-left:1.5px solid var(--bg)}
.gbar.sm{height:.55rem;width:6rem;flex:none}
.gn{font-family:var(--f-sans);font-size:.78rem;color:var(--muted);font-variant-numeric:tabular-nums;text-align:right}
.glegend{display:flex;flex-wrap:wrap;gap:.2rem 1rem;margin-top:.6rem}
.glegend i{display:inline-block;width:.65rem;height:.65rem;border-radius:2px;margin-right:.3rem;vertical-align:-.05rem}
@media (max-width:560px){.grow{grid-template-columns:1fr 1.6rem;gap:.15rem .6rem}.grow .gname{grid-column:1/-1}}
.sys{display:grid;grid-template-columns:1fr 1fr;gap:0 2rem;margin:1rem 0 2rem}
@media (max-width:700px){.sys{grid-template-columns:1fr}}
.sysc{border-top:1.5px solid var(--ink);padding:.7rem 0 1.2rem;min-width:0;scroll-margin-top:1rem}
.sysc h3{font-size:1.02rem;margin:0 0 .35rem}.sysc h3 a{color:var(--ink);text-decoration:none}
.sysc .take{font-size:.95rem;line-height:1.45;margin:.2rem 0 .5rem}
.sysc .sublabel{margin:.7rem 0 .15rem;font-size:.72rem;text-transform:uppercase;letter-spacing:.04em}
.kf{list-style:none;padding:0;margin:0;font-family:var(--f-sans);font-size:.84rem;line-height:1.35}
.kf li{margin:.22rem 0}.kf a{color:var(--ink);text-decoration:none;display:flex;gap:.45rem;align-items:baseline}
.kf a:hover{color:var(--link)}.kf .pill{flex:none}
.sysc .foot{margin:.6rem 0 0}
details.more{margin:.6rem 0 0;font-family:var(--f-sans);font-size:.82rem}
details.more>summary{color:var(--link);list-style:none}details.more>summary::-webkit-details-marker{display:none}
details.more>summary::before{content:"▸ ";color:var(--muted)}details.more[open]>summary::before{content:"▾ "}
details.more>p{font-family:var(--f-serif);font-size:.92rem;line-height:1.5;margin:.4rem 0 0}
.hjump{font-family:var(--f-sans);font-size:.85rem}
.hstat{display:flex;flex-wrap:wrap;gap:.2rem .8rem;align-items:baseline;margin:2rem 0 .3rem;font-size:1rem}
.hstat .pill{font-size:.95rem}
details.hyp{border-top:1px solid var(--rule)}
details.hyp>summary{list-style:none;display:grid;grid-template-columns:1fr auto;gap:.3rem 1rem;align-items:center;padding:.65rem 0}
details.hyp>summary::-webkit-details-marker{display:none}
.hname{font-family:var(--f-sans);font-weight:600;font-size:.95rem;min-width:0}
.hname::before{content:"▸";display:inline-block;width:1rem;color:var(--muted);font-weight:400}
details.hyp[open] .hname::before{content:"▾"}
.hname .oneline{padding-left:1rem}
.hbar{display:flex;align-items:center;gap:.5rem;justify-self:end}
.tally{font-family:var(--f-sans);font-size:.75rem;white-space:nowrap}.tally b{font-weight:600}
.tf{color:var(--rep)}.ta{color:var(--ref)}
.nolink{font-style:italic;white-space:nowrap}
@media (max-width:560px){details.hyp>summary{grid-template-columns:1fr}.hbar{justify-self:start;padding-left:1rem}}
.hbody{padding:.2rem 0 1.4rem 1rem;border-left:2px solid var(--rule);margin-left:.25rem;min-width:0}
.verdict{background:var(--sunk);border-left:3px solid var(--ink);padding:.6rem .8rem;margin:.4rem 0 .8rem;font-size:.98rem;line-height:1.5}
.claim{font-size:.95rem;color:var(--ink-2);margin:.4rem 0}
.hmeta{display:grid;gap:.3rem;margin:.7rem 0;font-family:var(--f-sans);font-size:.82rem}
.hmeta>div{display:grid;grid-template-columns:7.5rem 1fr;gap:.6rem}
.ml{color:var(--muted);font-weight:600}
@media (max-width:560px){.hmeta>div{grid-template-columns:1fr;gap:0}}
.chips{display:flex;flex-wrap:wrap;gap:.3rem .4rem}
.chip{font-family:var(--f-sans);font-size:.78rem;border:1px solid var(--rule);border-radius:3px;padding:.05rem .4rem;text-decoration:none;color:var(--link);white-space:nowrap}
.chip:hover{border-color:var(--link)}
.hsec{margin:1rem 0 0}
.bh{font-family:var(--f-sans);font-weight:600;font-size:.8rem;margin:0 0 .3rem;color:var(--ink)}
.cnt{color:var(--muted);font-weight:400;margin-left:.2rem}
.bal{display:grid;grid-template-columns:1fr 1fr;gap:1.2rem}
.settle{display:grid;grid-template-columns:1fr 1fr;gap:1.2rem}
@media (max-width:640px){.bal,.settle{grid-template-columns:1fr;gap:.6rem}}
.bal>div,.settle>div{min-width:0;padding-top:.45rem;border-top:3px solid var(--rule)}
.bal .for{border-top-color:var(--rep)}.bal .ag{border-top-color:var(--ref)}
.hbody ul{margin:0;padding-left:0;list-style:none;font-size:.9rem;line-height:1.45}
.hbody li{margin:.35rem 0;padding-left:1.1rem;position:relative}
.hbody .bal li::before,.hbody .settle li::before{position:absolute;left:0;font-family:var(--f-sans);font-weight:700}
.bal .for li::before{content:"+";color:var(--rep)}.bal .ag li::before{content:"−";color:var(--ref)}
.settle li::before{content:"○";color:var(--muted);font-weight:400!important;font-size:.8em;top:.15em}
.settle .tests li::before{content:"→";color:var(--link)}
.lsum{display:flex;align-items:center;gap:.6rem;margin:.1rem 0 .5rem}.lsum .gbar{width:8rem;flex:none}
.hbody ul.lf{font-family:var(--f-sans);font-size:.82rem;line-height:1.35}
.hbody ul.lf li{padding-left:0;display:grid;grid-template-columns:1fr auto;gap:.8rem;margin:.3rem 0}
.lf a{color:var(--ink);text-decoration:none;display:flex;gap:.45rem;align-items:baseline}.lf a:hover{color:var(--link)}.lf .pill{flex:none}
.ld{color:var(--faint);font-size:.74rem;white-space:nowrap}
@media (max-width:560px){.hbody ul.lf li{grid-template-columns:1fr;gap:0}.ld{padding-left:1rem}}
.notes li{font-size:.86rem}
.refs{margin-top:1rem}
.sevfig{margin:1rem 0 2rem}
.srow{display:grid;grid-template-columns:15rem 1fr 6.5rem;gap:.8rem;align-items:center;padding:.35rem 0;border-bottom:1px solid var(--rule);font-family:var(--f-sans);font-size:.84rem}
.srow .gbar{height:.8rem}
.sname small,.sn small{display:block;color:var(--muted);font-size:.72rem;line-height:1.3}
.sn{text-align:right;font-variant-numeric:tabular-nums;line-height:1.25}
@media (max-width:560px){.srow{grid-template-columns:1fr 5.5rem}.srow .sname{grid-column:1/-1}}
.ttl{margin-top:1.6rem;font-size:1rem}
.ttab td{font-size:.84rem;vertical-align:top}.ttab td:first-child{min-width:16rem}.ttab td:first-child a{color:var(--ink);text-decoration:none}.ttab td:first-child a:hover{color:var(--link)}
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
