#!/usr/bin/env python3
"""Render the linked parts of the site from data/graph.json: the step, arrow and constraint sections of pages/model.html,
the open questions on pages/questions.html and the proposed experiments on pages/index.html. Each generated block sits
between <!--gen:NAME--> and <!--/gen:NAME--> markers in the page fragment; everything outside the markers is hand-written.

usage: python3 -I make_graph.py <site_src_dir> <atlas_research_dir>
"""
import glob, html, json, os, re, sys

SRC, RES = sys.argv[1], sys.argv[2]
e = lambda s: html.escape(str(s or ""), quote=True)
G = json.load(open(os.path.join(SRC, "data", "graph.json"), encoding="utf8"))
T = json.load(open(os.path.join(SRC, "data", "treatments.json"), encoding="utf8"))["treatments"]
F = {}
for p in glob.glob(os.path.join(RES, "*.json")):
    for f in json.load(open(p, encoding="utf8")).get("findings", []):
        F[f["id"]] = f
H = {h["id"]: h for h in json.load(open(os.path.join(RES, "hypotheses.json"), encoding="utf8"))["hypotheses"]}
STEP = {s["n"]: s for s in G["steps"]}
Q = {q["id"]: q for q in G["questions"]}
QNUM = {q["id"]: i + 1 for i, q in enumerate(G["questions"])}
X = {x["id"]: x for x in G["experiments"]}
XNUM = {x["id"]: i + 1 for i, x in enumerate(G["experiments"])}

GCLS = {"replicated": "rep", "several": "sev", "early": "early", "lc-only": "early", "contested": "cont", "inferred": "cont", "refuted": "ref"}
GLAB = {**G["grades"], "refuted": "Refuted"}
EDGE = {"supported": ("rep", "Some support"), "partly": ("early", "Partly supported"), "untested": ("cont", "Untested")}
FIT = {"good": ("rep", "Fit: good"), "partial": ("early", "Fit: partial"), "poor": ("cont", "Fit: poor"), "none": ("ref", "Not explained")}
VCLS = {"try": "rep", "doctor": "sev", "specialist": "early", "trial-only": "cont", "avoid": "ref"}

def pill(grade):
    return f'<span class="pill {GCLS.get(grade, "cont")}">{e(GLAB.get(grade, grade))}</span>'

def short(text, n=125):
    text = re.split(r"(?<=[a-z0-9)%])\. ", text)[0].rstrip(".")
    return text if len(text) <= n else text[:n].rsplit(" ", 1)[0].rstrip(",;:") + "…"

def flink(fid):
    f = F[fid]
    return f'<a href="/atlas#f-{e(fid)}">{e(short(f["finding"]))}</a> {pill(f["grade"])}'

def steplink(n, name=True):
    return f'<a href="/model#step-{n}">{n}{" · " + e(STEP[n]["name"]) if name else ""}</a>'

def qlink(qid):
    return f'<a href="/questions#q-{e(qid)}">{e(Q[qid]["title"])}</a>'

def xlink(xid):
    return f'<a href="/#e-{e(xid)}">{XNUM[xid]}. {e(X[xid]["title"])}</a>'

def hlink(hid):
    return f'<a href="/atlas#h-{e(hid)}">{e(H[hid]["name"])}</a>'

def tlink(t):
    return f'<a class="{VCLS[t["verdict"]]}" href="/treatments#{e(t["id"])}">{e(t["name"])}</a>'

def joined(items):
    return ", ".join(items)

# ---------- model page
def model_steps():
    out = []
    for s in G["steps"]:
        n = s["n"]
        qs = [q["id"] for q in G["questions"] if n in q["steps"]]
        xs = [x["id"] for x in G["experiments"] if n in x["steps"] and len(x["steps"]) < 8]
        ts = [t for t in T if n in (t.get("links") or [])]
        order = list(VCLS)
        ts.sort(key=lambda t: order.index(t["verdict"]))
        ev = "".join(f"<li>{flink(i)}</li>" for i in s["evidence"])
        weak = "".join(f"<li>{flink(i)}</li>" for i in s.get("weak", []) + s.get("tension", []))
        links = []
        if qs: links.append(f'<div><b>Open questions</b>{joined(qlink(q) for q in qs)}</div>')
        if xs: links.append(f'<div><b>Proposed experiments</b>{joined(xlink(x) for x in xs)}</div>')
        if s.get("hypotheses"): links.append(f'<div><b>Related hypotheses in the atlas</b>{joined(hlink(h) for h in s["hypotheses"])}</div>')
        if ts: links.append(f'<div><b>Treatments that act here, coloured by verdict</b>{joined(tlink(t) for t in ts)}</div>')
        out.append(f'''<section class="step" id="step-{n}">
<h3><span class="sn">{n}</span>{e(s["name"])} {pill(s["found"])}</h3>
<p class="claim">{e(s["claim"])}</p>
<p class="fnote">{e(s["found_note"])}</p>
<div class="sev2"><div><p class="sublabel">Evidence</p><ul>{ev}</ul></div>{f'<div><p class="sublabel">Weaker, or against</p><ul>{weak}</ul></div>' if weak else ""}</div>
<p class="kill"><b>Would be disproved if:</b> {e(s["disproof"])}</p>
<div class="xl">{"".join(links)}</div>
</section>''')
    return "\n".join(out)

def model_edges():
    rows = []
    for ed in G["edges"]:
        a = "Infection" if ed["from"] == "infection" else steplink(ed["from"])
        cls, lab = EDGE[ed["status"]]
        rows.append(f'<tr><td class="num">{a} → {steplink(ed["to"])}</td><td><span class="pill {cls}">{lab}</span></td><td>{e(ed["note"])}</td></tr>')
    return f'<div class="tbl"><table class="edges"><thead><tr><th>Arrow</th><th>Status</th><th>Evidence</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>'

def model_constraints():
    rows = []
    for c in G["constraints"]:
        cls, lab = FIT[c["fit"]]
        refs = " ".join(f'<a class="fref" href="/atlas#f-{e(i)}">{pill(F[i]["grade"]) if F[i]["grade"] != "refuted" else '<span class="pill ref">Null trial</span>'}</a>' for i in c["findings"])
        rows.append(f'<tr><td>{e(c["fact"])} {refs}</td><td><span class="pill {cls}">{lab}</span></td><td>{e(c["how"])}</td></tr>')
    return f'<div class="tbl"><table class="cons"><thead><tr><th>Observation</th><th>Fit</th><th>How the model accounts for it</th></tr></thead><tbody>{"".join(rows)}</tbody></table></div>'

# ---------- open questions
def questions():
    out = []
    for q in G["questions"]:
        n = QNUM[q["id"]]
        steps = joined(steplink(s) for s in q["steps"]) if q["steps"] else "Not explained by the model"
        exp = f' Part of <a href="/#e-{e(q["experiment"])}">proposed experiment {XNUM[q["experiment"]]}</a>.' if q.get("experiment") else ""
        fl = joined(f'<a href="/atlas#f-{e(i)}">{e(short(F[i]["finding"], 70))}</a>' for i in q["findings"])
        extra = "".join(f'<a href="{e(u)}">{e(t)}</a>' for t, u in q.get("links", []))
        out.append(f'''<div class="lead" id="q-{e(q["id"])}">
<div class="num">{n}</div>
<div class="body">
<h3>{e(q["title"])} {pill(q["grade"])}</h3>
<p>{e(q["body"])}</p>
<div class="facts">
<div><b>Still unknown</b>{e(q["unknown"])}</div>
<div><b>Experiment</b>{q["test"]}.{exp}</div>
</div>
<p class="qlinks"><b>Model steps</b> {steps} · <b>Atlas</b> {fl}</p>
{f'<div class="src">{extra}</div>' if extra else ""}
</div>
</div>''')
    small = "".join(f'<li><b>{e(s["title"])}</b>: {s["text"]} <a href="/atlas#f-{e(s["finding"])}">Atlas</a></li>' for s in G["small_questions"])
    out.append(f'<details><summary>Smaller open leads ({len(G["small_questions"])})</summary><ul>{small}</ul></details>')
    legend = "".join(f'<span>{pill(k)} {d}</span>' for k, d in [("replicated", "seen in independent cohorts"), ("several", "more than one study"), ("early", "one or two studies"), ("contested", "conflicting results")])
    return f'<div class="legend">{legend}</div>\n' + "\n".join(out)

# ---------- proposed experiments
def experiments():
    out = []
    for x in G["experiments"]:
        steps = "all steps" if len(x["steps"]) == 8 else joined(steplink(s, False) for s in x["steps"])
        qs = joined(qlink(q) for q in x["questions"])
        hy = " · ".join(hlink(h) for h in x["hypotheses"])
        out.append(f'''<li id="e-{e(x["id"])}"><h3>{e(x["title"])}</h3>
<p>{e(x["body"])}</p>
<p>{e(x["decides"])}</p>
<p class="dec">Tests model steps {steps} · Answers: {qs}{f" · Separates: {hy}" if hy else ""}</p></li>''')
    return '<ol class="exps">' + "\n".join(out) + "</ol>"

GEN = {"model.html": {"constraints": model_constraints, "steps": model_steps, "edges": model_edges},
       "questions.html": {"questions": questions}, "index.html": {"experiments": experiments}}
for page, blocks in GEN.items():
    p = os.path.join(SRC, "pages", page)
    s = open(p, encoding="utf8").read()
    for name, fn in blocks.items():
        pat = re.compile(rf"(<!--gen:{name}-->).*?(<!--/gen:{name}-->)", re.S)
        if not pat.search(s):
            sys.exit(f"{page}: missing <!--gen:{name}--> markers")
        s = pat.sub(lambda m: m.group(1) + "\n" + fn() + "\n" + m.group(2), s)
    open(p, "w", encoding="utf8").write(s)
    print("graph:", page, ", ".join(blocks))
