#!/usr/bin/env python3
"""Write pages/treatments.html from data/treatments.json. usage: python3 -I make_treatments.py <site_src_dir>"""
import html, json, os, sys
SRC = sys.argv[1]
D = json.load(open(os.path.join(SRC, "data", "treatments.json"), encoding="utf8"))
e = lambda s: html.escape(str(s or ""), quote=True)
V = {"try": ("rep", "Reasonable to try", "Some controlled evidence of benefit or a direct fit to the model, low risk, low cost, no prescription needed."),
     "doctor": ("sev", "With a doctor", "Needs a GP prescription or test."),
     "specialist": ("early", "Specialist", "Needs a specialist, an infusion or a clinic."),
     "trial-only": ("cont", "Trial only", "Evidence of benefit too thin, or risk or cost too high, outside a trial."),
     "avoid": ("ref", "Avoid", "Evidence of harm, or a reliable null result.")}
EV = {"rct-positive": "Positive RCT", "rct-mixed": "Mixed RCTs", "rct-negative": "Negative RCT", "open-label": "Open-label only",
      "case-series": "Case series or surveys", "mechanistic": "Mechanism only", "none": "No direct evidence"}
GRAPH = json.load(open(os.path.join(SRC, "data", "graph.json"), encoding="utf8"))
LINKS = {s["n"]: s["name"] for s in GRAPH["steps"]}
T = D["treatments"]
counts = {k: sum(1 for t in T if t["verdict"] == k) for k in V}

def card(t):
    cls, lab, _ = V[t["verdict"]]
    links = "".join(f'<a class="lk" href="/model#step-{n}" title="{e(LINKS.get(n))}">{n}</a>' for n in t.get("links") or [])
    refs = "; ".join((f'<a href="{e(r["url"])}" rel="noopener">{e(r["cite"])}</a>' if r.get("url") else e(r["cite"])) + ("" if r.get("verified") else ' <span class="unv">unverified</span>') for r in t.get("refs") or [])
    rows = [("Why it might help, in the model", t.get("model_rationale")), ("What it does", t.get("what_it_does")), ("Who it might suit", t.get("who_it_might_suit")), ("Risks", t.get("risks")),
            ("Severe patients", t.get("severe_note")), ("Access and cost", " · ".join(x for x in [t.get("access"), t.get("cost")] if x)), ("Pending", t.get("pending"))]
    body = "".join(f'<div><b>{k}</b><p>{e(v)}</p></div>' for k, v in rows if v)
    search = " ".join([t["name"], t.get("model_rationale", ""), t.get("what_it_does", ""), t.get("evidence_line", ""), t.get("who_it_might_suit", "")]).lower()
    return f'''<details class="tx" id="{e(t["id"])}" data-v="{e(t["verdict"])}" data-s="{e(search)}">
<summary><span class="tn">{e(t["name"])}</span><span class="pill {cls}">{lab}</span>
<span class="tl">{e(t.get("evidence_line"))}</span>
<span class="tm"><span class="ev">{e(EV.get(t.get("evidence"), t.get("evidence")))}{(" · " + e(t["evidence_caveat"])) if t.get("evidence_caveat") else ""}</span>{f'<span class="lks">Steps {links}</span>' if links else ""}</span></summary>
<div class="txb">{body}{f'<p class="refs"><b>Sources</b> {refs}</p>' if refs else ""}</div></details>'''

board = []
for g in D["groups"]:
    items = [t for t in T if t["group"] == g["id"]]
    if not items: continue
    chips = "".join(f'<a class="chip {V[t["verdict"]][0]}" href="#{e(t["id"])}">{e(t["name"])}</a>' for t in items)
    board.append(f'<div class="brow"><span class="bg">{e(g["title"])}</span><span class="bchips">{chips}</span></div>')

sections = []
for g in D["groups"]:
    items = [t for t in T if t["group"] == g["id"]]
    if not items: continue
    order = list(V)
    items.sort(key=lambda t: order.index(t["verdict"]))
    sections.append(f'<section class="tg" id="g-{e(g["id"])}"><h2>{e(g["title"])}</h2><p class="gi">{e(g.get("intro"))}</p>{"".join(card(t) for t in items)}</section>')

legend = "".join(f'<button type="button" data-v="{k}" class="vb"><span class="pill {c}">{l}</span> <span class="n">{counts[k]}</span></button>' for k, (c, l, _) in V.items())
defs = "".join(f'<li><span class="pill {c}">{l}</span> {d}</li>' for k, (c, l, d) in V.items())
meta = {"title": "Treatments people are trying for ME/CFS", "description": f"{len(T)} treatments people with ME/CFS are trying, graded by mechanism, evidence and risk, with notes for severe patients.", "nav": "treatments"}
page = f'''<!--meta {json.dumps(meta)}-->
<style>
.p-treatments .wrap{{max-width:52rem}}
.rules2{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,20rem),1fr));gap:.8rem;margin:1.2rem 0}}
.rules2 div{{border-top:1px solid var(--rule);padding:.6rem 0 0;font-size:.92rem}}
.rules2 b{{display:block;font-family:var(--f-sans);font-size:.88rem;margin-bottom:.15rem}}
.board{{display:grid;gap:.55rem;margin:1rem 0 .5rem}}
.brow{{display:grid;grid-template-columns:10rem 1fr;gap:.8rem;align-items:start}}
.bg{{font-family:var(--f-sans);font-weight:600;font-size:.8rem;color:var(--ink);padding-top:.1rem}}
.bchips{{display:flex;flex-wrap:wrap;gap:.3rem}}
.chip{{font-family:var(--f-sans);font-size:.8rem;font-weight:500;padding:.05rem 0;margin-right:.6rem;border-bottom:2px solid currentColor;text-decoration:none;max-width:100%;overflow-wrap:anywhere}}
.chip:hover{{outline:1.5px solid currentColor}}
@media (max-width:560px){{.brow{{grid-template-columns:1fr;gap:.2rem}}}}
.vdefs{{list-style:none;padding:0;font-size:.88rem;color:var(--ink-2);display:grid;gap:.25rem}}
.tfilter{{position:sticky;top:env(safe-area-inset-top,0px);z-index:5;background:var(--bg);padding:.6rem 0;border-bottom:1px solid var(--rule);margin-top:1.5rem}}
.tfilter input{{width:100%;font:inherit;padding:.55rem .7rem;border:1px solid var(--rule);border-radius:6px;background:var(--surface);color:var(--ink)}}
.vbs{{display:flex;flex-wrap:wrap;gap:.35rem;margin-top:.5rem}}
.vb{{font:inherit;background:none;border:1px solid var(--rule);border-radius:3px;padding:.2rem .5rem;cursor:pointer;display:inline-flex;align-items:center;gap:.25rem}}
.vb .n{{font-family:var(--f-mono);font-size:.72rem;color:var(--muted)}}
.vb.off{{opacity:.35}}
.tg h2{{margin-top:2.4rem}}
.gi{{color:var(--ink-2);font-size:.95rem}}
details.tx{{border-top:1px solid var(--rule);margin:0}}
details.tx>summary{{list-style:none;padding:.75rem 0;display:grid;grid-template-columns:1fr auto;gap:.25rem .8rem;align-items:start}}
details.tx>summary .pill{{align-self:start;justify-self:end}}
details.tx>summary::-webkit-details-marker{{display:none}}

.tn{{font-family:var(--f-sans);font-weight:600;font-size:.98rem}}
.tl{{grid-column:1/-1;font-size:.9rem;color:var(--ink-2);line-height:1.45}}
.tm{{grid-column:1/-1;display:flex;flex-wrap:wrap;gap:.3rem .9rem;font-family:var(--f-sans);font-size:.76rem;color:var(--muted);align-items:center}}
.lk{{display:inline-grid;place-items:center;width:1.2rem;height:1.2rem;border-radius:50%;border:1px solid var(--faint);color:var(--ink);text-decoration:none;margin-left:.15rem;font-size:.7rem;font-weight:600}}
.txb{{padding:0 0 .9rem;display:grid;gap:.1rem}}
.txb b{{font-family:var(--f-sans);font-size:.76rem;color:var(--muted);font-weight:600}}
.txb p{{margin:.1rem 0 .55rem;font-size:.95rem}}
.refs{{font-size:.84rem!important;word-break:break-word}}
.unv{{font-family:var(--f-sans);font-size:.7rem;color:var(--muted);font-style:italic}}
details.tx:target{{outline:2px solid var(--accent);outline-offset:2px}}
</style>
<div class="wrap">
<h1>Treatments people are trying</h1>
<p class="deck">No drug is approved for ME/CFS, and many patients try treatments that are plausible but unproven. This page lists {len(T)} of them. Each entry says which step of the <a href="/model">working model</a> it acts on and why the model predicts it could help, then grades the evidence that it does, the risks and the cost.</p>
<p class="byline">Updated {e(D.get("updated"))}</p>
<div class="box warn"><p>Not medical advice. Prescription items need a doctor, and severe patients can be harmed by things that are safe for others.</p></div>
<h2 id="how">Trying a treatment</h2>
<div class="rules2">
<div><b>One change at a time</b>Start new things 1–2 weeks apart, so if something makes you worse you know what to drop.</div>
<div><b>Track the night</b>Overnight resting heart rate and HRV from a wearable, plus a one-minute 0–10 daily score.</div>
<div><b>Set a stop rule first</b>If resting heart rate stays 5+ bpm higher, or HRV 20% lower, for three nights after a start, pause the newest thing.</div>
<div><b>Pacing first</b>Repeated crashes make the effect of anything else impossible to judge, and may worsen the illness.</div>
</div>
<h2 id="board">Summary</h2>
<ul class="vdefs">{defs}</ul>
<div class="board">{"".join(board)}</div>
<div class="tfilter"><input id="tq" type="search" placeholder="Search treatments, mechanisms, symptoms" aria-label="Search treatments"><div class="vbs">{legend}</div></div>
{"".join(sections)}
<p class="small" style="margin-top:2rem">Numbers in each entry are steps in the <a href="/model">working model</a> the treatment acts on. The verdict weighs evidence of benefit against risk and cost: low risk alone is not enough for "Reasonable to try". Sources marked unverified have not yet been checked against the paper. Verdicts are updated as trials report; see the <a href="/method#predictions">predictions</a> for pending trials.</p>
</div>
<script>(function(){{var q=document.getElementById('tq'),off={{}},bs=document.querySelectorAll('.vb');
function run(){{var s=q.value.trim().toLowerCase();document.querySelectorAll('details.tx').forEach(function(d){{d.hidden=!!off[d.dataset.v]||(s&&d.dataset.s.indexOf(s)<0&&d.querySelector('.tn').textContent.toLowerCase().indexOf(s)<0)}});document.querySelectorAll('.tg').forEach(function(g){{g.hidden=!g.querySelector('details.tx:not([hidden])')}})}}
q.addEventListener('input',run);bs.forEach(function(b){{b.addEventListener('click',function(){{off[b.dataset.v]=!off[b.dataset.v];b.classList.toggle('off',!!off[b.dataset.v]);run()}})}});
function openHash(){{var t=location.hash&&document.getElementById(location.hash.slice(1));if(t&&t.tagName==='DETAILS')t.open=true}}addEventListener('hashchange',openHash);openHash();}})()</script>
'''
open(os.path.join(SRC, "pages", "treatments.html"), "w", encoding="utf8").write(page)
print("treatments:", len(T), counts)
