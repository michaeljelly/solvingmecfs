#!/usr/bin/env python3
"""Wrap every page fragment in site-src/pages/ in the shared head, header and footer, and write it to site/.

usage: python3 -I site_build.py <site_src_dir> <site_dir>
A fragment starts with <!--meta {json}--> (title, description, nav, optional head), then optional <style>, then body HTML.
"""
import hashlib, json, os, re, sys, datetime, shutil

SRC, SITE = sys.argv[1], sys.argv[2]
NAV = [("primer", "/primer", "Primer"), ("model", "/model", "Working model"), ("atlas", "/atlas", "Evidence"), ("questions", "/questions", "Open questions"), ("research", "/research", "Projects"), ("trials", "/trials", "Trials"),
       ("treatments", "/treatments", "Treatments"), ("funding", "/funding", "Funding"), ("data", "/data", "Data"), ("method", "/method", "Method")]
LOGO = '<svg viewBox="0 0 24 24" aria-hidden="true"><circle cx="12" cy="12" r="8.5" fill="none" stroke="var(--accent)" stroke-width="2.6" stroke-dasharray="44 9.4" stroke-linecap="round" transform="rotate(-62 12 12)"/><circle cx="18.9" cy="7.2" r="2.1" fill="var(--warm)"/></svg>'
FAVICON = "data:image/svg+xml," + ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><circle cx="12" cy="12" r="8.5" fill="none" stroke="%230d6a6f" stroke-width="2.6" stroke-dasharray="44 9.4" stroke-linecap="round" transform="rotate(-62 12 12)"/><circle cx="18.9" cy="7.2" r="2.1" fill="%23b4532a"/></svg>').replace("<", "%3C").replace(">", "%3E").replace('"', "'")
SUN = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 1.5a6.5 6.5 0 1 0 0 13 6.5 6.5 0 0 0 0-13Zm0 1.6v9.8a4.9 4.9 0 0 1 0-9.8Z" fill="currentColor"/></svg>'
FONTS = "https://fonts.googleapis.com/css2?family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&family=JetBrains+Mono:wght@400;600&display=swap"
THEME_JS = "<script>(function(){var r=document.documentElement;try{var s=localStorage.getItem('theme');if(s)r.setAttribute('data-theme',s)}catch(e){}window.__toggleTheme=function(){var d=r.getAttribute('data-theme')||(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light');var n=d==='dark'?'light':'dark';r.setAttribute('data-theme',n);try{localStorage.setItem('theme',n)}catch(e){}}})()</script>"

os.makedirs(os.path.join(SITE, "assets"), exist_ok=True)
css = open(os.path.join(SRC, "site.css"), encoding="utf8").read()
ver = hashlib.sha1(css.encode()).hexdigest()[:8]
open(os.path.join(SITE, "assets", "site.css"), "w", encoding="utf8").write(css)
shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), "og", "og.png"), os.path.join(SITE, "assets", "og.png"))
today = datetime.date.today()
built = f"{today.day} {today.strftime('%b %Y')}"

def svg(name):
    p = os.path.join(SRC, "svg", name)
    if not os.path.exists(p): return f"<!-- missing {name} -->"
    t = open(p, encoding="utf8").read()
    t = re.sub(r"<metadata.*?</metadata>", "", t, flags=re.S)
    return re.sub(r'\s+xmlns:c2pa="[^"]*"', "", t).strip()

def header(active):
    links = "".join(f'<a href="{h}"{" aria-current=\"page\"" if k == active else ""}>{t}</a>' for k, h, t in NAV)
    return f'''<header class="site-head"><div class="wide">
<a class="brand" href="/">Solving ME/CFS</a>
<nav class="site-nav" aria-label="Site">{links}</nav>
<button class="theme-btn" type="button" onclick="__toggleTheme()" aria-label="Switch light or dark theme">{SUN}</button>
</div></header>
<script>(function(){{var n=document.querySelector('.site-nav'),a=n&&n.querySelector('[aria-current]');if(a&&a.offsetLeft+a.offsetWidth>n.clientWidth)n.scrollLeft=a.offsetLeft-24}})()</script>'''

REPO = os.environ.get("SITE_REPO", "https://github.com/michaeljelly/solvingmecfs")
FOOT = f'''<footer class="site-foot"><div class="l-body">
<p>Solving ME/CFS is written by <a href="https://michaeljelly.com">Michael Jelly</a>, with AI agents assisting with literature review and analysis. <a href="mailto:hello@solvingmecfs.org">hello@solvingmecfs.org</a>{f' · <a href="{REPO}">Source and analysis code</a>' if REPO else ""} · <a href="/updates">What's new</a> (<a href="/feed.xml">feed</a>) · Updated {built}</p>
</div></footer>'''

# live numbers and the evidence map, from the atlas and catalogue data
STATS, EVMAP = {}, ""
ap, cp = os.path.join(SITE, "atlas-data.json"), os.path.join(SITE, "data", "catalogue.json")
if os.path.exists(ap):
    A = json.load(open(ap, encoding="utf8"))
    fs = [(d, f) for d in A["domains"] for f in d["findings"]]
    STATS["n_findings"] = len(fs)
    STATS["n_studies"] = sum(len(f.get("studies") or []) + len(f.get("contradicting") or []) for _, f in fs)
    STATS["n_hyp"] = len(A["hypotheses"])
    STATS["n_domains"] = len(A["domains"])
    STATS["n_replicated"] = sum(1 for _, f in fs if f.get("grade") == "replicated")
    G = [("replicated", "rep", "Replicated"), ("several", "sev", "Several studies"), ("early", "early", "Early"), ("contested", "cont", "Contested"), ("refuted", "ref", "Refuted")]
    SHORT = {"vascular": "Blood flow and vessels", "muscle-exertion": "Muscle and exertion", "neuro": "Brain and nerves", "immune": "Immune system", "infection": "Infection triggers", "omics": "Blood omics", "endocrine": "Hormones", "gut": "Gut", "genetics": "Genetics", "trials": "Treatment trials"}
    mx = max(len(d["findings"]) for d in A["domains"])
    rows = []
    for d in A["domains"]:
        c = {g: sum(1 for f in d["findings"] if f.get("grade") == g) for g, _, _ in G}
        segs = "".join(f'<span class="em-s em-{cls}" style="flex:{c[g]}" title="{c[g]} {lab.lower()}">{c[g]}</span>' for g, cls, lab in G if c[g])
        rows.append(f'<a class="em-row" href="/atlas#d-{d["domain"]}"><span class="em-name">{SHORT.get(d["domain"], d["title"])}</span><span class="em-bar" style="width:{100*len(d["findings"])/mx:.1f}%">{segs}</span><span class="em-n">{len(d["findings"])}</span></a>')
    leg = "".join(f'<span><i class="em-{cls}"></i>{lab}</span>' for _, cls, lab in G)
    EVMAP = f'<div class="evmap" role="figure" aria-label="Findings per body system by evidence grade">{"".join(rows)}<div class="em-leg">{leg}</div></div>'
tp, fp = os.path.join(SRC, "data", "trials.json"), os.path.join(SRC, "data", "funding.json")
if os.path.exists(tp): STATS["n_trials"] = len(json.load(open(tp, encoding="utf8"))["trials"])
if os.path.exists(fp): STATS["n_calls"] = len([c for c in json.load(open(fp, encoding="utf8"))["calls"] if c["id"] != "amar"])
if os.path.exists(cp):
    C = json.load(open(cp, encoding="utf8"))
    STATS["n_datasets"] = len(C)
    STATS["n_rehosted"] = sum(1 for c in C if c.get("status") == "rehosted")
    STATS["n_analyses"] = sum(1 for c in C if c["id"].startswith("analysis-"))

pages = sorted(f for f in os.listdir(os.path.join(SRC, "pages")) if f.endswith(".html"))
for f in pages:
    raw = open(os.path.join(SRC, "pages", f), encoding="utf8").read()
    m = re.match(r"\s*<!--meta (\{.*?\})-->\s*", raw, re.S)
    meta = json.loads(m.group(1)); body = raw[m.end():]
    body = re.sub(r"\{\{svg:([\w.\-]+)\}\}", lambda mm: svg(mm.group(1)), body)
    body = body.replace("{{built}}", built).replace("{{evidence_map}}", EVMAP)
    body = re.sub(r"\{\{(n_\w+)\}\}", lambda mm: f"{STATS.get(mm.group(1), '?'):,}" if isinstance(STATS.get(mm.group(1)), int) else "?", body)
    used = set(re.findall(r'id="([^"]+)"', body))
    def _h2id(mm):
        slug = re.sub(r"[^a-z0-9]+", "-", re.sub(r"<[^>]+>", "", mm.group(2)).lower()).strip("-")[:40]
        if not slug or slug in used: return mm.group(0)
        used.add(slug); return f'<h2{mm.group(1)} id="{slug}">{mm.group(2)}</h2>'
    body = re.sub(r'<h2((?:\s+class="[^"]*")?)>(.*?)</h2>', _h2id, body)
    styles = "".join(re.findall(r"<style>.*?</style>", body, re.S))
    body = re.sub(r"<style>.*?</style>\s*", "", body, flags=re.S)
    title = meta["title"] if meta.get("nav") == "home" else f'{meta["title"]} · Solving ME/CFS'
    url = "https://solvingmecfs.org/" + ("" if f == "index.html" else f[:-5])
    out = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{title}</title>
<link rel="alternate" type="application/atom+xml" title="Solving ME/CFS: what's new" href="/feed.xml">
<meta name="description" content="{meta["description"]}">
<meta property="og:title" content="{meta["title"]}">
<meta property="og:description" content="{meta["description"]}">
<meta property="og:type" content="website">
<meta property="og:url" content="{url}">
<meta property="og:site_name" content="Solving ME/CFS">
<meta property="og:image" content="https://solvingmecfs.org/assets/og.png?v=3">
<meta property="og:image:width" content="1200">
<meta property="og:image:height" content="630">
<meta property="og:image:alt" content="ME/CFS ruins lives and has no cure. 20+ million people suffer from it. Help us solve it now.">
<meta name="twitter:card" content="summary_large_image">
<meta name="twitter:title" content="{meta["title"]}">
<meta name="twitter:description" content="{meta["description"]}">
<meta name="twitter:image" content="https://solvingmecfs.org/assets/og.png?v=3">
<link rel="canonical" href="{url}">
<link rel="icon" href="{FAVICON}">
<meta name="theme-color" content="#f5f6f4" media="(prefers-color-scheme: light)">
<meta name="theme-color" content="#101518" media="(prefers-color-scheme: dark)">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="{FONTS}">
<link rel="stylesheet" href="/assets/site.css?v={ver}">
{THEME_JS}
{styles}
{meta.get("head", "")}
</head>
<body class="p-{meta.get("nav", "page")}">
{header(meta.get("nav"))}
<main id="main">
{body.strip()}
</main>
{FOOT}
</body>
</html>
'''
    open(os.path.join(SITE, f), "w", encoding="utf8").write(out)
    print("built", f, len(out))
