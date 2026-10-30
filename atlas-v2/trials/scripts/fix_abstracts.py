import json, glob, re, sys, html
base = sys.argv[1]
TAG = re.compile(r"</?(h\d|i|b|u|sup|sub|p|br|em|strong|span|div|abstracttext|sc|inf|list|list-item|title|sec|bold|italic)\b[^>]*>", re.I)
def clean(s): return html.unescape(TAG.sub(" ", s or "")).strip()
raw = {}
for f in glob.glob(f"{base}/raw/epmc/*.json") + glob.glob(f"{base}/raw/epmc_extra/*.json"):
    for r in json.load(open(f))["resultList"]["result"]:
        for k in (r.get("pmid"), r.get("doi"), r.get("id")):
            if k: raw[k] = r
P = json.load(open(f"{base}/work/pubs_selected.json"))
n = 0
for p in P:
    r = raw.get(p["pmid"]) or raw.get(p["doi"])
    if r:
        p["abstract"] = clean(r.get("abstractText")); p["title"] = clean(r.get("title")); n += 1
json.dump(P, open(f"{base}/work/pubs_selected.json", "w"))
print(n, len(P))
