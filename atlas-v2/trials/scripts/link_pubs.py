import json, re, sys
from difflib import SequenceMatcher
base = sys.argv[1]
P = json.load(open(f"{base}/work/pubs_selected.json"))
CT = {}
for f in ("ct_mecfs_flat", "ct_lc_flat", "ct_mecfsterm_flat", "ct_lcterm_flat"):
    for r in json.load(open(f"{base}/work/{f}.json")): CT[r["nct"]] = r
REG = re.compile(r"\b(NCT\s?\d{8}|ISRCTN\s?\d{8}|DRKS\s?\d{8}|ACTRN\s?\d{14}|ChiCTR[-\w]{6,20}|EudraCT\s*(?:No\.?|number)?\s*:?\s*\d{4}-\d{6}-\d{2}|\d{4}-\d{6}-\d{2}|IRCT\d{6,}N\d+|CTRI/\d{4}/\d+/\d+|NTR\s?\d{3,5}|NL\s?\d{4,5}|jRCT\w+|UMIN\d{9}|KCT\d{7}|EUCTR[\w-]+|RBR-\w+|PACTR\d+|TCTR\d+|SLCTR/[\d/]+|CTIS?\s?\d{4}-\d{6}-\d{2}-\d{2}|\d{4}-\d{6}-\d{2}-\d{2})", re.I)
pm2nct = {}
for n, r in CT.items():
    for pm, typ, cit in r["refs"]:
        if pm and typ in ("RESULT", "DERIVED"): pm2nct.setdefault(pm, set()).add(n)
def norm(t): return re.sub(r"[^a-z0-9 ]", " ", (t or "").lower())
links = {}
for p in P:
    txt = p["title"] + " " + p["abstract"]
    ids = sorted({re.sub(r"\s", "", m.group(0)).upper() for m in REG.finditer(txt)})
    p["reg_ids"] = ids
    ncts = {i for i in ids if i.startswith("NCT")}
    ncts |= pm2nct.get(p["pmid"], set())
    p["nct_links"] = sorted(ncts)
json.dump(P, open(f"{base}/work/pubs_selected.json", "w"))
print(sum(1 for p in P if p["nct_links"]), sum(1 for p in P if p["reg_ids"]))
from collections import Counter
print(Counter(i[:4] for p in P for i in p["reg_ids"]))
