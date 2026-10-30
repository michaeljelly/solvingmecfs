import json, sys, re
base = sys.argv[1]
T = json.load(open(f"{base}/work/trials_skeleton.json"))
P = {p["pid"]: p for p in json.load(open(f"{base}/work/pubs_selected.json"))}
CT = {}
for f in ("ct_mecfs_flat", "ct_lc_flat", "ct_mecfsterm_flat", "ct_lcterm_flat"):
    for r in json.load(open(f"{base}/work/{f}.json")): CT[r["nct"]] = r
def excerpt(ab, n=650):
    m = re.search(r"\b(Results?|Findings|RESULTS?)\b[:\s]", ab)
    s = ab[m.start():] if m else ab[-n:]
    s = re.sub(r"\s+", " ", s)
    if len(s) > n:
        c = re.search(r"\b(Conclusions?|Interpretation|CONCLUSIONS?|Discussion)\b", s)
        if c and c.start() > 0:
            head = s[:max(0, n - (len(s) - c.start()) )] if len(s) - c.start() < n else ""
            s = (s[:n - min(n, len(s) - c.start()) - 3] + " … " if len(s) - c.start() < n else "") + s[c.start():c.start() + n]
        else: s = s[:n]
    return s
mode = sys.argv[2]
out = []
for k, t in sorted(T.items()):
    if mode == "pubs" and t["pubs"]:
        r = CT.get(t["nct"]) if t["nct"] else None
        hdr = f"## {k}" + (f" | {r['brief'][:80]} | {r['status']} | n={r['enrollment']} | res={'Y' if r['has_results'] else 'N'}" if r else "")
        out.append(hdr)
        for pid in t["pubs"]:
            p = P[pid]
            out.append(f"  {pid} {p['year']} {p['title'][:140]}")
            out.append(f"    > {excerpt(p['abstract'])}")
    if mode == "ctres" and not t["pubs"] and t["nct"] and CT[t["nct"]]["has_results"]:
        r = CT[t["nct"]]
        out.append(f"## {k} | {r['brief'][:90]} | {r['status']} | n={r.get('results_n') or r['enrollment']} | {r['allocation']} {r['masking']}")
        for o in r.get("results_primary", [])[:2]:
            out.append(f"  PO: {o['title'][:120]} [{o['unit']}]")
            for a in o["analyses"][:3]: out.append(f"    AN: p={a[0]} {a[1]} {a[2]}={a[3]} CI[{a[4]},{a[5]}] {a[6][:80]}")
            out.append("    M: " + "; ".join(f"{m[0][:20]}/{m[1][:25]}={m[2]}±{m[3]}" for m in o["meas"][:6]))
open(f"{base}/work/review_{mode}.txt", "w").write("\n".join(out))
print(len(out))
