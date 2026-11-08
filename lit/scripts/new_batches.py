#!/usr/bin/env python3
"""Write batches of papers not yet screened (uncited primary/trial/meta, any year) to lit/batches/n<date>-NN.jsonl.
usage: python3 -I new_batches.py   -> prints the batch files written (none if nothing new)"""
import json, os, glob, datetime
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = json.load(open(os.path.join(H, "library.json"), encoding="utf8"))
done = set()
for f in glob.glob(os.path.join(H, "screen", "*.json")):
    done |= {str(x.get("id")) for x in json.load(open(f, encoding="utf8")) if isinstance(x, dict)}
todo = [r for r in L if r["kind"] in ("primary", "trial", "meta") and not r["cited"] and str(r["id"]) not in done]
tag = datetime.date.today().strftime("n%Y%m%d")
for i in range(0, len(todo), 250):
    p = os.path.join(H, "batches", f"{tag}-{i // 250 + 1:02d}.jsonl")
    with open(p, "w", encoding="utf8") as o:
        for r in todo[i:i + 250]:
            o.write(json.dumps({k: r.get(k) for k in ("id", "title", "year", "journal", "pop", "kind")} | {"abstract": r.get("abstract") or ""}, ensure_ascii=False) + "\n")
    print(p, len(todo[i:i + 250]))
print("unscreened", len(todo))
