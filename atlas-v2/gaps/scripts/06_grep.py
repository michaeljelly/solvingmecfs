import json,sys,re
d=json.load(open('/mnt/project-files/mecfs-research/site/atlas-data.json'))
# build per-finding text
items=[]
for dm in d['domains']:
    for f in dm['findings']: items.append((f['id'],json.dumps(f).lower()))
    for k in ('key_reviews','gaps','hypotheses_notes','glossary','comprehensive_studies'):
        items.append((dm['domain']+':'+k,json.dumps(dm.get(k,'')).lower()))
items.append(('HYPOTHESES',json.dumps(d['hypotheses']).lower()))
for t in sys.argv[1:]:
    hits=[i for i,txt in items if re.search(t.lower(),txt)]
    print(f"{t}: {len(hits)} {hits[:6]}")
