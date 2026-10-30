import json,sys,re
r=json.load(open('/mnt/project-files/mecfs-research/atlas-v2/gaps/reference_set.json'))
pats=[l.strip() for l in open(sys.argv[1]) if l.strip()]
seen=set()
for p in pats:
    for x in r:
        if re.search(re.escape(p),x['title'] or '',re.I) and x['idx'] not in seen:
            seen.add(x['idx'])
            a=re.sub(r'<[^>]+>','',x['abstract'])
            print(f"## [{x['idx']}] {x['au'][:30]} {x['year']} {x['journal']} | doi:{x['doi']} | cites {x['cites']}\nT: {x['title']}\nA: {a[:int(sys.argv[2])]}\n")
            break
    else:
        if not any(re.search(re.escape(p),x['title'] or '',re.I) for x in r): print('NOTFOUND',p)
