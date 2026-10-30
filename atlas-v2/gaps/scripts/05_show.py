import json,sys
r=json.load(open('/mnt/project-files/mecfs-research/atlas-v2/gaps/reference_set.json'))
pref=sys.argv[1]
for x in r:
    if any(s.startswith(pref) for s in x['sets']):
        m=('M:'+x['match']+':'+','.join(x['fids'])[:60]) if x['match'] else ''
        print(f"{x['idx']}|{x['year']}|{x['cites']}|{x['au'][:18]}|{(x['title'] or '')[:115]}|{m or x['sugg'][0][0]+' '+str(x['sugg'][0][1])}")
