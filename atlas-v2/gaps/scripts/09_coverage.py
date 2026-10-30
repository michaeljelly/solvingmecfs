import json, re, collections
G='/mnt/project-files/mecfs-research/atlas-v2/gaps/'
R=json.load(open(G+'reference_set.json')); A=json.load(open(G+'atlas_index.json'))
nodoi={(r.get('fid')) for r in A if not r.get('doi')}
prop=set()
for f in json.load(open(G+'new_findings.json')):
    for s in f['studies']+f['contradicting']: prop.add(s['doi'].lower())
for a in json.load(open(G+'add_studies.json')):
    for s in a['add_studies']+a['add_contradicting']: prop.add(s['doi'].lower())
OOS=re.compile(r'case report|\bmice\b|\bmouse\b|\brats?\b|traditional chinese|acupunct|herbal|hypothes|perspective|questionnaire|survey|rehabilitat|nursing|caregiver|qualitative|stigma|psychometric|validation of the|cost|economic|interview|occupational|patient experience|school|bibliometric', re.I)
cnt=collections.Counter(); per=collections.defaultdict(collections.Counter)
pref=lambda fid: ('review(key)' if fid=='KEY_REVIEW' else fid.split('-')[0]) if fid else 'unassigned'
for r in R:
    d=(r.get('doi') or '').lower(); m=r.get('match')
    fids=r.get('fids') or []
    if m=='doi' or (m and fids and any(f in nodoi for f in fids)): c='direct'
    elif d and d in prop: c='gap_proposed'
    elif OOS.search(r.get('title','')): c='out_of_scope'
    else: c='theme_or_unscreened'
    cnt[c]+=1
    dom = pref(fids[0]) if fids else (pref(r['sugg'][0][0]) if r.get('sugg') else 'unassigned')
    per[dom][c]+=1
    for s in r.get('sets',[]): per['SET:'+s.split('#')[0]][c]+=1
print(len(R), dict(cnt))
for k in sorted(per): print(k, dict(per[k]), sum(per[k].values()))
