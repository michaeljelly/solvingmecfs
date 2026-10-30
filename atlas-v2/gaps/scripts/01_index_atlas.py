import json,re,sys
d=json.load(open('/mnt/project-files/mecfs-research/site/atlas-data.json'))
rows=[]
def doi_of(u):
    if not u: return None
    m=re.search(r'(10\.\d{4,9}/[^\s?#"]+)',u)
    return m.group(1).lower().rstrip('.') if m else None
def pmid_of(u):
    m=re.search(r'pubmed\.ncbi\.nlm\.nih\.gov/(\d+)',u or '') or re.search(r'PMC?(\d+)',u or '')
    return m.group(1) if m else None
for dm in d['domains']:
    for f in dm['findings']:
        for role in ('studies','contradicting'):
            for s in f.get(role,[]) or []:
                rows.append(dict(domain=dm['domain'],fid=f['id'],grade=f.get('grade'),role=role,cite=s.get('cite',''),url=s.get('url',''),doi=doi_of(s.get('url','')),pmid=pmid_of(s.get('url','')),n=s.get('n',''),result=s.get('result','')))
    for r in dm.get('key_reviews',[]) or []:
        if isinstance(r,dict):
            rows.append(dict(domain=dm['domain'],fid='KEY_REVIEW',grade='',role='review',cite=r.get('cite',r.get('title','')),url=r.get('url',''),doi=doi_of(r.get('url','')),pmid=pmid_of(r.get('url','')),n='',result=''))
    for r in dm.get('comprehensive_studies',[]) or []:
        if isinstance(r,dict):
            rows.append(dict(domain=dm['domain'],fid='COMPREHENSIVE',grade='',role='comp',cite=r.get('cite',''),url=r.get('url',''),doi=doi_of(r.get('url','')),pmid=None,n='',result=''))
json.dump(rows,open(sys.argv[1],'w'),indent=0)
print(len(rows), sum(1 for r in rows if r['doi']), sum(1 for r in rows if r['pmid']))
# findings list
fl=[]
for dm in d['domains']:
    for f in dm['findings']:
        fl.append(dict(domain=dm['domain'],id=f['id'],sub=f.get('subdomain'),grade=f.get('grade'),finding=f['finding'],nst=len(f.get('studies',[])),ncon=len(f.get('contradicting',[]) or [])))
json.dump(fl,open(sys.argv[2],'w'),indent=0)
print(json.dumps(d['domains'][0]['key_reviews'][0])[:300])
