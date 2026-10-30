# usage: python3 -I 02_fetch_epmc.py outfile max_records sort "query"
import sys,json,urllib.request,urllib.parse,time
out,maxn,sort,q=sys.argv[1],int(sys.argv[2]),sys.argv[3],sys.argv[4]
base='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'
cur='*';res=[]
while len(res)<maxn:
    p=dict(query=q,format='json',resultType='core',pageSize=min(1000,maxn-len(res)),cursorMark=cur)
    if sort!='none': p['sort']=sort
    url=base+urllib.parse.urlencode(p)
    for attempt in range(4):
        try:
            r=json.load(urllib.request.urlopen(url,timeout=120));break
        except Exception as e:
            print('retry',e,file=sys.stderr);time.sleep(5)
    rl=r.get('resultList',{}).get('result',[])
    if not rl: break
    for x in rl:
        res.append({k:x.get(k) for k in ('id','source','pmid','pmcid','doi','title','authorString','journalTitle','pubYear','citedByCount','abstractText','firstPublicationDate')}|{'pubType':(x.get('pubTypeList') or {}).get('pubType')})
    nc=r.get('nextCursorMark')
    if not nc or nc==cur: break
    cur=nc
print(out, r.get('hitCount'), len(res))
json.dump(res,open(out,'w'))
