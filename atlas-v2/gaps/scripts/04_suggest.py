import json
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
G='/mnt/project-files/mecfs-research/atlas-v2/gaps/'
d=json.load(open('/mnt/project-files/mecfs-research/site/atlas-data.json'))
F=[];T=[]
for dm in d['domains']:
    for f in dm['findings']:
        txt=' '.join([f['finding'],f.get('plain',''),f.get('detail',''),f.get('subdomain','')]+[s.get('result','')+' '+s.get('cite','') for s in f.get('studies',[])+(f.get('contradicting') or [])])
        F.append(f['id']);T.append(txt)
refs=json.load(open(G+'reference_set.json'))
R=[(r['title'] or '')+' '+(r['title'] or '')+' '+r['abstract'][:1500] for r in refs]
v=TfidfVectorizer(stop_words='english',ngram_range=(1,2),min_df=2,sublinear_tf=True).fit(T+R)
S=cosine_similarity(v.transform(R),v.transform(T))
for i,r in enumerate(refs):
    o=S[i].argsort()[::-1][:3]
    r['idx']=i
    r['sugg']=[(F[j],round(float(S[i,j]),2)) for j in o]
json.dump(refs,open(G+'reference_set.json','w'),indent=0)
