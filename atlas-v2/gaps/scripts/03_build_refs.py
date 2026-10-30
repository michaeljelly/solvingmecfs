import json,re,sys,unicodedata
G='/mnt/project-files/mecfs-research/atlas-v2/gaps/'
D=G+'downloads/'
L=lambda f: json.load(open(D+f))
atlas=json.load(open(G+'atlas_index.json'))
def norm(s): return unicodedata.normalize('NFKD',s or '').encode('ascii','ignore').decode().lower()
adoi={}
for r in atlas:
    if r['doi']: adoi.setdefault(r['doi'],[]).append(r)
aay={}
for r in atlas:
    c=norm(r['cite'])
    m=re.match(r"\s*([a-z'\- ]+?)(?: et al\.?| and | & |,| \d)",c)
    y=re.search(r'(19|20)\d\d',c)
    if m and y:
        sur=m.group(1).strip().split()[-1] if m.group(1).strip() else None
        if sur: aay.setdefault((sur,int(y.group(0))),[]).append(r)
NONBIO=re.compile(r"cognitive.behavio|\bcbt\b|graded exercise therapy|psychotherap|questionnaire|(?<!-)\bscale\b|prevalence|qualitative|experience|perspective|interview|survey|stigma|quality of life|economic|\bcosts?\b|school|nurs|general practi|primary care|criteria|case definition|psychiatr|depressive disorder|anxiety|personality|illness (perception|belief)|attribution|self-management|mindfulness|yoga|tai chi|qigong|acupunct|herbal|chinese medicine|moxibust|lightning process|occupational|social|caregiver|employment|disabilit|reliability|validity|validation|telehealth|telemedicine|telerehab|educat|knowledge|attitude|online|lived|patient-reported|symptom (cluster|burden)|clinical guideline|guideline|parents|family|children's|adolescents' |sick leave|\bwork\b|workplace|return to work|insurance|psycholog|coping|trauma|childhood adversity|nice |consensus|protocol for|study protocol|rehabilitat|counsel|hypnos|placebo response in psych|cohort profile",re.I)
def first_sur(a):
    a=norm(a).split(',')[0].strip()
    parts=a.split()
    if not parts: return ''
    # drop initials at end
    toks=[p for p in parts if not (p.isalpha() and p.isupper())]
    sur=[p for p in parts[:-1]] if len(parts)>1 else parts
    return (sur[-1] if sur else parts[0]).lower()
def match(p):
    d=(p.get('doi') or '').lower()
    if d and d in adoi: return 'doi',[x['fid'] for x in adoi[d]]
    s=first_sur(p.get('authorString','')); y=int(p.get('pubYear') or 0)
    for yy in (y,y-1,y+1):
        if (s,yy) in aay: return 'author-year',[x['fid'] for x in aay[(s,yy)]]
    return None,[]
def key(p): return (p.get('doi') or '').lower() or norm(p.get('title'))[:80]
out={}
def add(p,setname,rank=None):
    k=key(p)
    if k not in out:
        how,fids=match(p)
        out[k]=dict(sets=[],doi=p.get('doi'),pmid=p.get('pmid'),title=p.get('title'),au=(p.get('authorString') or '')[:60],year=p.get('pubYear'),journal=p.get('journalTitle'),cites=p.get('citedByCount') or 0,src=p.get('source'),match=how,fids=sorted(set(fids)),abstract=(p.get('abstractText') or '')[:3000])
    out[k]['sets'].append(setname if rank is None else f'{setname}#{rank}')
isrev=lambda p: any('review' in t.lower() for t in (p.get('pubType') or []))
# (a) reviews: top 60 cited biology + 2024-26 recent biology + top 25 LC mechanism reviews
revs=[p for p in L('a_reviews_cited.json') if not NONBIO.search(p['title'] or '')]
for i,p in enumerate(revs[:60]): add(p,'a_mecfs_review_cited',i+1)
rec=[p for p in L('a_reviews_recent.json') if not NONBIO.search(p['title'] or '') and p.get('source')=='MED']
rec.sort(key=lambda p:-(p.get('citedByCount') or 0))
for p in rec[:40]: add(p,'a_mecfs_review_recent')
lcr=[p for p in L('a_lc_reviews.json') if not NONBIO.search(p['title'] or '')]
for i,p in enumerate(lcr[:25]): add(p,'a_lc_review',i+1)
# (b) top 300 cited primary biology
prim=[p for p in L('b_primary_cited.json') if not NONBIO.search(p['title'] or '') and not isrev(p) and not re.search(r'letter|comment|editorial',' '.join(p.get('pubType') or []),re.I)]
for i,p in enumerate(prim[:300]): add(p,'b_primary_cited',i+1)
# (c) all 2024-26 primary biology
BIO=re.compile(r"blood|plasma|serum|cell|immun|cytokin|antibod|gene|genom|metabol|prote|mitochond|muscle|brain|mri|pet\b|eeg|cerebr|autonom|orthostat|tilt|heart|cardi|vascul|endotheli|microbio|gut|virus|viral|epstein|herpes|infect|cortisol|hormon|biomark|exercise|exertion|cpet|oxygen|lactate|rna|dna|methyl|platelet|clot|coagul|complement|nk |natural killer|t cell|b cell|monocyte|neutrophil|mast|sleep|pain|nerve|fibre|fiber|trial|placebo|randomi|treatment|drug|capillar|erythro|red blood|iron|vesicle|lipid|fatty|amino|kynuren|tryptophan|neuro|glia|inflamm|oxidat|mouse|mice|model|sex|estrogen|pituitar|vasopressin|thyroid|transcript|single-cell|machine learning",re.I)
cr=[p for p in L('c_recent.json') if not NONBIO.search(p['title'] or '') and BIO.search((p['title'] or '')+' '+(p.get('abstractText') or '')[:600]) and not re.search(r'letter|comment|editorial|erratum|correction',' '.join(p.get('pubType') or []),re.I) and not re.search(r'^(erratum|correction|retraction)',p['title'] or '',re.I)]
for p in cr: add(p,'c_recent_2024_26')
# (d) LC PEM/exercise/fatigue mechanism
MECH=re.compile(r"exercise|cpet|cardiopulmonary|exertion|muscle|mitochond|oxygen|lactate|orthostat|autonomic|pots\b|tachycardia|tilt|small.fib|cerebral|brain|mri|biomark|metabol|proteom|immun|cortisol|myalgic|ME/CFS|chronic fatigue syndrome|endotheli|microvasc|capillar|preload|hyperventil|breathing pattern|heart rate variab|chronotrop",re.I)
lc=[p for p in L('d_lc_pem.json') if not NONBIO.search(p['title'] or '') and not isrev(p)]
lc=[p for p in lc if (MECH.search(p['title'] or '') or re.search('post-exertional',p.get('abstractText') or '',re.I)) and not re.search(r'training|program|intervention|prevalence|incidence|cohort study of symptoms',p['title'] or '',re.I)]
for i,p in enumerate(lc[:250]): add(p,'d_lc_pem_mech',i+1)
lcrec=[p for p in lc[250:] if (p.get('pubYear') or '0')>='2024' and re.search('post-exertional|PEM|exercise|exertion',(p['title'] or '')+(p.get('abstractText') or ''),re.I)]
for p in lcrec: add(p,'d_lc_pem_mech_recent')
refs=list(out.values())
json.dump(refs,open(G+'reference_set.json','w'),indent=0)
from collections import Counter
c=Counter(); m=Counter()
for r in refs:
    for s in r['sets']:
        s=s.split('#')[0]; c[s]+=1; m[s]+= bool(r['match'])
for s in c: print(s,c[s],m[s])
print(len(refs), sum(bool(r['match']) for r in refs))
