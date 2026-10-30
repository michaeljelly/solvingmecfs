"""Harmonise every metabolite across cohorts.
Nodes = (study, feature). Edges when two features in different studies share an unambiguous identifier:
Metabolon COMP_ID, InChIKey first block (skeleton), HMDB, KEGG, PubChem CID, ChEBI, or normalised name.
An identifier value that points to >1 feature inside one study is treated as ambiguous for that study and not used.
Union-find gives harmonised metabolite groups ("hkey"); if a group holds >1 feature of one study, the most-detected one is used.
Hoel 2021 (Metabolon, summary statistics only) is added as an annotation-only study.
Output: harmonised_map.csv (hkey, label, study, feature, link_types), /tmp/work_mw/hmap.pkl
"""
import re, pickle
import numpy as np, pandas as pd
from collections import defaultdict

S = pickle.load(open('/tmp/work_mw/studies_full.pkl', 'rb'))
OUT = '/mnt/project-files/mecfs-research/deep2/metabolome_wide/'
RAW = '/mnt/project-files/mecfs-research/downloads/mwb/raw/'

# Hoel summary table (Metabolon HD4 + xenobiotics sheet)
h = pd.read_excel(RAW + 'Hoel2021/sdd1.xlsx', '3.FC & stats', header=1)
h = h.rename(columns={'MEall': 'FC', 'MEall.1': 'p'})
x = pd.read_excel(RAW + 'Hoel2021/sdd1.xlsx', '11.HD4 Xenobiotics', header=2).rename(columns={'Biochemical': 'Met', 'Sub_pathway': 'SUB_PATHWAY', 'MEall': 'FC', 'MEall.1': 'p'})
hx = pd.concat([h, x]).dropna(subset=['FC', 'p'])
hx['Met'] = hx.Met.astype(str).str.strip(); hx = hx.drop_duplicates('Met').reset_index(drop=True)
print('Hoel cols', hx.columns.tolist()[:12])
hann = pd.DataFrame({'comp_id': hx.get('COMP_ID'), 'kegg': hx.get('KEGG'), 'hmdb': hx.get('HMDB_ID'), 'sub': hx.SUB_PATHWAY}).fillna('').astype(str); hann.index = hx.Met.values
hann['comp_id'] = [str(int(float(v))) if str(v).replace('.0', '').isdigit() and str(v) not in ('0', '0.0') else '' for v in hann.comp_id]
for c in ['inchikey', 'pubchem', 'chebi', 'super']: hann[c] = ''
hx[['Met', 'SUB_PATHWAY', 'FC', 'p']].to_csv('/tmp/work_mw/hoel_summary.csv', index=False)


def split(v): return [t for t in re.split(r'[,;|\s]+', str(v)) if t and t.lower() not in ('nan', '0', 'none')]


def ids_of(r, name):
    out = []
    for t in split(r.get('comp_id', '')):
        t = t.replace('.0', '')
        if t.isdigit(): out.append('COMP:' + t)
    for t in split(r.get('inchikey', '')):
        if re.fullmatch(r'[A-Z]{14}(-[A-Z]{10}-[A-Z])?', t): out.append('IK:' + t[:14])
    for t in split(r.get('hmdb', '')):
        mm = re.fullmatch(r'HMDB0*(\d+)', t, re.I)
        if mm: out.append('HMDB:%07d' % int(mm.group(1)))
    for t in split(r.get('kegg', '')):
        if re.fullmatch(r'C\d{5}', t): out.append('KEGG:' + t)
    for t in split(r.get('pubchem', '')):
        t = t.replace('.0', '')
        if t.isdigit() and int(t) > 1: out.append('CID:' + t)
    for t in split(r.get('chebi', '')):
        mm = re.search(r'(\d+)', t)
        if 'CHEBI' in t.upper() and mm: out.append('CHEBI:' + mm.group(1))
    out.append('NAME:' + nname(name))
    lk = lipid_key(name)
    if lk: out.append('LIPID:' + lk)
    return out


def nname(s):
    s = str(s).lower().strip()
    s = re.sub(r'_pm_a\d+$', '', s).replace('_', ' ')
    s = s.replace('*', '').strip()
    s = re.sub(r'\s*\([^()]*\)\s*$', '', s)
    s = re.sub(r'\s+[ab]$', '', s)
    s = re.sub(r'^(l|d|dl|\(s\)|\(r\)|\(\+/-\))-', '', s)
    s = re.sub(r'\bacid\b', '', s).strip()
    s = re.sub(r'ic\b', 'ate', s)
    return re.sub(r'[^a-z0-9]', '', s)


CLS = {'gpc': 'PC', 'gpe': 'PE', 'gpi': 'PI', 'gps': 'PS', 'gpg': 'PG'}
def lipid_key(name):
    n = str(name).strip().replace('*', '')
    if re.search(r'\b(p-|o-|P-|O-)\d', n) or 'enyl' in n: return None
    m = re.search(r'-(GP[CEISG])\s*\(([\d:/n\-]+)\)', n)
    if m:
        ch = re.findall(r'(\d+):(\d+)', m.group(2))
        if not ch: return None
        c = sum(int(a) for a, b in ch); d = sum(int(b) for a, b in ch)
        return ('L' if len(ch) == 1 else '') + CLS[m.group(1).lower()] + f'{c}:{d}'
    m = re.search(r'sphingomyelin\s*\(d18:1/(\d+):(\d+)\)', n, re.I) or re.fullmatch(r'SM\s*\(d18:1/(\d+):(\d+)\)', n)
    if m: return f'SM{18 + int(m.group(1))}:{1 + int(m.group(2))}'
    m = re.fullmatch(r'SM\s*\(d(\d+):(\d+)\)(\s*[AB])?', n)
    if m: return f'SM{m.group(1)}:{m.group(2)}'
    m = re.search(r'-sphingosine\s*\(d18:1/(\d+):(\d+)\)', n) or re.fullmatch(r'Ceramide\s*\(d18:1/(\d+):(\d+)\)', n)
    if m: return f'Cer{18 + int(m.group(1))}:{1 + int(m.group(2))}'
    m = re.fullmatch(r'Ceramide\s*\(d(\d+):(\d+)\)(\s*[AB])?', n)
    if m: return f'Cer{m.group(1)}:{m.group(2)}'
    m = re.fullmatch(r'(L?P[CEISG])\s*\(?\s*(\d+):(\d+)\s*\)?\s*[AB]?;?', n)
    if m: return f'{m.group(1)}{m.group(2)}:{m.group(3)}'
    m = re.fullmatch(r'(L?P[CEISG])\s*\((\d+):(\d+)/(\d+):(\d+)\)', n)
    if m: return f'{m.group(1)}{int(m.group(2)) + int(m.group(4))}:{int(m.group(3)) + int(m.group(5))}'
    return None


studies = {k: v['ann'].loc[v['L'].columns] for k, v in S.items()}
studies['Hoel2021'] = hann
node_ids = {}
for st, a in studies.items():
    rows = {f: ids_of(r, f) for f, r in a.iterrows()}
    holders = defaultdict(list)
    for f, ii in rows.items():
        for i in set(ii): holders[i].append(f)
    for f, ii in rows.items():
        keep = []
        for i in set(ii):
            hs = holders[i]
            if len(hs) == 1: keep.append(i)
            elif i.startswith('LIPID:'): continue            # isomers A/B etc: ambiguous, drop
            elif len({nname(x) for x in hs}) == 1 or i.startswith('NAME:'):
                if f == hs[0]: keep.append(i)                # same compound measured twice (platforms): first listed
        node_ids[(st, f)] = keep

# union-find
par = {n: n for n in node_ids}
def fnd(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
by_id = defaultdict(list)
for n, ii in node_ids.items():
    for i in ii: by_id[i].append(n)
# avoid KEGG / PubChem generic links pulling in wrong compounds: only link via an id type if the two
# nodes do not contradict on a higher-priority id type (COMP, IK)
def conflict(a, b):
    for pre in ('COMP:', 'IK:'):
        A = {i for i in node_ids[a] if i.startswith(pre)}; B = {i for i in node_ids[b] if i.startswith(pre)}
        if A and B and not (A & B): return True
    return False
linktype = defaultdict(set)
for i, ns in by_id.items():
    if len(ns) < 2: continue
    for a in ns[1:]:
        b = ns[0]
        if a[0] == b[0] or conflict(a, b): continue
        ra, rb = fnd(a), fnd(b)
        if ra != rb: par[ra] = rb
        linktype[a].add(i.split(':')[0]); linktype[b].add(i.split(':')[0])

groups = defaultdict(list)
for n in node_ids: groups[fnd(n)].append(n)
PRI = ['BioMapAI', 'Germain2022', 'Hoel2021', 'Germain2018', 'Germain2020', 'Che2022', 'NagySzakal2018', 'Naviaux2016', 'NCNP2026', 'Armstrong2015']
rows = []
for gi, (root, ns) in enumerate(groups.items()):
    sts = sorted({s for s, _ in ns}, key=PRI.index)
    label = [f for s, f in ns if s == sts[0]][0]
    for st in sts:
        fs = [f for s, f in ns if s == st]
        if len(fs) > 1:
            det = studies[st].loc[fs, 'detect'] if 'detect' in studies[st] else pd.Series(0, index=fs)
            fs = [det.astype(float).fillna(0).idxmax()] if det.notna().any() else fs[:1]
        rows.append(dict(hkey=f'H{gi:05d}', label=label, study=st, feature=fs[0], n_group=len(ns), links=','.join(sorted(linktype[(st, fs[0])]))))
Mp = pd.DataFrame(rows)
# annotate with Metabolon super/sub pathway where available
sub = {}
for st in ['BioMapAI', 'Germain2022', 'Germain2018', 'Germain2020', 'Hoel2021']:
    a = studies[st]
    for _, r in Mp[Mp.study == st].iterrows():
        sub.setdefault(r.hkey, (a.loc[r.feature].get('super', ''), a.loc[r.feature].get('sub', '')))
Mp['super'] = Mp.hkey.map(lambda k: sub.get(k, ('', ''))[0]); Mp['sub'] = Mp.hkey.map(lambda k: sub.get(k, ('', ''))[1])
Mp['k_studies'] = Mp.groupby('hkey').study.transform('nunique')
Mp.to_csv(OUT + 'harmonised_map.csv', index=False)
pickle.dump(Mp, open('/tmp/work_mw/hmap.pkl', 'wb'))
k = Mp.drop_duplicates('hkey').k_studies
print('groups', Mp.hkey.nunique(), 'in >=2 studies', (k >= 2).sum(), '>=3', (k >= 3).sum(), '>=5', (k >= 5).sum())
print(pd.crosstab(Mp.study, Mp.k_studies >= 2))
# how many non-Metabolon features map to a Metabolon group
met = set(Mp[Mp.study.isin(['BioMapAI', 'Germain2022', 'Germain2018', 'Germain2020', 'Hoel2021'])].hkey)
for st in ['Che2022', 'NagySzakal2018', 'Naviaux2016', 'NCNP2026', 'Armstrong2015']:
    sub_ = Mp[Mp.study == st]; print(st, len(sub_), 'mapped to Metabolon:', sub_.hkey.isin(met).sum())
print('largest groups', Mp.groupby('hkey').size().sort_values().tail(3))
