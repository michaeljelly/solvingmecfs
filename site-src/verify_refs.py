#!/usr/bin/env python3
"""Check references marked unverified in data/treatments.json and the atlas research JSON against Crossref, PubMed and
ClinicalTrials.gov. A DOI or PubMed reference is marked verified when the record exists and its first author's surname
and year match the citation text; a PubMed search URL that returns exactly one matching paper is replaced by that PMID;
a trial registration is verified when the NCT record exists. Everything else is left as it was and reported.

usage: python3 -I verify_refs.py <site_src_dir> <atlas_research_dir> [--write]
"""
import glob, json, os, re, sys, time, unicodedata, urllib.parse, urllib.request
SRC, RES = sys.argv[1], sys.argv[2]
WRITE = "--write" in sys.argv
UA = {"User-Agent": "solvingmecfs-refcheck (mailto:hello@solvingmecfs.org)"}

def get(url):
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=20) as r:
                return r.status, r.read().decode("utf8", "replace")
        except urllib.error.HTTPError as e:
            if e.code == 404: return 404, ""
        except Exception:
            pass
        time.sleep(1.5 * (i + 1))
    return 0, ""

def norm(s):
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()

def cite_bits(cite):
    m = re.match(r"\s*((?:(?:van|de|der|den|von|da|di|le|la)\s+)*[A-Za-zÀ-ÿ'\-]+)", cite, re.I)
    y = re.search(r"\b(19|20)\d\d\b", cite)
    return (norm(m.group(1)) if m else ""), (y.group(0) if y else "")

def check_doi(doi, cite):
    st, body = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))
    if st != 200: return None, f"crossref {st}"
    w = json.loads(body)["message"]
    names = [norm(a.get("family", "")) for a in w.get("author", [])]
    years = {str(p["date-parts"][0][0]) for k in ("published-print", "published-online", "issued", "created") if (p := w.get(k)) and p.get("date-parts") and p["date-parts"][0][0]}
    au, yr = cite_bits(cite)
    ok_au = bool(au) and (au in names[:1] or any(au == n for n in names))
    ok_yr = (not yr) or yr in years or any(abs(int(yr) - int(x)) <= 1 for x in years)
    return (ok_au and ok_yr), f"crossref: {names[:2]} {sorted(years)}"

def check_pmid(pmid, cite):
    st, body = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json")
    if st != 200: return None, f"pubmed {st}"
    r = json.loads(body)["result"].get(str(pmid), {})
    names = [norm(a["name"].split()[0]) for a in r.get("authors", [])]
    au, yr = cite_bits(cite)
    yrs = r.get("pubdate", "")[:4]
    ok = bool(au) and au in names and ((not yr) or abs(int(yr) - int(yrs or 0)) <= 1)
    return ok, f"pubmed: {names[:2]} {yrs}"

def check(ref):
    url, cite = (ref.get("url") or "").strip(), ref.get("cite", "")
    if m := re.search(r"doi\.org/(10\.[^\s?#]+)", url):
        return check_doi(m.group(1), cite) + (None,)
    if m := re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", url):
        return check_pmid(m.group(1), cite) + (None,)
    if m := re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/\?term=([^&#]+)", url):
        st, body = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term=" + m.group(1))
        ids = json.loads(body)["esearchresult"]["idlist"] if st == 200 else []
        if len(ids) == 1:
            ok, note = check_pmid(ids[0], cite)
            return ok, note, (f"https://pubmed.ncbi.nlm.nih.gov/{ids[0]}/" if ok else None)
        return None, f"pubmed search: {len(ids)} hits", None
    if m := re.search(r"(NCT\d{8})", url + " " + cite):
        st, _ = get("https://clinicaltrials.gov/api/v2/studies/" + m.group(1))
        return (st == 200), f"ctgov {st}", None
    if "nice.org.uk/guidance/ng206" in url:
        st, body = get(url)
        return (st == 200 and "NG206" in body), f"nice {st}", None
    return None, "not checkable automatically", None

refs = []
T = json.load(open(os.path.join(SRC, "data", "treatments.json"), encoding="utf8"))
for t in T["treatments"]:
    refs += [(r, "treatment " + t["id"]) for r in t.get("refs") or [] if not r.get("verified")]
docs = {}
for p in glob.glob(os.path.join(RES, "*.json")):
    d = docs[p] = json.load(open(p, encoding="utf8"))
    for f in d.get("findings", []):
        for k in ("studies", "contradicting"):
            refs += [(r, "finding " + f["id"]) for r in f.get(k) or [] if not r.get("verified")]
    for h in d.get("hypotheses", []):
        refs += [(r, "hypothesis " + h["id"]) for r in h.get("key_refs") or [] if not r.get("verified")]

n_ok = 0
for r, where in refs:
    ok, note, newurl = check(r)
    flag = "OK " if ok else ("BAD" if ok is False else "-- ")
    print(f"{flag} {where} | {r.get('cite')} | {note}")
    if ok:
        n_ok += 1
        r["verified"] = True
        if newurl: r["url"] = newurl
print(f"{n_ok} of {len(refs)} unverified references now verified")
if WRITE:
    json.dump(T, open(os.path.join(SRC, "data", "treatments.json"), "w", encoding="utf8"), ensure_ascii=False, indent=1)
    for p, d in docs.items():
        json.dump(d, open(p, "w", encoding="utf8"), ensure_ascii=False, indent=1)
