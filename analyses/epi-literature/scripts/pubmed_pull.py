"""Fetch all PubMed records for ME/CFS (title/abstract/MeSH) as XML batches, then parse to JSONL."""
import requests, time, os, re, json
import xml.etree.ElementTree as ET
OUT = "/tmp/work_epi/pubmed"; os.makedirs(OUT, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
term = ('"fatigue syndrome, chronic"[MeSH Terms] OR "chronic fatigue syndrome"[tiab] OR "myalgic encephalomyelitis"[tiab] '
        'OR "myalgic encephalopathy"[tiab] OR "ME/CFS"[tiab] OR "systemic exertion intolerance"[tiab] OR "post-viral fatigue"[tiab] '
        'OR "postviral fatigue"[tiab] OR "post-infectious fatigue"[tiab] OR "postinfectious fatigue"[tiab]')
r = requests.get(E + "esearch.fcgi", params={"db": "pubmed", "term": term, "usehistory": "y", "retmax": 0}, timeout=120)
x = ET.fromstring(r.text); n = int(x.find("Count").text); we = x.find("WebEnv").text; qk = x.find("QueryKey").text
print("count", n, flush=True)
for start in range(0, n, 500):
    fn = f"{OUT}/batch_{start:06d}.xml"
    if os.path.exists(fn) and os.path.getsize(fn) > 1000: continue
    for attempt in range(5):
        try:
            rr = requests.get(E + "efetch.fcgi", params={"db": "pubmed", "WebEnv": we, "query_key": qk, "retstart": start, "retmax": 500, "retmode": "xml"}, timeout=300)
            if rr.status_code == 200 and "<PubmedArticleSet" in rr.text:
                open(fn, "w").write(rr.text); break
        except Exception as e:
            print("err", e)
        time.sleep(3)
    print(start, flush=True); time.sleep(0.4)
# parse
recs = []
for fn in sorted(os.listdir(OUT)):
    if not fn.endswith(".xml"): continue
    root = ET.parse(f"{OUT}/{fn}").getroot()
    for a in root.findall(".//PubmedArticle"):
        pmid = a.findtext(".//PMID")
        title = "".join(a.find(".//ArticleTitle").itertext()) if a.find(".//ArticleTitle") is not None else ""
        ab = " ".join("".join(t.itertext()) for t in a.findall(".//Abstract/AbstractText"))
        year = a.findtext(".//PubDate/Year") or (a.findtext(".//PubDate/MedlineDate") or "")[:4]
        au = a.find(".//AuthorList/Author")
        first = (au.findtext("LastName") or au.findtext("CollectiveName") or "") if au is not None else ""
        journal = a.findtext(".//Journal/ISOAbbreviation") or ""
        ptypes = [p.text for p in a.findall(".//PublicationType")]
        mesh = [m.findtext("DescriptorName") for m in a.findall(".//MeshHeading")]
        recs.append(dict(pmid=pmid, year=year, first=first, journal=journal, title=title, abstract=ab, ptypes=ptypes, mesh=mesh))
with open("/tmp/work_epi/pubmed_mecfs.jsonl", "w") as f:
    for r_ in recs: f.write(json.dumps(r_) + "\n")
print("parsed", len(recs))
