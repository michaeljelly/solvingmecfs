"""Generic PubMed fetch: python3 pubmed_generic.py NAME 'query' -> /tmp/work_epi/pm_NAME.json"""
import requests, time, sys, json
import xml.etree.ElementTree as ET
name, term = sys.argv[1], sys.argv[2]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
x = ET.fromstring(requests.get(E + "esearch.fcgi", params={"db": "pubmed", "term": term, "usehistory": "y", "retmax": 0}, timeout=120).text)
n = int(x.find("Count").text); we = x.find("WebEnv").text; qk = x.find("QueryKey").text
print(name, "count", n, flush=True)
recs = []
for start in range(0, n, 400):
    for attempt in range(5):
        try:
            rr = requests.get(E + "efetch.fcgi", params={"db": "pubmed", "WebEnv": we, "query_key": qk, "retstart": start, "retmax": 400, "retmode": "xml"}, timeout=300)
            root = ET.fromstring(rr.text); break
        except Exception as e:
            print("retry", e); time.sleep(3)
    for a in root.findall(".//PubmedArticle"):
        t = a.find(".//ArticleTitle"); aus = a.findall(".//AuthorList/Author")
        nm = lambda x: (x.findtext("LastName") or x.findtext("CollectiveName") or "") if x is not None else ""
        recs.append(dict(pmid=a.findtext(".//PMID"), year=a.findtext(".//PubDate/Year") or (a.findtext(".//PubDate/MedlineDate") or "")[:4],
            first=nm(aus[0]) if aus else "", last=nm(aus[-1]) if aus else "", journal=a.findtext(".//Journal/ISOAbbreviation") or "",
            title="".join(t.itertext()) if t is not None else "", abstract=" ".join("".join(z.itertext()) for z in a.findall(".//Abstract/AbstractText")),
            ptypes=[p.text for p in a.findall(".//PublicationType")], mesh=[m.findtext("DescriptorName") for m in a.findall(".//MeshHeading")]))
    time.sleep(0.4)
json.dump(recs, open(f"/tmp/work_epi/pm_{name}.json", "w")); print(name, "saved", len(recs))
