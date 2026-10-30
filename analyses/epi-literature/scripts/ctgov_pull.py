"""Pull all interventional studies for ME/CFS, long COVID, POTS from ClinicalTrials.gov API v2."""
import requests, json, time, sys
OUT = "/tmp/work_epi/ctgov"
import os; os.makedirs(OUT, exist_ok=True)
queries = {
  "mecfs": "chronic fatigue syndrome OR myalgic encephalomyelitis OR ME/CFS OR systemic exertion intolerance disease",
  "longcovid": "long COVID OR post-acute COVID-19 syndrome OR post COVID-19 condition OR PASC OR post-COVID syndrome OR post-acute sequelae of SARS-CoV-2",
  "pots": "postural orthostatic tachycardia syndrome OR orthostatic intolerance",
}
for name, q in queries.items():
    studies = []; tok = None
    while True:
        params = {"query.cond": q, "query.term": "AREA[StudyType]INTERVENTIONAL", "pageSize": 200, "format": "json"}
        if tok: params["pageToken"] = tok
        r = requests.get("https://clinicaltrials.gov/api/v2/studies", params=params, timeout=120)
        r.raise_for_status(); j = r.json()
        studies += j.get("studies", []); tok = j.get("nextPageToken")
        print(name, len(studies), flush=True)
        if not tok: break
        time.sleep(0.5)
    json.dump(studies, open(f"{OUT}/{name}.json", "w"))
