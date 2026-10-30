"""Download Risteys endpoint pages and parse key figures, CodeWAS and LabWAS tables to CSV.
Also pull legacy R10 survival associations (before/after with HR) where available."""
import requests, sys, os, re, json, csv, time
from bs4 import BeautifulSoup
OUT = "/mnt/project-files/mecfs-research/deep2/epi_literature/risteys"; os.makedirs(OUT, exist_ok=True)
eps = sys.argv[1:]
for ep in eps:
    fn = f"/tmp/work_epi/ris/{ep}.page.html"
    if not os.path.exists(fn):
        open(fn, "w").write(requests.get(f"https://risteys.finngen.fi/endpoints/{ep}", timeout=180).text)
    soup = BeautifulSoup(open(fn).read(), "lxml")
    for sec in ["codewas", "labwas"]:
        box = soup.find(id=sec)
        # section header anchor id; the table is the next role=table / rowgroup after it
        rows = []
        if box is not None:
            art = box.find_parent("article") or box.parent
            want = 7 if sec == "codewas" else 13
            for r in art.select('div[role="rowgroup"] > div[role="row"]'):
                cells = [re.sub(r"\s+", " ", c.get_text(" ", strip=True)) for c in r.find_all("div", recursive=False)]
                if len(cells) == want and "▲" not in "".join(cells) and cells[0] not in ("Code", "(ID) Name"): rows.append(cells)
        with open(f"{OUT}/{ep}_{sec}.csv", "w", newline="") as f:
            w = csv.writer(f)
            if sec == "codewas": w.writerow(["code", "vocab", "description", "OR", "mlog10p", "n_cases", "n_controls"])
            else: w.writerow(["concept", "n_cases", "n_controls", "OR_measured", "mlog10p_measured", "meanN_cases", "meanN_controls", "mean_cases", "mean_controls", "unit", "mlog10p_value", "n_cases_v", "n_controls_v"])
            w.writerows(rows)
        print(ep, sec, len(rows))
    txt = re.sub(r"\s+", " ", soup.get_text(" "))
    m = re.search(r"-FinRegistry- Key figures All Female Male Number of individuals Whole population (\S+) (\S+) (\S+)", txt)
    m2 = re.search(r"-FinnGen- Key figures All Female Male Number of individuals (\S+) (\S+) (\S+)", txt)
    m3 = re.search(r"Median age at first event \(years\) Whole population (\S+) (\S+) (\S+)", txt)
    json.dump({"FR_n": m.groups() if m else None, "FG_n": m2.groups() if m2 else None, "FR_median_age": m3.groups() if m3 else None}, open(f"{OUT}/{ep}_keyfigures.json", "w"))
    for rel in ["r10", "r9"]:
        try:
            j = requests.get(f"https://{rel}.risteys.finngen.fi/api/endpoints/{ep}/assocs.json", timeout=120).json()
            json.dump(j, open(f"{OUT}/{ep}_assocs_{rel}.json", "w"))
        except Exception as e: print("assoc fail", ep, rel, e)
    time.sleep(1)
