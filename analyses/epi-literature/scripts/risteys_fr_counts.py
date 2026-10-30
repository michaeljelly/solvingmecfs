"""Get FinRegistry whole-population case counts for endpoints (from Risteys pages) and compute O/E overlap with G6_POSTVIRFAT in FinRegistry."""
import requests, re, sys, json, os, time, pandas as pd, numpy as np
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests
eps = sys.argv[1:]
out = []
for ep in eps:
    fn = f"/tmp/work_epi/ris/{ep}.page.html"
    if not os.path.exists(fn):
        r = requests.get(f"https://risteys.finngen.fi/endpoints/{ep}", timeout=180); open(fn, "w").write(r.text); time.sleep(0.5)
    t = re.sub(r"<[^>]+>", " ", open(fn).read()); t = re.sub(r"\s+", " ", t)
    m = re.search(r"-FinRegistry- Key figures All Female Male Number of individuals Whole population (\d+)", t)
    out.append((ep, int(m.group(1)) if m else np.nan))
c = pd.DataFrame(out, columns=["endpoint", "n_FR"])
rel = pd.read_csv("/mnt/project-files/mecfs-research/deep2/epi_literature/risteys/POSTVIRFAT_overlap_enrichment.csv")
x = c.merge(rel[["endpoint", "longname", "o_FR", "o_FG", "num_cases", "OE_FG", "q_FG"]], on="endpoint", how="left")
N = 7166416; A = 1622
x["E_FR"] = A * x.n_FR / N; x["OE_FR"] = x.o_FR / x.E_FR
x["p_FR"] = [poisson.sf(o - 1, e) if o == o and e == e else np.nan for o, e in zip(x.o_FR, x.E_FR)]
x.to_csv("/mnt/project-files/mecfs-research/deep2/epi_literature/risteys/POSTVIRFAT_selected_FR_vs_FG_enrichment.csv", index=False)
pd.set_option("display.width", 250)
print(x[["endpoint", "longname", "o_FG", "OE_FG", "o_FR", "n_FR", "E_FR", "OE_FR", "p_FR"]].round(3).to_string(index=False))
