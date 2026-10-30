"""Mine ME/CFS abstracts for biomarker findings (patients vs controls): concept, direction, sample size.
Builds a replication catalogue across independent groups (distinct last authors)."""
import json, re, sys, collections
import pandas as pd, numpy as np
sys.path.insert(0, "/mnt/project-files/mecfs-research/deep2/epi_literature/scripts")
from biomarker_dict import D
OUT = "/mnt/project-files/mecfs-research/deep2/epi_literature/literature"
R = json.load(open("/tmp/work_epi/pubmed_mecfs.json"))
EXCL = {"Review", "Systematic Review", "Meta-Analysis", "Editorial", "Comment", "Letter", "Published Erratum", "News", "Case Reports", "Practice Guideline", "Guideline", "Retracted Publication", "Retraction of Publication", "Congress"}
papers = [r for r in R if r["abstract"] and not (set(r["ptypes"]) & EXCL) and not ("Animals" in r["mesh"] and "Humans" not in r["mesh"])]
print("original research with abstract:", len(papers))
UP = r"\b(higher|increased|increase|elevated|elevation|greater|raised|upregulated|up-regulated|enhanced|excess\w*|more (frequent|common|prevalent)|overexpress\w*|augmented|exaggerated|prolonged|above (normal|expected|reference)|abnormally high|excessive)\b"
DOWN = r"\b(lower|decreased|decrease|reduced|reduction|diminished|impaired|deficien\w+|deficit\w*|less|downregulated|down-regulated|depleted|depletion|attenuated|fewer|loss|low|blunted|declin\w+|suppressed|shorter|smaller|hypo\w+|subnormal|below (normal|expected|reference|predicted)\w*|abnormally low)\b"
NULL = r"(no (significant |statistically significant |apparent |major )?(difference|differences|change|changes|association|alteration|abnormalit\w+|evidence)|did not differ|not differ|were not (significantly )?(different|altered|elevated|reduced|increased|decreased|changed)|was not (significantly )?(different|altered|elevated|reduced|increased|decreased|changed)|similar (in|between|to|across)|comparable|not significantly|no group difference|within (the )?normal (range|limits)|unchanged|not detected|non-?significant)"
PAT = r"(patients?|ME/CFS|\bCFS\b|\bME\b|cases|sufferers|subjects with|individuals with|participants with|SEID|people with)"
CTRL = r"(controls?\b|healthy|\bHCs?\b|\bHV\b|compared|versus|\bvs\.?|than\b|sedentary)"
INV = r"(controls?|HCs?|healthy (subjects|volunteers|individuals|participants|people))\s+(had|showed|displayed|exhibited|demonstrated|presented)\s"
SKIP = r"\b(aims?|objectives?|purpose|hypothes\w+|whether|to determine|to investigate|to examine|to compare|to assess|to evaluate|we (examined|investigated|measured|assessed|evaluated|studied|explored|sought)|previous(ly)?|prior (studies|reports|work)|has been (reported|suggested|proposed|implicated)|have been (reported|suggested|proposed|implicated)|literature|might|remains? (unclear|unknown|controversial)|is thought|are thought)\b"
skip_re = re.compile(SKIP, re.I)
within_re = re.compile(r"(than|compared (to|with)|versus|vs\.?) (those|patients|subjects|participants|individuals|cases)? ?(without|with(?! ME| CFS| chronic))|patients with (high|low|severe|mild|moderate)", re.I)
other_re = re.compile(r"fibromyalgia|depress\w+|multiple sclerosis|\bMS\b|\bPBC\b|rheumatoid|cancer|Q[- ]fever|\bQFS\b|lupus|Gulf War|\bGWI\b|endometriosis|schizophren\w+|rats?\b|mice\b|mouse", re.I)
mecfs_re = re.compile(r"ME/CFS|\bCFS\b|\bME\b|chronic fatigue|myalgic|SEID", re.I)
cre = {k: re.compile(v[1], re.I) for k, v in D.items()}
up_re, dn_re, null_re, pat_re, ctrl_re, inv_re = [re.compile(x, re.I) for x in (UP, DOWN, NULL, PAT, CTRL, INV)]
def nsize(ab):
    # patients: number before patient words; controls: number before control words
    nums = lambda pat: [int(x.replace(",", "")) for x in re.findall(pat, ab, re.I)]
    p = nums(r"(\d[\d,]*)\s+(?:\w+\s){0,4}?(?:patients|participants|subjects|individuals|cases|women|adults|adolescents|people|persons)\b(?!\s+(?:without|free))")
    p += nums(r"(?:ME/CFS|CFS|ME|patients)\s*\(\s*n\s*=\s*(\d[\d,]*)")
    c = nums(r"(\d[\d,]*)\s+(?:\w+\s){0,3}?(?:healthy controls|controls|healthy (?:subjects|volunteers|individuals|participants))\b")
    c += nums(r"(?:controls?|HCs?)\s*\(\s*n\s*=\s*(\d[\d,]*)")
    p = [x for x in p if 3 <= x <= 200000 and not (1900 <= x <= 2030)]
    c = [x for x in c if 3 <= x <= 200000 and not (1900 <= x <= 2030)]
    return (max(p) if p else np.nan, max(c) if c else np.nan)
rows = []
for r in papers:
    ab = r["abstract"]; text = r["title"] + ". " + ab
    npat, nctl = nsize(ab)
    has_ctrl_design = bool(re.search(r"control|healthy|HCs?\b", ab, re.I))
    for sent in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", text):
        if not pat_re.search(sent) or not ctrl_re.search(sent): continue
        if skip_re.search(sent): continue
        if within_re.search(sent): continue
        if other_re.search(sent) and not mecfs_re.search(sent): continue
        # split into clauses on ; and contrastive conjunctions
        for cl in re.split(r";|,? whereas |,? while |,? but |,? although ", sent):
            hits = [(k, m.start()) for k in cre for m in [cre[k].search(cl)] if m]
            if not hits: continue
            if not (up_re.search(cl) or dn_re.search(cl) or null_re.search(cl)): continue
            inv = bool(inv_re.search(cl))
            for k, pos in hits:
                if null_re.search(cl):
                    d = "null"
                else:
                    ups = [m.start() for m in up_re.finditer(cl)]; dns = [m.start() for m in dn_re.finditer(cl)]
                    du = min([abs(u - pos) for u in ups], default=9999); dd = min([abs(x - pos) for x in dns], default=9999)
                    if min(du, dd) > 120: continue
                    d = "up" if du < dd else "down"
                    if inv: d = {"up": "down", "down": "up"}[d]
                rows.append(dict(pmid=r["pmid"], year=r["year"], first=r["first"], last=r["last"], journal=r["journal"], concept=k, category=D[k][0],
                                 direction=d, n_pat=npat, n_ctl=nctl, ctrl_design=has_ctrl_design, clause=cl.strip()[:400], title=r["title"][:150]))
X = pd.DataFrame(rows)
# one call per paper x concept: majority direction (non-null beats null if any directional claim)
def agg(g):
    c = collections.Counter(g.direction)
    dirs = {k: v for k, v in c.items() if k != "null"}
    if dirs:
        d = max(dirs, key=dirs.get)
        if len(dirs) > 1 and list(dirs.values()).count(max(dirs.values())) > 1: d = "mixed"
    else: d = "null"
    return pd.Series(dict(direction=d, n_up=c.get("up", 0), n_down=c.get("down", 0), n_null=c.get("null", 0), year=g.year.iloc[0], first=g["first"].iloc[0], last=g["last"].iloc[0],
                          n_pat=g.n_pat.iloc[0], n_ctl=g.n_ctl.iloc[0], ctrl_design=g.ctrl_design.iloc[0], category=g.category.iloc[0], example=g.clause.iloc[0], title=g.title.iloc[0]))
P = X.groupby(["concept", "pmid"]).apply(agg).reset_index()
P = P[P.ctrl_design]
X.to_csv(f"{OUT}/biomarker_mentions_all.csv", index=False)
P.to_csv(f"{OUT}/biomarker_paper_calls.csv", index=False)
# catalogue: independent groups = distinct last authors
cat = []
for k, g in P.groupby("concept"):
    grp = lambda d: g[g.direction == d]["last"].str.lower().nunique()
    pm = lambda d: (g.direction == d).sum()
    up, dn, nu, mx = grp("up"), grp("down"), grp("null"), grp("mixed")
    tot_n = g.n_pat.fillna(0).sum()
    status = ("replicated up" if up >= 3 and up >= 2 * dn else "replicated down" if dn >= 3 and dn >= 2 * up else
              "contradictory" if up >= 2 and dn >= 2 else "single/limited" if up + dn <= 2 else "mixed")
    cat.append(dict(concept=k, category=D[k][0], papers=len(g), groups_up=up, groups_down=dn, groups_null=nu, groups_mixed=mx,
                    papers_up=pm("up"), papers_down=pm("down"), papers_null=pm("null"), sum_patients=int(tot_n), status=status,
                    first_year=int(pd.to_numeric(g.year, errors="coerce").min()) if len(g) else None,
                    last_year=int(pd.to_numeric(g.year, errors="coerce").max()) if len(g) else None))
C = pd.DataFrame(cat).sort_values(["category", "papers"], ascending=[True, False])
C.to_csv(f"{OUT}/biomarker_replication_catalogue.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
print(C.to_string(index=False))
