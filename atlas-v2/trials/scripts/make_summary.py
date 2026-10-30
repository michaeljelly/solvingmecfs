"""Write summary.md and coverage_check.json from trials.json. Usage: python3 -I make_summary.py <base>"""
import json, sys, collections, re

base = sys.argv[1].rstrip("/")
D = json.load(open(base + "/trials.json"))
T, M = D["trials"], D["meta"]
C = collections.Counter
I = {r["id"]: r for r in T}


def table(counter, head):
    rows = ["| %s | Trials |" % head, "|---|---:|"]
    rows += ["| %s | %d |" % (k, v) for k, v in counter.most_common()]
    return "\n".join(rows)


def rate(rs):
    c = C(r["result"] for r in rs)
    n = len(rs)
    return n, (round(100 * c["positive"] / n) if n else 0), (round(100 * c["null"] / n) if n else 0)


out = []
w = out.append
w("# ME/CFS and long COVID treatment trials: database summary\n")
w("Built %s. Data file: `trials.json` (%d trials). Scripts are in `scripts/`; raw downloads are in `raw/`.\n" % (M["generated"], len(T)))
w("Sources: ClinicalTrials.gov API v2 (all interventional ME/CFS records, plus long COVID records screened by target symptom), "
  "Europe PMC (published trials, matched to registry records), ISRCTN and EU CTIS. "
  "`verified: true` means the result was read in a published abstract or in ClinicalTrials.gov posted results. "
  "No result was taken from memory or from press releases.\n")

w("## Counts\n")
w("### By population\n" + table(C(r["population"] for r in T), "Population") + "\n")
w("### By design\n" + table(C(r["design"] for r in T), "Design") + "\n")
w("\"RCT-unspecified\" covers ISRCTN or CTIS records that are randomised but do not name the comparator.\n")
st = C()
for r in T:
    s = r["status"]
    if s.startswith("CTIS"):
        s = "EU CTIS (status code not mapped)"
    elif r["registry"] == "ISRCTN":
        s = "ISRCTN (status not given)"
    elif s == "published":
        s = "published only (not on CT.gov, ISRCTN or CTIS)"
    st[s] += 1
w("### By status\n" + table(st, "Status") + "\n")
w("### By result\n" + table(C(r["result"] for r in T), "Result") + "\n")
w("Two values were added to the requested set. **feasibility-only** means the trial reported only recruitment, adherence or safety. "
  "**not-assessable** means a paper or posted result exists but efficacy cannot be read from it, for example because there is no abstract or the trial was terminated with unusable data. "
  "**withdrawn** marks trials withdrawn before enrolment.\n")
w("### By intervention class\n" + table(C(r["intervention_class"] for r in T), "Class") + "\n")
w("### By model step (one trial can count under several steps)\n")
names = M["model_steps_key"]
sc = C(s for r in T for s in r["model_steps"])
w("| Step | Trials |\n|---|---:|")
for s in range(1, 9):
    w("| %d %s | %d |" % (s, names[str(s)], sc[s]))
w("| no step assigned (for example exercise, vitamins, acupuncture) | %d |\n" % sum(1 for r in T if not r["model_steps"]))
w("### Verified results by population\n")
vt = C((r["population"], r["result"]) for r in T if r["verified"])
w("| Population | positive | null | mixed | negative/harm | feasibility-only |\n|---|---:|---:|---:|---:|---:|")
for p in ("ME/CFS", "long COVID", "both"):
    w("| %s | %d | %d | %d | %d | %d |" % (p, vt[(p, "positive")], vt[(p, "null")], vt[(p, "mixed")], vt[(p, "negative/harm")], vt[(p, "feasibility-only")]))
w("")

# design vs positivity
v = [r for r in T if r["verified"] and r["result"] in ("positive", "null", "mixed", "negative/harm")]
w("### How often trials are positive, by design (verified efficacy results only)\n")
w("| Group | Trials | % positive | % null |\n|---|---:|---:|---:|")
for lab, rs in [
    ("All verified efficacy results", v),
    ("Uncontrolled (single-arm, open-label, n-of-1)", [r for r in v if r["design"] in ("single-arm", "open-label", "n-of-1")]),
    ("Randomised against usual care or wait list", [r for r in v if r["design"] == "RCT-usual-care"]),
    ("Blinded placebo or sham (parallel or crossover)", [r for r in v if r["design"] in ("RCT-placebo", "crossover") and r["blinded"]]),
    ("Blinded placebo or sham with n of 100 or more", [r for r in v if r["design"] in ("RCT-placebo", "crossover") and r["blinded"] and (r["n_actual"] or 0) >= 100]),
    ("TCM, herbal, acupuncture or moxibustion", [r for r in v if re.search(r"TCM|herbal|acupunct|moxib|qigong|tuina|chinese|granule|decoction", r["intervention"] + r["title"], re.I)]),
]:
    n, p, nl = rate(rs)
    w("| %s | %d | %d | %d |" % (lab, n, p, nl))
w("\nThe more tightly a trial is controlled, the less often it is positive. This is the strongest pattern in the database.\n")

w("## Exclusions\n")
w("ClinicalTrials.gov records screened out (ME/CFS and long COVID searches combined, after deduplication):\n")
def exl(k):
    k2 = k.lower()
    for pat, lab in [("dropped", "dropped on manual review"), ("treatment trial|diagnostic", "not a treatment trial (mechanistic, diagnostic, observational)"),
                     ("population", "not an ME/CFS or long COVID population"), ("lung|breath", "primary outcome lung/breathlessness only"),
                     ("psych|sleep", "primary outcome psychiatric or sleep only"), ("smell", "primary outcome smell/taste only"),
                     ("pain", "primary outcome pain only"), ("cardiac|organ|hair", "primary outcome other organ (heart, skin, GI, eye, etc.)"),
                     ("acute", "acute COVID, prevention or post-discharge recovery"), ("vaccine", "vaccine trial")]:
        if re.search(pat, k2):
            return lab
    return k
ex = C()
for k, n in M["exclusions_ctgov"].items():
    ex[exl(k)] += n
w(table(ex, "Reason") + "\n")
w("The long COVID screen kept trials aimed at fatigue, PEM, exercise tolerance, cognition, autonomic or orthostatic symptoms, or global function. "
  "It excluded trials of lung function or breathlessness only, smell or taste only, vaccines, acute COVID, and mental health only. "
  "Further exclusions after manual review: 3 diagnostic or mechanistic studies, 1 CT.gov record, and 56 Europe PMC publications that were not treatment trials or were out of scope (52 on abstract review plus 4 earlier). Also 6 duplicate publications or registry entries were merged into existing records, and the 2 platform parent records (RECOVER-AUTONOMIC, RECOVER-ENERGIZE) were split into their arms.\n")
for reg, c in M["registry_screening"].items():
    w("- %s: %s" % (reg, "; ".join("%s %d" % (k, n) for k, n in c.items())))
w("")

w("## Overdue: completed more than 12 months ago with no results\n")
od = [r for r in T if r["result"] == "unreported-after-completion"]
w("%d trials have a primary completion date before 2025-10-29, have no results posted on ClinicalTrials.gov, and have no results paper that we could find in Europe PMC. "
  "They are %d ME/CFS, %d long COVID and %d both. "
  "The %d with UNKNOWN status are registry records that have not been updated past their completion date. Most of these probably finished, but some may never have run.\n"
  % (len(od), sum(r["population"] == "ME/CFS" for r in od), sum(r["population"] == "long COVID" for r in od), sum(r["population"] == "both" for r in od), sum(r["status"] == "UNKNOWN" for r in od)))
w("### Blinded or randomised trials, n of 50 or more (actual or planned), completed or terminated\n")
w("| Trial | Intervention | Population | n | Primary completion | Sponsor |\n|---|---|---|---:|---|---|")
big = [r for r in od if r["status"] in ("COMPLETED", "TERMINATED") and r["design"].startswith("RCT") and (r["n_actual"] or r["n_planned"] or 0) >= 50]
for r in sorted(big, key=lambda r: -(r["n_actual"] or r["n_planned"] or 0)):
    w("| %s | %s | %s | %s | %s | %s |" % (r["id"], r["intervention"][:50], r["population"], r["n_actual"] or r["n_planned"], r["primary_completion"], (r["sponsor"] or "")[:40]))
w("\n### All other overdue trials\n")
rest = [r for r in od if r not in big]
w("| Trial | Intervention | Status | n | Primary completion |\n|---|---|---|---:|---|")
for r in sorted(rest, key=lambda r: (r["status"], r["primary_completion"] or "")):
    w("| %s | %s | %s | %s | %s |" % (r["id"], r["intervention"][:45], r["status"], r["n_actual"] or r["n_planned"] or "", r["primary_completion"]))
w("")

w("## Largest and best-designed trials\n")
w("Blinded placebo- or sham-controlled trials with a verified result, sorted by size (n of 100 or more):\n")
w("| Trial | Population | Intervention | n | Result | Finding |\n|---|---|---|---:|---|---|")
best = [r for r in v if r["design"] in ("RCT-placebo", "crossover") and r["blinded"] and (r["n_actual"] or 0) >= 100]
for r in sorted(best, key=lambda r: -(r["n_actual"] or 0)):
    w("| %s | %s | %s | %d | %s | %s |" % (r["id"], r["population"], r["intervention"][:40], r["n_actual"], r["result"], r["result_summary"].replace("|", "/")[:170]))
w("\nLarge unblinded RCTs (n of 200 or more) with a verified result:\n")
w("| Trial | Population | Intervention | n | Result | Finding |\n|---|---|---|---:|---|---|")
for r in sorted([r for r in v if r["design"] in ("RCT-usual-care", "RCT-active") and (r["n_actual"] or 0) >= 200], key=lambda r: -(r["n_actual"] or 0)):
    w("| %s | %s | %s | %d | %s | %s |" % (r["id"], r["population"], r["intervention"][:40], r["n_actual"], r["result"], r["result_summary"].replace("|", "/")[:170]))
w("\nLarge trials that have finished or are about to report and are worth watching: RECOVER-AUTONOMIC IVIG (NCT06305793, n=200) and ivabradine (NCT06305806, n=181), "
  "the probiotic trial PURE-LC (NCT06643299, n=194), the tVNS trial NCT07190105 (n=546), the LIFT trial of pyridostigmine and LDN (NCT06366724, n=160), "
  "the Mount Sinai sirolimus trial (NCT06960928, n=80), the ResetME daratumumab RCT (CTIS 2024-520094-13-00, n=66), and the GLP-1 trial NCT07128082 (n=1,000).\n")

w("## Registries\n")
w("| Registry | Reached? | What was done |\n|---|---|---|")
w("| ClinicalTrials.gov (API v2) | yes | Pulled in full: condition and term searches for ME/CFS and long COVID |")
w("| Europe PMC (REST, core) | yes | Searched for trial publications; abstracts were read for the judgements |")
w("| ISRCTN | yes (XML API) | Searched by condition keyword only: 71 unique records; 21 merged, 22 added, 28 excluded |")
w("| EU CTIS | yes (POST search API; GET returns 403) | 5 keyword searches: 28 records; 10 merged, 6 added, 12 excluded. Status codes are numbers and were not mapped |")
w("| DRKS (Germany) | page loads, but no data | Search results are drawn by JavaScript and are not in the HTML; the guessed API URL returned 404. German trials appear here only when they are also on CT.gov, CTIS or in a paper |")
w("| WHO ICTRP | page loads, but cannot be queried | The search is an ASP.NET form (__VIEWSTATE), so curl cannot run it. Chinese, Iranian, Indian and Japanese registries were reached only through IDs cited in papers |")
w("| ANZCTR | no (HTTP 403) | Not searched |")
w("| ChiCTR | no (HTTP 405) | Not searched; ChiCTR trials appear only when published |")
w("")

# ---------------- coverage check ----------------
cov = {
  "page_out_of_date_or_wrong": [
    {"treatment": "pacing-heart-rate", "issue": "The page says the only controlled test was PACE's APT arm. RECOVER-ENERGIZE structured pacing vs usual care (NCT06404073, n=298, long COVID) reduced how often PEM happened (-0.46, CI -0.65 to -0.27) and how severe it was. It is the first positive randomised trial of pacing.", "trials": ["NCT06404073"]},
    {"treatment": "low-dose-naltrexone", "issue": "The page calls the Nacul trial 'presented May 2026, unpublished'. It is now a medRxiv preprint, Oct 2026 (DOI 10.64898/2026.10.22.26363582; NCT05430152): no difference in fatigue (FSS), pain, symptoms, steps or quality of life at 16 weeks. The cited reference should be updated.", "trials": ["NCT05430152"]},
    {"treatment": "rintatolimod", "issue": "The page cites a 2024 company release for the long COVID trial. ClinicalTrials.gov now has posted results (NCT05592418, n=80): the change in PROMIS fatigue was essentially the same on drug and placebo. The 1994 placebo-controlled RCT (n=92, PMID 8148460, positive on Karnofsky score and treadmill work) is not cited.", "trials": ["NCT05592418", "PMID:8148460"]},
    {"treatment": "tvns", "issue": "PAuSing-pCF is listed as a pending preprint. It is now published in NIHR Open Research (DOI 10.3310/nihropenres.14236.1, n=114): fatigue at 8 weeks did not differ from sham or placebo. It should move from 'pending' into the evidence line. Also pending: NCT07190105 (n=546, finished Dec 2025, not reported).", "trials": ["DOI:10.3310/nihropenres.14236.1", "NCT07190105"]},
    {"treatment": "mito-supplements", "issue": "The page says 'creatine has only a 12-person long-COVID RCT'. That is out of date: a 3-arm placebo-controlled RCT (n=67, PMID 41958909) found that 6 g/day improved fatigue and muscle strength (18 g did not). There are also two other small long COVID creatine RCTs (PMID 38099612, 38684388) and an open-label ME/CFS MRS study (PMID 39408275, n=11). NCT05703074 (nicotinamide riboside plus mind-body therapy, n=310, completed Sep 2024) has not reported.", "trials": ["PMID:41958909", "PMID:37970399", "PMID:38099612", "PMID:38684388", "PMID:39408275", "NCT05703074"]},
    {"treatment": "immunoadsorption", "issue": "The page says plasma exchange has only case reports. There is a sham-controlled plasma exchange RCT (NCT05445674, n=50, posted results): safe, but no difference on any efficacy outcome. That adds to the negative evidence. On BC007, the page says it failed a 114-person unpublished company trial. The database has the 30-person reCOVer crossover published as positive (CTIS 2022-001781-35) and the phase 2 RCT NCT05911009 (n=119, completed Jun 2024) with no posted results and no paper. The page should cite the registry record and make clear the failure is not yet in a publication or posted results.", "trials": ["NCT05445674", "2022-001781-35", "NCT05911009", "NCT05841498", "NCT05710770"]},
    {"treatment": "stimulants", "issue": "Only methylphenidate is cited. Missing: lisdexamfetamine vs placebo (NCT01071044, n=26), which improved executive function and fatigue, and solriamfetol vs placebo (NCT04622293, n=44), which improved fatigue at week 8 but not week 6. Also NCT01966276 (methylphenidate plus nutrient formula, n=134 planned, completed 2014): no results ever reported.", "trials": ["NCT01071044", "NCT04622293", "NCT01966276"]},
    {"treatment": "cbt-supportive", "issue": "The evidence line rests on PACE alone. The database has about 40 other randomised CBT or psychological trials in ME/CFS: FITNET-NHS (n=~240), FINE (n=296, null at 70 weeks), Dutch group, internet and stepped-care CBT, and CBSM, which helped only the high-PEM subgroup. Most are unblinded and use self-rated outcomes, and many are positive. The page's conclusion can stand, but it should say that the wider body exists and why it is discounted (unblinded subjective outcomes, few objective measures).", "trials": ["ISRCTN18020851", "PMID:20418251", "NCT01650636", "PMID:11265953"]},
    {"treatment": "graded-exercise", "issue": "The page could cite MAGENTA (ISRCTN23962803, paediatric, GET no better than activity management) and GETSET (ISRCTN22975026, guided self-help GET, a self-rated gain). PACE alone is not the full trial base.", "trials": ["ISRCTN23962803", "ISRCTN22975026"]},
    {"treatment": "ivig", "issue": "This matches the database. Add the terminated subcutaneous IgPro20 trial in post-COVID POTS (NCT06524739, n=16, uninterpretable). RECOVER-AUTONOMIC IVIG (NCT06305793, n=200) has finished, with no results posted or published.", "trials": ["NCT06524739", "NCT06305793"]},
    {"treatment": "ivabradine", "issue": "The page cites the RECOVER-AUTONOMIC ACC 2026 presentation. The database has no paper or posted results for NCT06305806 yet, so this could not be checked. NCT05481177 (ivabradine nested RCT, n=250, completed Sep 2024) has not reported.", "trials": ["NCT06305806", "NCT05481177"]},
    {"treatment": "steroids", "issue": "There are two unreported long COVID steroid RCTs: prednisolone plus B vitamins (NCT05638633, n=321, completed Jan 2025) and methylprednisolone vs placebo (NCT05986422, terminated, n=96).", "trials": ["NCT05638633", "NCT05986422"]},
    {"treatment": "fmt", "issue": "FMT matches the database. The probiotics entry (one 39-person RCT) is missing several placebo RCTs that have not reported: NCT04741841 (GutMagnific, ME/CFS, n=74, completed 2021), NCT05975034 (long COVID, n=240 planned, completion Dec 2023, status UNKNOWN), NCT05813899 (L. paracasei PS23, terminated). Also PURE-LC (NCT06643299, n=194), completed Jul 2026.", "trials": ["NCT04741841", "NCT05975034", "NCT05813899", "NCT06643299"]},
    {"treatment": "herpes-antivirals", "issue": "The valacyclovir plus celecoxib pilot RCT (NCT06316843, n=59) was completed Oct 2024 and has not reported. The page mentions only the planned trial NCT07597902.", "trials": ["NCT06316843"]},
    {"treatment": "sleep-medicines", "issue": "Add the sodium oxybate placebo RCTs (SAFFE NCT02055898, n=13, more slow-wave sleep but no daytime gain; NCT00498485, terminated, n=17). A second melatonin plus zinc RCT (NCT05454683, n=106 planned, completion Sep 2024) has not reported.", "trials": ["NCT02055898", "NCT00498485", "NCT05454683"]},
    {"treatment": "iv-saline", "issue": "The page says there is no randomised trial. Two have not reported: NCT04740736 (ME/CFS, saline infusion vs placebo in a PEM study, n=106, completed Jul 2025) and NCT02854683 (oral rehydration vs IV saline crossover, status UNKNOWN, completion 2021).", "trials": ["NCT04740736", "NCT02854683"]},
    {"treatment": "eecp", "issue": "This matches the database. NCT05668039 is still not reported (status UNKNOWN, primary completion Apr 2024).", "trials": ["NCT05668039"]},
    {"treatment": "stellate-ganglion-block", "issue": "The page cites STAR-CO NCT06253806, which is not in the database (it was not returned by the ME/CFS or long COVID searches; check it is the right ID). The database has two placebo-controlled RCTs that are recruiting (NCT07374562; a cervicothoracic sympathetic chain block vs saline, NCT07468604) and an open-label dual-block study, NCT05638620 (n=20, COMPASS-31 down 38%).", "trials": ["NCT07374562", "NCT07468604", "NCT05638620"]},
  ],
  "treatments_with_trials_but_no_page_entry": [
    {"treatment": "fluvoxamine", "evidence": "Two placebo RCTs, both favourable: an adaptive fluvoxamine/metformin trial (NCT06128967, n=399), where fluvoxamine reduced fatigue at days 60 and 90, and NCT05874037 (n=137), with a bigger fall in symptom score but no p-value posted. Another RCT is recruiting (NCT07359482).", "result": "positive/mixed", "trials": ["NCT06128967", "NCT05874037", "NCT07359482"]},
    {"treatment": "metformin (treatment of established long COVID)", "evidence": "Null in the adaptive trial (NCT06128967) and in a Korean RCT (KCT0009342, n=396). NCT07280572 (n=1,250) is enrolling.", "result": "null", "trials": ["NCT06128967", "KCT0009342", "NCT07280572"]},
    {"treatment": "colchicine", "evidence": "Null in a placebo RCT (CTRI/2021/11/038234, n=346). In STIMULATE-ICP (open-label) it reduced fatigue slightly. NCT07280572 is ongoing.", "result": "null/mixed", "trials": ["CTRI/2021/11/038234", "ISRCTN10665760"]},
    {"treatment": "temelimab (anti-HERV-W)", "evidence": "Placebo RCT, n=203: no improvement in fatigue.", "result": "null", "trials": ["NCT05497089"]},
    {"treatment": "efgartigimod (FcRn blocker)", "evidence": "Placebo RCT in post-COVID POTS, n=53: no benefit (COMPASS-31 -14.4 vs -16.0). The extension study was terminated.", "result": "null", "trials": ["NCT05633407", "NCT05918978"]},
    {"treatment": "cell therapy (cord blood, MSC)", "evidence": "REGENECYTE cord blood vs placebo (NCT05682560, n=30): positive on fatigue, but small and industry-run. NCT06492798 (n=76) is overdue. Five MSC trials were withdrawn.", "result": "positive (small)", "trials": ["NCT05682560", "NCT06492798"]},
    {"treatment": "Prospekta (Russian ultra-low-dose antibody product)", "evidence": "Placebo RCT for post-COVID asthenia, n=676: FSS fell 16.5 vs 14.0. The difference is significant but small, and the trial was industry-sponsored. This is the largest positive blinded trial in the database.", "result": "positive (small effect)", "trials": ["NCT05074888"]},
    {"treatment": "vortioxetine", "evidence": "Placebo RCT, n=149: no significant gain on the primary cognitive measure (DSST).", "result": "mixed", "trials": ["NCT05047952"]},
    {"treatment": "lithium (low dose)", "evidence": "Placebo RCT, n=52: no effect on fatigue or brain fog.", "result": "null", "trials": ["NCT05618587"]},
    {"treatment": "Pycnogenol and other polyphenols", "evidence": "PYCNOVID (n=153) was null. Smaller positive RCTs: quercetin phytosome (n=78), cocoa flavonoids (FLALOC, n=46), high-polyphenol chocolate (n=10). The page's fibre-polyphenols entry says there is no treatment trial.", "result": "mixed", "trials": ["NCT05890534", "PMID:37716119", "NCT06064838", "PMID:21092175"]},
    {"treatment": "cognitive rehabilitation and tDCS", "evidence": "RECOVER-NEURO (NCT05965739, n=328): no benefit from BrainHQ, PASC-CoRE or tDCS. Single-centre tDCS trials are split (NCT04876417 null; NCT05389592 small gains).", "result": "null", "trials": ["NCT05965739", "NCT04876417", "NCT05389592"]},
    {"treatment": "amantadine", "evidence": "Open-label RCT positive on FSS (PMID 38228731). The placebo RCT NCT06055244 (n=64) finished in Dec 2025 and has not reported. A 1997 ME/CFS RCT found it poorly tolerated.", "result": "mixed", "trials": ["PMID:38228731", "NCT06055244", "PMID:9018019"]},
    {"treatment": "galantamine", "evidence": "Placebo RCT in ME/CFS, n=434: no benefit on any outcome.", "result": "null", "trials": ["PMID:15353532"]},
    {"treatment": "LAU-7b (fenretinide), deupirfenidone", "evidence": "Both null in placebo RCTs of long COVID (n=272 and n=177).", "result": "null", "trials": ["NCT05999435", "NCT04652518"]},
    {"treatment": "structured rehabilitation programmes for long COVID", "evidence": "REGAIN online group rehab (ISRCTN11466448, n=585) improved quality of life. SIPCOV (n=314) and personalised resistance exercise (n=233) were positive, and PEM did not increase. All were unblinded. The page has no long COVID rehab entry separate from GET.", "result": "positive (unblinded)", "trials": ["ISRCTN11466448", "NCT05196451", "NCT04900961"]},
    {"treatment": "plitidepsin", "evidence": "A phase 3 long COVID RCT is registered in CTIS (2023-504087-42-00, n=260 planned). No results yet.", "result": "pending", "trials": ["2023-504087-42-00"]},
  ],
}
json.dump(cov, open(base + "/coverage_check.json", "w"), indent=1, ensure_ascii=False)
w("## Coverage check against the treatments page\n")
w("This compares the database with `site-src/data/treatments.json` (47 treatments, updated 2026-10-28), which was read but not edited. The same content is in `coverage_check.json`.\n")
w("### Page entries whose evidence is out of date, wrong or incomplete\n")
for x in cov["page_out_of_date_or_wrong"]:
    w("- **%s**: %s" % (x["treatment"], x["issue"]))
w("\n### Treatments with trial evidence and no page entry\n")
for x in cov["treatments_with_trials_but_no_page_entry"]:
    w("- **%s** (%s): %s" % (x["treatment"], x["result"], x["evidence"]))
w("\nThese entries match the database and need no change: fludrocortisone, clonidine, rituximab, cyclophosphamide, daratumumab, covid-antivirals, oxaloacetate, coq10-nadh, hbot-clinical, red-light, mast-cell, microclot-therapy, rapamycin, pyridostigmine. "
  "Entries with no trials to check: upright-time, oral-fluids-salt, hbot-home, bpc-157, nicotine-patch, thyroid-testosterone, desmopressin, guanfacine-nac, midodrine (the database holds only a terminated droxidopa study), and glp1 (the trials are all ongoing).\n")

w("## Data-quality caveats\n")
w("- Screening used keyword rules on registry fields, followed by manual review of borderline and long COVID records. Some trials may be wrongly in or out, especially long COVID trials with mixed respiratory and fatigue aims.")
w("- Results were judged from abstracts and posted results only, not full texts. A 'positive' label means the abstract reports benefit on its main outcome; it does not mean the trial was free of bias. Many positive trials are small, unblinded or industry-run.")
w("- %d verified efficacy results come from TCM, herbal, acupuncture or moxibustion trials. These are mostly Chinese, unblinded, and use 'total effective rate' outcomes, and they push up the positive count for the 'supplement' and 'other' classes." % rate([r for r in v if re.search(r"TCM|herbal|acupunct|moxib|qigong|tuina|chinese|granule|decoction", r["intervention"] + r["title"], re.I)])[0])
w("- %d records rest only on preprints (`preprint_only: true`), and 12 publication records have no abstract (result not-assessable)." % sum(r.get("preprint_only") for r in T))
w("- The 'overdue' count may overstate non-reporting. A results paper may exist that Europe PMC did not match to the trial (for example, no registry ID in the abstract, or a non-English paper). UNKNOWN-status CT.gov records past their completion date are counted as overdue.")
w("- ISRCTN coverage is partial (condition keyword search only). CTIS status codes are not mapped. DRKS, ANZCTR, ChiCTR and ICTRP could not be searched, so non-US/EU trials appear only when published or cross-registered.")
w("- Intervention names and classes are normalised with a keyword dictionary (`scripts/build_trials.py`). Model steps are assigned by intervention, not by trial outcome, and combination trials take the first named active drug's class.")
w("- Population 'both' means the trial enrolled or targeted both ME/CFS and long COVID, including post-COVID ME/CFS.")
w("- n_actual comes from CT.gov ACTUAL enrolment, or from the abstract for publication-only records; abstract n is sometimes the number analysed rather than randomised.")
open(base + "/summary.md", "w").write("\n".join(out) + "\n")
print("ok", len(out))
