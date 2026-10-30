"""Build trials.json from CT.gov flats, publications, registries and hand judgements.
Usage: python3 -I build_trials.py <base>"""
import json, sys, re, ast, collections, datetime

base = sys.argv[1].rstrip("/")
W = base + "/work/"
TODAY = datetime.date(2026, 10, 29)
CUTOFF = "2025-10-29"  # completed before this = more than 12 months ago


def L(x):
    if isinstance(x, str) and x[:1] in "[{(":
        try:
            return ast.literal_eval(x)
        except Exception:
            return x
    return x if x is not None else []


def J(f):
    return json.load(open(W + f))


# ---------------- load ----------------
CT = {}
for f in ("ct_mecfs_flat", "ct_lc_flat", "ct_mecfsterm_flat", "ct_lcterm_flat"):
    for r in J(f + ".json"):
        CT[r["nct"]] = r
S = {x["nct"]: x for x in J("screen.json")}
P = {p["pid"]: p for p in J("pubs_selected.json")}
keep = J("epmc_keep.json")
k = keep[1079]
P["P0614"] = {"pmid": k.get("pmid"), "pmcid": k.get("pmcid"), "doi": k.get("doi"), "title": k.get("title"),
              "year": str(k.get("year") or k.get("pubYear") or "2025"), "journal": k.get("journal"),
              "authors": k.get("authors"), "pubtypes": k.get("pubtypes", []), "abstract": k.get("abstract", ""),
              "isme": False, "islc": True, "pid": "P0614", "reg_ids": [], "nct_links": []}
extra = set(J("extra_ct.json"))
manual = J("manual_links.json")
pub_drop = set(J("pub_drop.json"))

JD = {}
for f in ("judg_1.txt", "judg_2.txt", "judg_3.txt", "judg_4.txt"):
    for line in open(W + f):
        line = line.rstrip("\n")
        if not line.strip() or "|" not in line:
            continue
        parts = line.split("|")
        JD[parts[0].strip()] = parts  # later files override earlier

# ---------------- link overrides ----------------
LINK_OVR = {"P0221": "NCT03186027", "P0137": "NCT01071044", "P0230": "NCT05608629", "P0571": "NCT05877508",
            "P0393": "NCT05047952", "P0614": "NCT05035628"}
UNLINK = {"P0114", "P0118", "P0368"}
EXTRA_EXCL = {"NCT05817006": "diagnostic/mechanistic only", "NCT06391970": "diagnostic/mechanistic only",
              "NCT05238415": "diagnostic/mechanistic only"}
for kk, v in JD.items():
    if kk.startswith("NCT") and v[1] == "DROP":
        EXTRA_EXCL[kk] = "dropped on review: " + v[2][:120]
CT_DUP = {}
PUB_DUP = {}
for kk, v in JD.items():
    if v[1] != "DUP":
        continue
    m = re.findall(r"(NCT\d{8}|P\d{4})", v[2])
    if not m:
        continue
    (CT_DUP if kk.startswith("NCT") else PUB_DUP)[kk.replace("PUB:", "")] = m[0]
PLATFORM = {"NCT06305780": ["NCT06305793", "NCT06305806"], "NCT06404047": ["NCT06404073", "NCT06404060"]}

for pid, p in P.items():
    links = set(p.get("nct_links") or [])
    if pid in manual:
        links.add(manual[pid])
    if pid in LINK_OVR:
        links = {LINK_OVR[pid]}
    if pid in UNLINK:
        links = set()
    p["links"] = sorted(n for n in links if n in CT)

# ---------------- included NCT set ----------------
inc = {n for n, x in S.items() if x["decision"] == "in"} | extra
excl = collections.Counter()
excl_detail = collections.defaultdict(list)
for n, x in S.items():
    if x["decision"] == "out" and n not in extra:
        excl[x["reason"]] += 1
forced = []
for pid, p in P.items():
    if pid in pub_drop or JD.get("PUB:" + pid, [None, ""])[1] == "DROP":
        continue
    for n in p["links"][:1]:
        if n not in inc and n not in EXTRA_EXCL:
            inc.add(n)
            forced.append((n, pid))
            if n in S and S[n]["decision"] == "out":
                excl[S[n]["reason"]] -= 1
for n in list(inc):
    if n in EXTRA_EXCL:
        inc.discard(n)
        excl[EXTRA_EXCL[n] if "dropped" not in EXTRA_EXCL[n] else "dropped on manual review (not a treatment trial or out of scope)"] += 1
    elif n in CT_DUP:
        inc.discard(n)
for n in PLATFORM:
    inc.discard(n)
for n in [n for n in inc if n not in CT]:
    print("WARN not in CT flats", n)
    inc.discard(n)

# ---------------- intervention dictionary ----------------
IM, AV, AU, ME, SU, BR, DV, AP, MB, OT = ("immunomodulator", "antiviral", "autonomic/vasoactive",
    "metabolic/mitochondrial", "supplement", "behavioural/rehab", "device/neuromodulation",
    "apheresis/immunoadsorption", "microbiome", "other")
D = [
    (r"rintatolimod|ampligen|poly ?i[:.]? ?c", "rintatolimod", IM, [1, 7]),
    (r"\bbc ?-?007|rovunaptabin", "BC007 (rovunaptabin)", IM, [1, 4]),
    (r"immuno-?adsorption|plasmapheres|plasma exchange|apheresis|h\.e\.l\.p|inuspheresis", "immunoadsorption/apheresis", AP, [1, 4]),
    (r"valganciclovir", "valganciclovir", AV, [1]),
    (r"valacy|valaciclovir", "valacyclovir", AV, [1]),
    (r"\bacyclovir|\baciclovir", "acyclovir", AV, [1]),
    (r"famciclovir", "famciclovir", AV, [1]),
    (r"nirmatrelvir|paxlovid", "nirmatrelvir/ritonavir", AV, [1]),
    (r"ensitrelvir", "ensitrelvir", AV, [1]),
    (r"remdesivir", "remdesivir", AV, [1]),
    (r"molnupiravir", "molnupiravir", AV, [1]),
    (r"tenofovir|emtricitabine|truvada", "tenofovir/emtricitabine", AV, [1]),
    (r"plitidepsin", "plitidepsin", AV, [1]),
    (r"ivermectin", "ivermectin", AV, [1]),
    (r"maraviroc", "maraviroc", IM, [1, 6]),
    (r"leronlimab", "leronlimab", IM, [1, 6]),
    (r"inosine pranobex|isoprinosine", "inosine pranobex", IM, [1]),
    (r"immunoglobulin|\bivig\b|\bscig\b|gamunex|privigen|octagam|hizentra", "immunoglobulin (IVIG/SCIG)", IM, [1]),
    (r"rituximab", "rituximab", IM, [1]),
    (r"cyclophosphamide", "cyclophosphamide", IM, [1]),
    (r"daratumumab", "daratumumab", IM, [1]),
    (r"efgartigimod", "efgartigimod", IM, [1]),
    (r"temelimab", "temelimab", IM, [1, 7]),
    (r"tocilizumab", "tocilizumab", IM, [6, 7]),
    (r"anakinra", "anakinra", IM, [6, 7]),
    (r"baricitinib|tofacitinib|ruxolitinib", "JAK inhibitor", IM, [6, 7]),
    (r"abatacept", "abatacept", IM, [1]),
    (r"infliximab|etanercept|adalimumab|anti-?tnf", "anti-TNF", IM, [6, 7]),
    (r"larazotide", "larazotide", IM, [1, 7]),
    (r"naltrexone|\bldn\b", "low-dose naltrexone", IM, [6, 7]),
    (r"colchicine", "colchicine", IM, [6, 7]),
    (r"rapamycin|sirolimus", "rapamycin (sirolimus)", IM, [1, 6]),
    (r"hydrocortisone", "hydrocortisone", OT, [2]),
    (r"fludrocortisone|florinef", "fludrocortisone", AU, [3]),
    (r"prednis|dexamethasone|methylpred|corticoster|budesonide", "corticosteroid", IM, [6, 7]),
    (r"mesenchymal|stem cell|\bmsc\b|cord blood|exosome|regenecyte|stromal", "cell therapy (MSC/cord blood)", IM, [6]),
    (r"montelukast", "montelukast", IM, [6, 7]),
    (r"famotidine|loratadine|cetirizine|fexofenadine|antihistamin|cromolyn|cromoglic|ketotifen|rupatadine|desloratadine|mast cell", "antihistamine/mast-cell stabiliser", IM, [6, 7]),
    (r"interferon|staphylococc|staphypan", "immune stimulant (interferon/staph toxoid)", IM, [1]),
    (r"metformin", "metformin", ME, [5, 6]),
    (r"abrocitinib|upadacitinib", "JAK inhibitor", IM, [6, 7]),
    (r"inebilizumab", "inebilizumab", IM, [1]),
    (r"n-803|anktiva|il-15", "N-803 (IL-15 superagonist)", IM, [1]),
    (r"keytruda|pembrolizumab", "pembrolizumab", IM, [1]),
    (r"aer002|vyd2311|monoclonal antibod", "SARS-CoV-2 monoclonal antibody", AV, [1]),
    (r"ruconest|c1.?esterase", "C1-esterase inhibitor", IM, [1, 4]),
    (r"rslv-132", "RSLV-132", IM, [1]),
    (r"ne3107", "NE3107", IM, [6, 7]),
    (r"sonlicromanol", "sonlicromanol", ME, [5, 6]),
    (r"dapagliflozin|empagliflozin", "SGLT2 inhibitor", ME, [5, 6]),
    (r"meldonium|ethylmethylhydroxypyridine|ethyl methyl hydroxypyridine", "meldonium/EMHPS", ME, [5, 6]),
    (r"taurine", "taurine", ME, [5, 6]),
    (r"n-acetyl ?cysteine|\bnac\b", "N-acetylcysteine", SU, [6]),
    (r"pirfenidone|lyt-100", "(deu)pirfenidone", OT, [6]),
    (r"apabetalone", "apabetalone", OT, [4, 6]),
    (r"tonabersat", "tonabersat", OT, [2]),
    (r"fampridine", "fampridine", OT, [2]),
    (r"somatropin|growth hormone", "growth hormone", OT, [2]),
    (r"\bdhea\b", "DHEA", OT, [2]),
    (r"ondansetron|granisetron|tropisetron|5-ht3", "5-HT3 antagonist", OT, [2]),
    (r"osu ?6162", "OSU6162", OT, [2]),
    (r"doxycycline|azithromycin|antibiotic", "antibiotic", OT, [1]),
    (r"lau-?7b|fenretinide", "LAU-7b (fenretinide)", OT, [6, 7]),
    (r"cistanche", "herbal/TCM", SU, []),
    (r"tianeptine|venlafaxine|imipramine|nefazodone|tricyclic", "antidepressant", OT, [2]),
    (r"ivabradine", "ivabradine", AU, [8]),
    (r"guanfacine", "guanfacine (+/- NAC)", AU, [2, 8]),
    (r"clonidine", "clonidine", AU, [8]),
    (r"propranolol|metoprolol|bisoprolol|beta.?block|nebivolol|atenolol", "beta-blocker", AU, [8]),
    (r"stellate ganglion|sympathetic (chain )?block|cervicothoracic sympathetic", "stellate ganglion/sympathetic block", AU, [8]),
    (r"pyridostigmine|mestinon", "pyridostigmine", AU, [3, 5]),
    (r"midodrine", "midodrine", AU, [3, 4]),
    (r"droxidopa", "droxidopa", AU, [3, 4]),
    (r"desmopressin", "desmopressin", AU, [3]),
    (r"\bsaline\b|intravenous fluid|\biv fluid|volume expan", "IV saline/fluids", AU, [3]),
    (r"compression garment|compression stocking|abdominal binder|compression", "compression garments", AU, [3]),
    (r"sildenafil|tadalafil", "PDE5 inhibitor", AU, [4, 5]),
    (r"vericiguat|riociguat", "sGC stimulator", AU, [4, 5]),
    (r"counterpulsation|\beecp\b", "EECP", DV, [4, 5]),
    (r"pentoxifylline", "pentoxifylline", AU, [4]),
    (r"apixaban|clopidogrel|anticoag|heparin|enoxaparin|nattokinase|lumbrokinase|aspirin|microclot|triple therapy", "anticoagulant/antiplatelet", AU, [4]),
    (r"nitrate|beetroot", "dietary nitrate", SU, [4, 5]),
    (r"arginine|citrulline", "L-arginine/citrulline", SU, [4, 5]),
    (r"aripiprazole|abilify", "low-dose aripiprazole", OT, [2]),
    (r"lisdexamfetamine|methylphenidate|modafinil|armodafinil|amphetamine|solriamfetol|stimulant", "stimulant", OT, [2]),
    (r"fluvoxamine", "fluvoxamine", OT, [2]),
    (r"vortioxetine", "vortioxetine", OT, [2]),
    (r"duloxetine|sertraline|fluoxetine|citalopram|escitalopram|bupropion|amitriptyline|nortriptyline|phenelzine|moclobemide|antidepress|mirtazapine|agomelatine|trazodone|selegiline", "antidepressant", OT, [2]),
    (r"amisulpride", "amisulpride", OT, [2]),
    (r"lithium", "lithium", OT, [2]),
    (r"galantamine|donepezil|memantine", "cholinesterase inhibitor/memantine", OT, [2]),
    (r"amantadine", "amantadine", OT, [2]),
    (r"prospekta", "Prospekta", OT, [2]),
    (r"ketamine", "ketamine", OT, [2]),
    (r"cannabi|\bcbd\b|dronabinol", "cannabinoid", OT, [2]),
    (r"rimonabant", "rimonabant", OT, [2, 7]),
    (r"pregabalin|gabapentin", "gabapentinoid", OT, [8]),
    (r"sodium oxybate|xyrem|zolpidem|melatonin|suvorexant|daridorexant|lemborexant|insomnia|sleep", "sleep treatment", OT, [8]),
    (r"oxaloacetate|benagene", "oxaloacetate", ME, [5, 6]),
    (r"coenzyme ?q|coq ?10|ubiquin", "coenzyme Q10", ME, [5, 6]),
    (r"\bnadh\b|nicotinamide|\bnad\+?\b|niagen", "NADH/NAD+ precursor", ME, [5, 6]),
    (r"creatine|guanidinoacet", "creatine/GAA", ME, [5, 6]),
    (r"ribose", "D-ribose", ME, [5, 6]),
    (r"carnitine", "carnitine", ME, [5, 6]),
    (r"mitoquinone|mitoq", "mitoquinone", ME, [5, 6]),
    (r"axa ?1125|\baxa\b", "AXA1125", ME, [5, 6]),
    (r"molecular hydrogen|hydrogen[- ]rich|hydrogen water|hydrogen inhal|hydrogen gas", "molecular hydrogen", ME, [5, 6]),
    (r"semaglutide|tirzepatide|liraglutide|glp-?1", "GLP-1 agonist", ME, [6, 7]),
    (r"omega-?3|fish oil|\bepa\b|\bdha\b|polyphenol|resveratrol|quercetin|curcumin|flavon|cocoa|pycnogenol|pine bark|luteolin|palmitoylethanolamide|\bpea\b", "polyphenol/omega-3/anti-inflammatory supplement", SU, [6, 7]),
    (r"faecal transplant|fecal transplant|probiotic|synbiotic|prebiotic|lactobac|bifido|faecal microbiota|fecal microbiota|microbiota transplant|\bfmt\b|inulin|sim01|rifaximin|\bfib(re|er)\b", "microbiome therapy", MB, [1, 7]),
    (r"magnesium", "magnesium", SU, []),
    (r"vitamin|ascorbic|thiamine|cobalamin|\bb12\b|folic|folate|\bzinc\b|multivitamin|\biron\b|selenium", "vitamins/minerals", SU, []),
    (r"ginseng|ginkgo|herbal|traditional chinese|\btcm\b|decoction|granule|capsule.*chinese|chinese medicine|kampo|ayurved|astragal|\bqi\b|\btang\b|\bwan\b", "herbal/TCM", SU, []),
    (r"homeopath", "homeopathy", OT, []),
    (r"hyperbaric|\bhbot\b", "hyperbaric oxygen", OT, [5, 6]),
    (r"photobiomodulation|red light|near-infrared|vielight|low-level laser", "photobiomodulation", DV, [5]),
    (r"vagus|vagal|\btvns\b|\btavns\b", "vagus nerve stimulation", DV, [7, 8]),
    (r"\btdcs\b|transcranial direct|\brtms\b|transcranial magnetic|\btacs\b|transcranial|neuromodulation|neurofeedback", "transcranial stimulation/neurofeedback", DV, [2]),
    (r"biofeedback", "biofeedback", DV, [8]),
    (r"lightning process|brain retrain|gupta|\bdnrs\b|amygdala", "brain-retraining programme", BR, [2, 8]),
    (r"cognitive behavio|\bcbt\b|acceptance and commitment|mindfulness|psychotherap|counsel|psychoeducation|psycholog|self-management|coaching", "CBT/psychological therapy", BR, [2]),
    (r"cognitive rehab|cognitive training|brain training|computeri[sz]ed cognitive|brainhq|cognitive remediation", "cognitive rehabilitation", BR, [2]),
    (r"virtual reality|immersive", "immersive VR", BR, [2]),
    (r"pacing|energy management|activity management|energy envelope", "pacing/energy management", BR, [6]),
    (r"inspiratory muscle|respiratory muscle|breathing", "breathing/respiratory training", BR, []),
    (r"graded exercise|\bget\b", "graded exercise therapy", BR, []),
    (r"exercise|rehabilitat|physiotherap|aerobic|resistance training|yoga|tai chi|qigong|pilates|walking|physical activity|training", "exercise/rehabilitation", BR, []),
    (r"acupunct|moxibust|cupping|massage|osteopath|manual therap|chiropract", "acupuncture/manual therapy", OT, []),
    (r"ketogenic|fasting|\bdiet\b|dietary|nutrition", "diet", SU, []),
    (r"nicotine", "nicotine patch", OT, [7]),
    (r"thyroid|testosterone|levothyrox|liothyron", "thyroid/testosterone", OT, [2]),
    (r"bpc-?157", "BPC-157", OT, [6]),
    (r"statin\b|atorvastatin|rosuvastatin", "statin", OT, [6, 7]),
    (r"ozone", "ozone therapy", OT, []),
]
DC = [(re.compile(a, re.I), b, c, d) for a, b, c, d in D]
CTRL = re.compile(r"^(placebo|sham|usual care|standard care|standard of care|control|wait ?-?list|no interventions?|treatment as usual|tau)\b|placebo$|^sham", re.I)
TYPE_CLASS = {"DRUG": OT, "BIOLOGICAL": IM, "DEVICE": DV, "BEHAVIORAL": BR, "DIETARY_SUPPLEMENT": SU,
              "PROCEDURE": OT, "OTHER": OT, "RADIATION": OT, "GENETIC": OT, "COMBINATION_PRODUCT": OT}


def classify(texts, fallback_type=None, fallback_name=None):
    names, cls, steps = [], None, set()
    for t in texts:
        if not t:
            continue
        for rx, nm, c, st in DC:
            if rx.search(t) and nm not in names:
                names.append(nm)
                cls = cls or c
                steps |= set(st)
                break
    if len(names) > 1 and "IV saline/fluids" in names:
        names.remove("IV saline/fluids")
        steps = set()
        for nm2 in names:
            steps |= set(next(st for rx, nm, c, st in DC if nm == nm2))
        cls = next(c for rx, nm, c, st in DC if nm == names[0])
    if not names:
        return (fallback_name or "unclassified"), TYPE_CLASS.get(fallback_type or "", OT), []
    return " + ".join(names[:3]), cls, sorted(steps)


def ct_interventions(r):
    iv = [x for x in L(r.get("interventions")) if isinstance(x, (list, tuple))]
    act = [x for x in iv if not CTRL.search((x[1] or "").strip())]
    act = act or iv
    texts = [x[1] for x in act]
    name, cls, st = classify(texts, act[0][0] if act else None, "; ".join(texts[:2]) if texts else None)
    if name == (("; ".join(texts[:2])) if texts else None) or name == "unclassified":
        # try descriptions then title
        n2, c2, s2 = classify([x[2] for x in act if len(x) > 2] + [r.get("title", "")])
        if n2 != "unclassified":
            return n2, c2, s2, texts
    return name, cls, st, texts


# ---------------- field derivations ----------------
def design_ct(r):
    model = (r.get("model") or "").upper()
    alloc = (r.get("allocation") or "").upper()
    arms = [x for x in L(r.get("arms")) if isinstance(x, (list, tuple))]
    n = r.get("enrollment")
    if "CROSSOVER" in model:
        return "crossover"
    if "SINGLE" in model or alloc in ("NA", "", "NONE"):
        if "SINGLE" in model or len(arms) <= 1:
            return "n-of-1" if n == 1 else "single-arm"
    if alloc == "NON_RANDOMIZED":
        return "open-label"
    types = {a[0] for a in arms}
    labels = " ".join(a[1] or "" for a in arms).lower()
    if types & {"PLACEBO_COMPARATOR", "SHAM_COMPARATOR"} or re.search(r"placebo|sham", labels + " " + (r.get("title") or "").lower()):
        return "RCT-placebo"
    if "NO_INTERVENTION" in types or re.search(r"usual|standard|wait|treatment as usual|no treatment|no intervention", labels):
        return "RCT-usual-care"
    return "RCT-active"


def severity(text):
    if not text:
        return "not stated"
    low = text.lower()
    parts = re.split(r"exclusion criteria", low, maxsplit=1)
    inc_t, exc_t = parts[0], (parts[1] if len(parts) > 1 else "")
    rx = r"housebound|house-bound|bedbound|bed-bound|bed bound|severe(ly)? (me|cfs|affected|ill|disease)|very severe|severe myalgic"
    if re.search(rx, inc_t):
        return "severe included"
    if re.search(rx, exc_t):
        return "severe excluded"
    return "not stated"


MEX = re.compile(r"myalgic|chronic fatigue|\bcfs\b|me/cfs|seid|systemic exertion", re.I)
LCX = re.compile(r"long.?covid|post.?covid|post-acute|pasc|sars-cov-2|covid", re.I)


def pop_from(isme, islc):
    return "both" if isme and islc else ("ME/CFS" if isme else ("long COVID" if islc else "unclear"))


def pop_ct(n, r):
    if n in S:
        return pop_from(S[n]["isme"], S[n]["islc"])
    t = " ".join([r.get("title") or "", " ".join(L(r.get("conditions")) or [])])
    return pop_from(bool(MEX.search(t)), bool(LCX.search(t)))


RES = {"POS": "positive", "NULL": "null", "NEG": "negative/harm", "MIX": "mixed", "FEAS": "feasibility-only"}
POP_OVR = {"NCT06786832": "ME/CFS"}


def derive_status_result(status, pcomp):
    st = (status or "").upper()
    if st == "WITHDRAWN":
        return "withdrawn"
    if st in ("COMPLETED", "TERMINATED", "UNKNOWN") or (st in ("ACTIVE_NOT_RECRUITING", "SUSPENDED") and False):
        if pcomp and pcomp[:7] < CUTOFF[:7]:
            return "unreported-after-completion"
        return "pending"
    return "pending"


def pubinfo(pid):
    p = P[pid]
    return {"pid": pid, "pmid": p.get("pmid"), "doi": p.get("doi"), "year": p.get("year"), "title": p.get("title"),
            "journal": p.get("journal"),
            "preprint": bool(re.search(r"rxiv|research square|preprint", p.get("journal") or "", re.I) or re.match(r"10\.(1101|21203|64898)/", p.get("doi") or "")),
            "protocol": bool(re.search(r"protocol|design and rationale|rationale and design", p.get("title") or "", re.I))}


def best_pub(pids):
    ps = [pubinfo(x) for x in pids]
    np_ = [x for x in ps if not x["protocol"]] or ps
    np_.sort(key=lambda x: (x["year"] or "0"), reverse=True)
    return np_[0] if np_ else None


def judge_for(keys):
    """Return (code, summary, source_key) preferring a non-NORES judgement."""
    got = [(kk, JD[kk]) for kk in keys if kk in JD and JD[kk][1] not in ("DUP", "DROP")]
    for kk, v in got:
        if v[1] in RES:
            return v[1], v[2], kk
    for kk, v in got:
        return v[1], v[2], kk
    return None, None, None


# ---------------- build CT records ----------------
pubs_by_nct = collections.defaultdict(list)
for pid, p in P.items():
    if pid in pub_drop or JD.get("PUB:" + pid, [None, ""])[1] == "DROP" or pid in PUB_DUP:
        continue
    if p["links"]:
        n = p["links"][0]
        n = CT_DUP.get(n, n)
        pubs_by_nct[n].append(pid)
for d_, t_ in PUB_DUP.items():
    pass  # handled below

judg3 = {l.split("|")[0] for l in open(W + "judg_3.txt") if "|" in l}
records = {}
for n in sorted(inc):
    r = CT[n]
    dups = [d for d, t in CT_DUP.items() if t == n]
    all_ids = [n] + dups + [x for x in L(r.get("secondary_ids")) or [] if x] + ([r["org_id"]] if r.get("org_id") else [])
    iname, icls, isteps, raw_iv = ct_interventions(r)
    pids = pubs_by_nct.get(n, [])
    keys = [n] + dups + ["PUB:" + x for x in pids]
    code, summ, src = judge_for(keys)
    enr, et = r.get("enrollment"), (r.get("enroll_type") or "").upper()
    status = r.get("status")
    pcomp = r.get("completion")
    if code in RES:
        result = RES[code]
    else:
        result = derive_status_result(status, pcomp)
        if code == "NORES" and src and not src.startswith("PUB:") and str(r.get("has_results")) == "True":
            result = "not-assessable"
    bp = best_pub(pids) if pids else None
    if code in RES:
        if src.startswith("PUB:"):
            sp = P[src[4:]]
            source = "abstract PMID %s" % sp.get("pmid") if sp.get("pmid") else "abstract DOI %s" % sp.get("doi")
        elif src in judg3 or str(r.get("has_results")) == "True" and not pids:
            source = "ClinicalTrials.gov posted results"
        else:
            source = ("abstract PMID %s" % bp["pmid"]) if bp and bp.get("pmid") else (("abstract DOI %s" % bp["doi"]) if bp and bp.get("doi") else ("ClinicalTrials.gov posted results" if str(r.get("has_results")) == "True" else "publication abstract"))
        if str(r.get("has_results")) == "True" and "posted" not in source:
            source += "; ClinicalTrials.gov results posted"
        verified = True
    else:
        source = "registry status only"
        verified = False
        if summ is None:
            summ = {"withdrawn": "Withdrawn before enrolment", "pending": "No results yet",
                    "unreported-after-completion": "Primary completion %s; no results posted or published abstract found" % pcomp}[result]
    masking = (r.get("masking") or "").upper()
    phases = L(r.get("phases")) or []
    rec = {
        "id": n, "all_ids": list(dict.fromkeys(all_ids)), "title": r.get("title"), "acronym": r.get("acronym"),
        "intervention": iname, "intervention_raw": raw_iv, "intervention_class": icls,
        "population": POP_OVR.get(n) or pop_ct(n, r), "severity_included": severity(r.get("criteria")),
        "design": design_ct(r), "blinded": bool(masking and masking != "NONE"), "masking": masking or None,
        "phase": "/".join(phases) if isinstance(phases, list) else phases,
        "n_planned": enr if et == "ESTIMATED" else None, "n_actual": enr if et == "ACTUAL" else None,
        "start": r.get("start"), "primary_completion": pcomp, "completion": r.get("completion_full") or pcomp,
        "status": status, "why_stopped": r.get("why_stopped"),
        "country": L(r.get("countries")) or [], "sponsor": r.get("sponsor"),
        "primary_outcome": "; ".join(x[0] for x in (L(r.get("primary_outcomes")) or []) if isinstance(x, (list, tuple)))[:500],
        "result": result, "result_summary": summ, "result_source": source,
        "doi": bp["doi"] if bp else None, "pmid": bp["pmid"] if bp else None,
        "publications": [pubinfo(x) for x in pids],
        "has_results_posted": str(r.get("has_results")) == "True",
        "model_steps": isteps, "verified": verified, "registry": "ClinicalTrials.gov",
        "url": "https://clinicaltrials.gov/study/" + n,
    }
    if r.get("enrollment") is not None and et == "ACTUAL" and r.get("results_n"):
        rec["n_actual"] = r.get("results_n") or enr
    records[n] = rec

# platform parent ids onto children
for parent, kids in PLATFORM.items():
    for kd in kids:
        if kd in records:
            records[kd]["all_ids"].append(parent)
            records[kd]["platform"] = "RECOVER" if parent == "NCT06305780" else "RECOVER-ENERGIZE"

# ---------------- publication-only records ----------------
for pid, p in P.items():
    if pid in pub_drop or pid in PUB_DUP or p["links"]:
        continue
    j = JD.get("PUB:" + pid)
    if not j or j[1] == "DROP":
        if not j:
            print("WARN pub without judgement", pid, (p.get("title") or "")[:70])
        continue
    dups = [d for d, t in PUB_DUP.items() if t == pid]
    allp = [pid] + dups
    reg = [x for x in (p.get("reg_ids") or []) if not x.startswith("NCT")]
    if pid == "P0368":
        reg = ["NCT04960215", "2020-005961-16"] + reg
    rid = reg[0] if reg else ("PMID:" + p["pmid"] if p.get("pmid") else "DOI:" + (p.get("doi") or pid))
    base_rid, kq = rid, 2
    while rid in records:
        rid = "%s (paper %d)" % (base_rid, kq); kq += 1
    jint, jcls, jdes, jn, jpop = (j + [None] * 8)[3:8]
    name, cls, st = classify([jint or "", p.get("title") or ""])
    if name == "unclassified":
        name, cls = jint or "unclassified", jcls or OT
        st = []
    if jint and name != "unclassified" and len(jint) < 60 and " + " not in name:
        # keep hand-normalised name but regex class/steps
        name = jint
    ab = (p.get("abstract") or "").lower()
    code = j[1]
    res = RES.get(code)
    if res is None:
        res = "pending" if re.search(r"protocol", j[2], re.I) else "not-assessable"
    try:
        nn = int(re.sub(r"[^\d]", "", jn or "")) if jn else None
    except ValueError:
        nn = None
    records[rid] = {
        "id": rid, "all_ids": list(dict.fromkeys([rid] + reg + ["PMID:" + P[x]["pmid"] for x in allp if P[x].get("pmid")])),
        "title": p.get("title"), "acronym": None, "intervention": name, "intervention_raw": [jint] if jint else [],
        "intervention_class": cls, "population": jpop or pop_from(p.get("isme"), p.get("islc")),
        "severity_included": severity(ab) if re.search(r"housebound|bedbound|bed-bound|severe", ab) else "not stated",
        "design": jdes or "unclear", "blinded": bool(re.search(r"blind|sham|masked", ab)), "masking": None, "phase": None,
        "n_planned": None, "n_actual": nn, "start": None, "primary_completion": None, "completion": None,
        "status": "published", "why_stopped": None, "country": [], "sponsor": None, "primary_outcome": None,
        "result": res, "result_summary": j[2],
        "result_source": ("abstract PMID %s" % p["pmid"]) if p.get("pmid") else ("abstract DOI %s" % p.get("doi")),
        "doi": p.get("doi"), "pmid": p.get("pmid"), "publications": [pubinfo(x) for x in allp],
        "has_results_posted": False, "model_steps": st, "verified": code in RES and "no abstract" not in j[2].lower(),
        "registry": "none found" if not reg else "other registry", "year": p.get("year"),
        "url": ("https://europepmc.org/abstract/MED/" + p["pmid"]) if p.get("pmid") else ("https://doi.org/" + p["doi"] if p.get("doi") else None),
    }
    if "no abstract" in j[2].lower():
        records[rid]["result_source"] = "title only; result not verified"

# PUB DUPs that point into CT-linked pubs: attach
for d_, t_ in PUB_DUP.items():
    for r in records.values():
        if any(x["pid"] == t_ for x in r["publications"]) and not any(x["pid"] == d_ for x in r["publications"]):
            r["publications"].append(pubinfo(d_))

# ---------------- registries: ISRCTN, CTIS ----------------
idx = {}
for kk, r in records.items():
    for a in r["all_ids"]:
        idx[a.upper()] = kk
    for pb in r["publications"]:
        for rg in P[pb["pid"]].get("reg_ids") or []:
            idx.setdefault(rg.upper(), kk)
pid2rec = {pb["pid"]: kk for kk, r in records.items() for pb in r["publications"]}

ISR_MERGE = {"ISRCTN77567702": "P0173", "ISRCTN13629437": "P0296", "ISRCTN15414370": "P0446", "ISRCTN16033549": "P0539",
             "ISRCTN18015802": "P0299", "ISRCTN59878666": "NCT00893438", "ISRCTN72940450": "NCT05911906"}
ISR_ADD = ["ISRCTN" + x for x in "16025168 42784463 74088209 58998320 16045700 07518149 31632033 31455243 21288417 37461058 16012321 33340595 46454974 39864861 11563542 15068418 12368131 11868601 36407216 15931869 11050086 15106165".split()]
CTIS_MERGE = {"2024-516174-30-01": "NCT05967052", "2024-517713-33-00": "NCT05697640", "2024-510790-19-00": "NCT06012552",
              "2024-516378-31-00": "NCT06766825", "2024-512500-19-00": "P0287", "2025-521866-92-01": "NCT07298005",
              "2024-519540-32-00": "NCT07359482", "2023-508744-22-00": "NCT06524739", "2023-504826-21-00": "NCT05986422",
              "2025-522611-42-00": "NCT06928272"}
CTIS_ADD = ["2025-524443-12-00", "2024-511628-16-00", "2024-517199-39-00", "2023-507870-42-00", "2023-504087-42-00", "2024-520094-13-00"]


def target_key(t):
    if t.startswith("P") and not t.startswith("PMID"):
        return pid2rec.get(t)
    return t if t in records else None


def attach(regid, tgt):
    kk = target_key(tgt) if tgt else idx.get(regid.upper())
    if kk and kk in records:
        if regid not in records[kk]["all_ids"]:
            records[kk]["all_ids"].append(regid)
        return kk
    return None


reg_log = {"ISRCTN": collections.Counter(), "CTIS": collections.Counter()}
for r in J("isrctn_flat.json"):
    i = r["id"]
    if attach(i, ISR_MERGE.get(i)):
        reg_log["ISRCTN"]["merged into existing record"] += 1
        continue
    if i not in ISR_ADD:
        reg_log["ISRCTN"]["excluded (observational, diagnostic, mechanistic, breathlessness/mental-health only, or off-topic)"] += 1
        continue
    txt = " ".join([r.get("title", ""), r.get("conditions", "")])
    nm, cls, st = classify([r.get("intervention", ""), r.get("title", "")], None, None)
    end = r.get("end") or ""
    res = "unreported-after-completion" if end and end < CUTOFF else "pending"
    if (r.get("final_n") in ("0",) and end and end < CUTOFF):
        res = "unreported-after-completion"
    des = r.get("design") or ""
    dl = des.lower()
    design = ("crossover" if "cross" in dl else "RCT-placebo" if "placebo" in dl else
              "RCT-unspecified" if "random" in dl else "single-arm" if "single" in dl or "non-random" in dl else "unclear")
    records[i] = {
        "id": i, "all_ids": [i] + [x for x in r.get("ext_refs") or [] if isinstance(x, str)], "title": r.get("scientific_title") or r["title"],
        "acronym": r.get("acronym"), "intervention": nm if nm != "unclassified" else (r.get("acronym") or r["title"])[:80],
        "intervention_raw": [(r.get("intervention") or "")[:300]], "intervention_class": cls,
        "population": pop_from(bool(MEX.search(txt)), bool(LCX.search(txt))), "severity_included": "not stated",
        "design": design, "blinded": bool(re.search(r"blind", dl)), "masking": None, "phase": None,
        "n_planned": r.get("target_n"), "n_actual": r.get("final_n") if r.get("final_n") not in (None, "", "0") else None,
        "start": r.get("start"), "primary_completion": end, "completion": end,
        "status": r.get("overall_status_hint") or "see registry", "why_stopped": None,
        "country": sorted(set(r.get("countries") or [])), "sponsor": r.get("sponsor"),
        "primary_outcome": (r.get("primary_outcome") or "")[:500], "result": res,
        "result_summary": "Registry record only; no results read" + ("; outputs listed on registry" if r.get("outputs") else ""),
        "result_source": "registry status only", "doi": None, "pmid": None, "publications": [],
        "has_results_posted": False, "model_steps": st, "verified": False, "registry": "ISRCTN",
        "url": "https://www.isrctn.com/" + i, "registry_design": des,
    }
    reg_log["ISRCTN"]["added as new record"] += 1

for r in J("ctis_flat.json"):
    i = r["ctNumber"]
    if attach(i, CTIS_MERGE.get(i)):
        reg_log["CTIS"]["merged into existing record"] += 1
        continue
    if i not in CTIS_ADD:
        reg_log["CTIS"]["excluded (imaging, dyspnoea, acute COVID, prevention, off-topic)"] += 1
        continue
    txt = " ".join([r.get("ctTitle", ""), r.get("conditions", "")])
    nm, cls, st = classify([r.get("product", ""), r.get("ctTitle", "")])
    t = r.get("ctTitle", "").lower()
    design = "RCT-placebo" if "placebo" in t else ("RCT-unspecified" if "random" in t else "unclear")
    def iso(d):
        m = re.match(r"(\d\d)/(\d\d)/(\d{4})", d or "")
        return "%s-%s-%s" % (m.group(3), m.group(2), m.group(1)) if m else None
    records[i] = {
        "id": i, "all_ids": [i, r.get("shortTitle")] if r.get("shortTitle") else [i], "title": r.get("ctTitle"),
        "acronym": r.get("shortTitle"), "intervention": nm if nm != "unclassified" else (r.get("product") or "")[:80],
        "intervention_raw": [r.get("product")], "intervention_class": cls,
        "population": pop_from(bool(MEX.search(txt)), bool(LCX.search(txt))), "severity_included": "not stated",
        "design": design, "blinded": "blind" in t, "masking": None, "phase": r.get("trialPhase"),
        "n_planned": r.get("totalNumberEnrolled"), "n_actual": None, "start": iso(r.get("startDateEU")),
        "primary_completion": None, "completion": None,
        "status": "CTIS status code %s (unmapped)" % r.get("ctStatus"), "why_stopped": None,
        "country": [c.split(":")[0] for c in r.get("trialCountries") or []], "sponsor": r.get("sponsor"),
        "primary_outcome": (r.get("primaryEndPoint") or "")[:500], "result": "pending",
        "result_summary": "Registry record only (EU CTIS); results first received: %s" % r.get("resultsFirstReceived"),
        "result_source": "registry status only", "doi": None, "pmid": None, "publications": [],
        "has_results_posted": r.get("resultsFirstReceived") == "Yes", "model_steps": st, "verified": False,
        "registry": "EU CTIS", "url": "https://euclinicaltrials.eu/search-for-clinical-trials/?lang=en&EUCT=" + i,
    }
    reg_log["CTIS"]["added as new record"] += 1

for r in records.values():
    r["preprint_only"] = bool(r["publications"]) and all(x["preprint"] for x in r["publications"] if not x["protocol"]) and any(not x["protocol"] for x in r["publications"])

# ---------------- output ----------------
out = sorted(records.values(), key=lambda x: (x["population"], x["intervention_class"], x["id"]))
for r in out:
    for kk in ("n_planned", "n_actual"):
        try:
            r[kk] = int(r[kk]) if r[kk] not in (None, "") else None
        except (ValueError, TypeError):
            r[kk] = None
meta = {
    "generated": TODAY.isoformat(), "count": len(out),
    "result_values": ["positive", "null", "negative/harm", "mixed", "feasibility-only", "pending",
                      "unreported-after-completion", "withdrawn", "not-assessable"],
    "result_value_notes": {"feasibility-only": "trial reported only feasibility/acceptability, no efficacy test",
                           "not-assessable": "published or posted but efficacy cannot be read (no abstract, terminated with unusable data, outcomes not reported)",
                           "withdrawn": "registry status WITHDRAWN (never enrolled)"},
    "model_steps_key": {"1": "immune priming", "2": "brain/endocrine control", "3": "low blood volume",
                        "4": "vessels/capillaries", "5": "tissue under-perfusion", "6": "tissue injury/immune response",
                        "7": "sickness response", "8": "sympathetic arousal/poor sleep"},
    "notes": "verified=true only where the result was read in a published abstract or ClinicalTrials.gov posted results. "
             "unreported-after-completion = completed/terminated/unknown status with primary completion before " + CUTOFF + " and no results found.",
    "exclusions_ctgov": dict(excl), "registry_screening": {k: dict(v) for k, v in reg_log.items()},
    "forced_includes_from_publications": forced,
}
json.dump({"meta": meta, "trials": out}, open(base + "/trials.json", "w"), indent=1, ensure_ascii=False)
C = collections.Counter
print("TOTAL", len(out))
for f in ("population", "design", "status", "result", "intervention_class", "registry", "verified"):
    print(f, dict(C(str(r[f]) for r in out).most_common()))
print("unclassified", [(r["id"], r["intervention"]) for r in out if r["intervention"] in ("unclassified", "") or r["intervention_class"] == OT and not r["model_steps"]][:400])
print("REG", reg_log)
print("FORCED", forced)
