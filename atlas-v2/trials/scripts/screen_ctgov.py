import json, re, sys
from collections import Counter
base = sys.argv[1]
A = json.load(open(f"{base}/work/ct_mecfs_flat.json")); B = json.load(open(f"{base}/work/ct_lc_flat.json"))
allr = {}
for r in A + B: allr[r["nct"]] = r
ME = re.compile(r"myalgic|chronic fatigue syndrome|fatigue syndrome, chronic|me/cfs|\bcfs\b|\bseid\b|systemic exertion intolerance|post.?viral fatigue syndrome|\bme\b", re.I)
LC = re.compile(r"post.{0,3}covid|long.?corona|after (recovery|remission)|long.?covid|post.?acute|pasc\b|post.?covid|postcovid|long.?haul|post covid|persistent (covid|symptoms)|sequela|after covid|post-sars|chronic covid|post.?acute infect|covid.?19 survivors|following covid|post.?infectious", re.I)
OUT = {
 "fatigue": r"fatigu|tired|energy|vitality|chalder|\bfss\b|facit|mfi|cis\b|checklist individual strength|pem\b|post.?exertion|exertion|malaise|crash",
 "exercise": r"fitness|cardiorespiratory|exercise|cpet|cardiopulmonary exercise|vo2|peak oxygen|6.?min|six.?min|walk test|walking|physical capacity|functional capacity|aerobic capacity|endurance|step test|sit.?to.?stand|workload|exercise tolerance|muscle (strength|fatig)|handgrip|physical performance",
 "cognition": r"cogniti|brain fog|memory|attention|concentrat|executive|processing speed|moca|neuropsych|\bnih toolbox|trail making|stroop|digit symbol|montreal",
 "autonomic": r"orthostat|pots\b|postural|tachycardia|heart rate variab|\bhrv\b|autonomic|compass|tilt|standing test|lean test|dysautonom|syncope|baroreflex|blood volume|heart rate",
 "global": r"sf.?36|sf.?12|physical function|functional status|function(al)? (impairment|outcome)|quality of life|qol|eq.?5d|promis|global (impression|health|rating)|pgic|cgi|recovery|symptom (burden|severity|score)|disab|whodas|daily activit|activity level|steps|post.?covid functional status|pcfs|return to work|work ability|bell|karnofsky|well.?being|symptom",
}
OUTc = {k: re.compile(v, re.I) for k, v in OUT.items()}
EXC = {
 "smell_taste": r"smell|olfact|anosmia|parosmia|hyposmia|taste|gustat|dysgeus",
 "lung": r"dyspn|breath|pulmonar|lung|respirat|spirom|fvc|dlco|fibrosis|inspiratory|oxygen saturat|cough|ventilat|mmrc",
 "psych": r"ptsd|post.?traumatic|depress|anxiet|insomnia|sleep|psychiatr|mental health|stress|burnout|grief|loneliness|suicid",
 "pain": r"\bpain|headache|migraine|myalgia|arthralg",
 "hair_skin": r"hair|alopecia|skin|tinnitus|hearing|vision|ocular|eye|dental|oral|voice|dysphag|swallow",
 "cardiac": r"myocard|ejection|cardiac function|troponin|echocardi|endothel|arteri|vascular|thrombo|d.?dimer|coagul",
 "metab": r"diabet|glucose|weight|obesity|bmi|lipid",
}
EXCc = {k: re.compile(v, re.I) for k, v in EXC.items()}
VAX = re.compile(r"vaccin|immunogenicity", re.I)
ACUTE = re.compile(r"hospitali[sz]ed (patients|adults) with covid|acute covid|mild.?to.?moderate covid|early treatment|outpatients? with covid|within \d+ days of (symptom|diagnos|positive)|recently diagnosed|newly diagnosed covid|prevent(ion|ing)? (of )?(long|post|pasc)|icu|intensive care|mechanical ventil|critically ill|discharged|post.?discharge|pneumonia|post.?intensive", re.I)
rows = []
for nct, r in allr.items():
    tconds = " ".join([r["title"], r["brief"], r.get("acronym") or ""] + r["conditions"] + r["keywords"])
    isme = bool(ME.search(" ".join([r["title"], r["brief"]] + r["conditions"] + r["keywords"])))
    # guard "me" false positives: require specific
    _t = " ".join([r["title"], r["brief"]] + r["conditions"] + r["keywords"])
    isme = bool(re.search(r"myalgic|chronic fatigue syndrome|fatigue syndrome, chronic|me/cfs|\bcfs\b|\bseid\b|systemic exertion|post.?viral fatigue syndrome", _t, re.I) or re.search(r"\bME\b", _t))
    islc = bool(re.search(r"covid|sars.?cov|pasc|corona", tconds, re.I)) and bool(LC.search(tconds))
    po = " ; ".join(m or "" for m, t in r["primary_outcomes"])
    tgt = [k for k, c in OUTc.items() if c.search(po)]
    tgt_title = [k for k, c in OUTc.items() if c.search(r["title"] + " " + " ".join(r["conditions"])) and k != "global"]
    exc = [k for k, c in EXCc.items() if c.search(po + " " + r["title"] + " " + " ".join(r["conditions"]))]
    crit = r["criteria"][:2500]
    flags = []
    if VAX.search(r["title"] + " " + " ".join(r["interventions"][i][1] or "" for i in range(len(r["interventions"])))): flags.append("vaccine")
    if ACUTE.search(r["title"] + " " + r["summary"][:600] + " " + crit[:1200]): flags.append("acute?")
    if r["purpose"] in ("BASIC_SCIENCE", "DIAGNOSTIC", "SCREENING", "HEALTH_SERVICES_RESEARCH", "PREVENTION"): flags.append("purpose:" + str(r["purpose"]))
    rows.append(dict(nct=nct, isme=isme, islc=islc, tgt=tgt, tgt_title=tgt_title, exc=exc, flags=flags, po=po[:300], title=r["brief"], purpose=r["purpose"], conds=r["conditions"][:5]))
import os
ov = json.load(open(f"{base}/work/screen_overrides.json")) if os.path.exists(f"{base}/work/screen_overrides.json") else {}
PSY = re.compile(r"depress|anxiet|ptsd|post.?traumatic|insomnia|sleep|stress\b|loneliness|hamilton|phq|gad.?7|mood", re.I)
def dom(po):
    d=[k for k,c in EXCc.items() if c.search(po)]
    if "psych" in d and not PSY.search(po): d.remove("psych")
    return d
for x in rows:
    r = allr[x["nct"]]
    po = " ; ".join(m or "" for m, t in r["primary_outcomes"])
    if x["nct"] in ov:
        x["decision"], x["reason"] = ov[x["nct"]]; continue
    if not x["islc"] and not x["isme"]: x["decision"]="out"; x["reason"]="not ME/CFS or long COVID population"; continue
    if x["purpose"] in ("BASIC_SCIENCE","DIAGNOSTIC","SCREENING") and not x["tgt"]:
        x["decision"]="out"; x["reason"]="not a treatment trial (basic science/diagnostic)"; continue
    if x["isme"]: x["decision"]="in"; x["reason"]="ME/CFS"; continue
    if "vaccine" in x["flags"] and not re.search(r"post.?vacc|vaccine.?(injur|induc|related|syndrome)|long.?covid", r["title"], re.I):
        x["decision"]="out"; x["reason"]="vaccine trial"; continue
    if x["tgt"]: x["decision"]="in"; x["reason"]="target outcome: "+",".join(x["tgt"]); continue
    d=dom(po)
    if d: x["decision"]="out"; x["reason"]="non-target primary outcome: "+d[0]; continue
    if x["tgt_title"]: x["decision"]="in"; x["reason"]="target in title/conditions: "+",".join(x["tgt_title"]); continue
    x["decision"]="in"; x["reason"]="global/unspecified outcome in long COVID population"
json.dump(rows, open(f"{base}/work/screen.json", "w"))
c = Counter((x["isme"], x["islc"]) for x in rows); print(c)
