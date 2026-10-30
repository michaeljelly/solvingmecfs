import json, re, sys
from collections import Counter
base = sys.argv[1]
R = json.load(open(f"{base}/work/epmc_flat.json"))
TRIALPT = re.compile(r"randomi[sz]ed controlled trial|clinical trial|controlled clinical trial", re.I)
BADPT = re.compile(r"^(review|review-article|systematic review|meta-analysis|letter|comment|editorial|published erratum|retraction|clinical trial protocol|news|guideline|practice guideline|case reports?|case-report)$", re.I)
INTERV = re.compile(r"randomi[sz]|placebo|open.?label|pilot (study|trial)|single.?arm|were (treated|given|randomi|assigned|allocated)|received (a |an |the |\d|oral|intravenous|treatment)|treatment with|treated with|intervention group|control group|n.?of.?1|case series|before and after|pre.?post|crossover|cross-over|sham|waiting.?list|wait.?list|usual care|supplementation", re.I)
SECOND = re.compile(r"secondary analys|post.?hoc|follow.?up (study|of|analysis)|long.?term follow|mediat(ion|ors)|moderators?|cost.?effective|cost.?utility|economic evaluation|qualitative|process evaluation|re.?analys|reanalys|protocol|study design|rationale and design|statistical analysis plan|baseline characteristics|predictors of|recruitment|feasibility of recruit|participants'? (views|experiences)|interview", re.I)
ME = re.compile(r"chronic fatigue syndrome|myalgic|me/cfs|cfs/me|\bcfs\b|systemic exertion|chronic fatigue", re.I)
LC = re.compile(r"long.?covid|post.?covid|post.?acute (sequelae|covid)|pasc\b|postacute|long.?haul|persistent covid|post.?covid.?19 condition|post covid|postcovid|long-term (effects|sequelae|symptoms) of covid", re.I)
NONLC = re.compile(r"olfact|smell|anosmia|parosmia|hyposmia|taste|gustat|dysgeus|pulmonary fibrosis|lung fibrosis|interstitial lung|vaccin|hair loss|alopecia|tinnitus|dental|oral health|myocarditis|thrombo", re.I)
ACUTE = re.compile(r"hospitali[sz]ed (patients|adults) with covid|acute covid|mild.?to.?moderate covid|outpatients? with covid|early treatment|within \d+ days of|intensive care|icu\b|mechanically ventilated|critically ill|post.?discharge|after (hospital )?discharge|pneumonia", re.I)
TGT = re.compile(r"fatigu|exertion|pem\b|malaise|exercise (capacity|tolerance|intolerance)|walk(ing)? (test|distance)|6mwt|cpet|vo2|cogniti|brain fog|memory|attention|orthostat|pots\b|postural|autonom|dysautonom|heart rate variab|physical function|functional (status|capacity)|quality of life|sf-36|eq-5d|promis|symptom burden|well.?being|recovery|disab", re.I)
out = []
c = Counter()
for r in R:
    t = r["title"]; ab = r["abstract"]; pts = r["pubtypes"]; txt = t + " " + ab
    istrialpt = any(TRIALPT.search(p) for p in pts)
    bad = [p for p in pts if BADPT.match(p.strip())]
    reason = None
    if r["source"] in ("PAT", "ETH", "AGR", "CBA"): reason = "non-article source"
    elif bad and not istrialpt: reason = "pubtype:" + bad[0]
    elif re.search(r"\b(review|meta-analys|systematic|scoping|protocol for|study protocol|a protocol|erratum|correction|retraction|retracted|commentary|guideline|consensus|survey|case report)\b", t, re.I) and not re.search(r"randomi[sz]ed|trial\b", t, re.I): reason = "title: review/protocol/etc"
    elif re.search(r"protocol", t, re.I): reason = "title: protocol"
    elif not ab: reason = "no abstract" if not istrialpt else None
    if reason is None and not istrialpt and not INTERV.search(txt): reason = "no intervention language"
    isme = bool(ME.search(t) or ME.search(ab[:600]))
    islc = bool(LC.search(txt))
    if reason is None and not isme and not islc: reason = "population not ME/CFS or long COVID"
    if reason is None and islc and not isme:
        if NONLC.search(t) and not TGT.search(t): reason = "LC non-target domain (title)"
        elif ACUTE.search(t + " " + ab[:700]) and not re.search(r"long.?covid|post.?covid (syndrome|condition)|pasc", t, re.I): reason = "acute/post-discharge COVID"
        elif not TGT.search(txt): reason = "LC no target outcome"
    if reason is None and SECOND.search(t): reason = "secondary/follow-up/protocol (title)"
    c[reason or "KEEP"] += 1
    r["screen"] = reason or "KEEP"; r["isme"] = isme; r["islc"] = islc; r["trialpt"] = istrialpt
    out.append(r)
json.dump(out, open(f"{base}/work/epmc_screen.json", "w"))
for k, v in c.most_common(): print(v, k)
