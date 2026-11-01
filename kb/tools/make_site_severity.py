"""Write site/data/cohorts-severity.json: compact severity view of kb/cohorts.json for the data page."""
import json, os
KB = os.path.join(os.path.dirname(__file__), '..')
OUT = os.path.join(KB, '..', 'site', 'data', 'cohorts-severity.json')
F = ['n_mild', 'n_moderate', 'n_severe', 'n_very_severe', 'n_housebound', 'n_bedbound']
rows = []
for c in json.load(open(os.path.join(KB, 'cohorts.json'))):
    if not any(c.get(f) is not None for f in F):
        continue
    rows.append({'id': c['cohort_id'], 'name': c.get('name'),
                 'basis': c.get('severity_basis_type'), 'scale': c.get('severity_scale'),
                 'n_cases': c.get('n_cases'), 'sample_types': c.get('sample_types'),
                 **{f: c.get(f) for f in F}})
tot = {}
for r in rows:
    b = tot.setdefault(r['basis'] or 'unknown', {f: 0 for f in F})
    for f in F:
        b[f] += r[f] or 0
json.dump({'cohorts': rows, 'naive_totals_by_basis': tot}, open(OUT, 'w'), ensure_ascii=False)
print(len(rows), tot)
