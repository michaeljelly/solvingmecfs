"""Shared loaders for the BioMapAI cohort (all visits)."""
import numpy as np, pandas as pd

B = '/mnt/project-files/biomap/BioMapAI/'
IN = B + 'codes/AI/input/'
OUT = '/mnt/project-files/biomap/deep/stats_replication/'

SYMPTOMS = {'NLT.Orthostatic Intolerance': 'orthostatic_NLT', 'PSQI.sleeping_problem': 'sleep_PSQI',
            'RAND36.Energy/Fatigue': 'fatigue_energy', 'DANA.Cognitive Efficiency': 'cognition_DANA',
            'severity': 'overall_severity'}


def load_meta():
    m = pd.read_csv(IN + 'metadata.csv', index_col=0)
    m['person'] = m.sample_id_tp1
    m['case'] = (m.study_ptorhc == 'MECFS').astype(int)
    m['short'] = (m.illness_duration == 'Short').astype(int)
    m['long'] = (m.illness_duration == 'Long').astype(int)
    m['male'] = (m.gender == 'Male').astype(int)
    m['tp'] = m.timepoints.str[2:].astype(int)
    m['tp2'] = (m.tp == 2).astype(int); m['tp3'] = (m.tp >= 3).astype(int)
    # BMI: constant per person; 39/249 people missing -> sex-specific median imputation + missing indicator
    m['bmi_missing'] = m.bmi.isna().astype(int)
    med = m.drop_duplicates('person').groupby('gender').bmi.median()
    m['bmi_imp'] = m.bmi.fillna(m.gender.map(med))
    m['age_c'] = (m.age - 45) / 10
    m['bmi_c'] = (m.bmi_imp - 26) / 5
    s = pd.read_csv(IN + 'score.csv', index_col=0)
    s['RAND36.Role Limitations due to Physical Health'] /= 4
    s['RAND36.Role Limitations due to Emotional Problems'] /= 3
    s['severity'] = s.mean(axis=1)   # mean of 12 scaled scores, higher = worse
    return m.join(s, how='left')


def _int(X):
    """rank-based inverse normal transform per column"""
    from scipy.stats import norm
    R = X.rank(axis=0)
    return pd.DataFrame(norm.ppf((R - 0.5) / len(X)), index=X.index, columns=X.columns)


def clr(rel, prev=0.1, det=1e-6):
    """rel: samples x taxa relative abundance. Prevalence filter then CLR with pseudocount = half min nonzero."""
    rel = rel.div(rel.sum(axis=1), axis=0)
    keep = (rel > det).mean() >= prev
    r = rel.loc[:, keep]
    pc = r[r > 0].min().min() / 2
    lg = np.log(r + pc)
    return lg.sub(lg.mean(axis=1), axis=0)


def load_layer(L):
    if L == 'immune':      # raw flow percentages (the input table is a covariate residual)
        X = pd.read_csv(B + 'data/immuno/immune_percentage.csv', index_col=0).T
        X = X.loc[:, X.std() > 0]
        X = _int(X)
    elif L == 'species':   # full MetaPhlAn4 table, >=10% prevalence, CLR
        d = pd.read_csv(B + 'data/metagenomics/metaphlan4_species_full.csv', index_col=0)
        d = d.select_dtypes('number').T
        X = clr(d)
    elif L == 'kegg':      # gene counts-like; total-sum to CPM, log2, >=10% prevalence
        d = pd.read_csv(IN + 'kegg.csv', index_col=0).T
        d = d.loc[:, (d > 0).mean() >= 0.1]
        X = np.log2(d.div(d.sum(axis=1), axis=0) * 1e6 + 1)
        X = X.loc[:, X.std() > 0]
    elif L == 'metabolome':  # QC-normalised imputed Metabolon data, named compounds, log2
        q = pd.read_csv(B + 'data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv')
        q = q[q.TYPE == 'NAMED']
        q.index = q.CHEMICAL_NAME
        S = [c for c in q.columns if c.startswith('SAM')]
        X = np.log2(q[S].astype(float)).T
        # drop compounds detected (above the imputed minimum) in <20% of samples (mostly drugs/rare xenobiotics)
        X = X.loc[:, (X > X.min() + 1e-9).mean() >= 0.2]
        X = X.loc[:, X.std() > 0]
    elif L == 'quest':     # routine labs as provided (already residualised by the authors)
        X = pd.read_csv(IN + 'quest.csv', index_col=0).T
        X = _int(X)
    return X


def metab_meta():
    q = pd.read_csv(B + 'data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv')
    q = q[q.TYPE == 'NAMED']
    return q.set_index('CHEMICAL_NAME')[['SUPER_PATHWAY', 'SUB_PATHWAY', 'INCHIKEY', 'COMP_ID', 'HMDB', 'KEGG']]


def annotate(L, feats):
    if L == 'metabolome':
        mm = metab_meta()
        return [f"{mm.SUPER_PATHWAY.get(f, '')} / {mm.SUB_PATHWAY.get(f, '')}" for f in feats]
    if L == 'kegg':
        km = pd.read_csv(IN + 'feature_meta/kegg_feature_meta.csv', index_col=0).definition
        return [str(km.get(f, ''))[:70] for f in feats]
    return [''] * len(feats)


LAYERS = ['immune', 'species', 'kegg', 'metabolome', 'quest']
