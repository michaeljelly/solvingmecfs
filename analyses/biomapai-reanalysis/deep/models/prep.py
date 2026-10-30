"""Feature preparation and model pipelines (all learned steps live inside sklearn pipelines -> fit within folds)."""
import sys; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from common import *
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import SelectKBest, f_classif, f_regression
from sklearn.linear_model import LogisticRegression, ElasticNet
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
import lightgbm as lgb
import warnings; warnings.filterwarnings('ignore')


class CompFilter(BaseEstimator, TransformerMixin):
    """For compositional layers: per-sample CLR (or log) transform, then keep features whose prevalence (raw>0)
    in the TRAINING data is >= min_prev. mode: 'clr' | 'log' | 'none'."""
    def __init__(self, mode='clr', min_prev=0.1, pseudo=1e-6):
        self.mode = mode; self.min_prev = min_prev; self.pseudo = pseudo

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        self.keep_ = (X > 0).mean(0) >= self.min_prev
        if self.keep_.sum() < 2: self.keep_[:] = True
        return self

    def transform(self, X):
        X = np.asarray(X, float)
        if self.mode == 'none': return X[:, self.keep_]
        L = np.log10(X + self.pseudo)
        if self.mode == 'clr': L = L - L.mean(1, keepdims=True)
        return L[:, self.keep_]


def load_layer(L):
    """Returns DataFrame samples x features ready for pipelines (compositional layers stay raw here;
    species as relative abundance, KEGG converted to per-sample relative abundance)."""
    if L == 'demo':
        d = meta[['age', 'male', 'bmi']].copy(); return d
    if L == 'metabqc':  # less-processed Metabolon QC-normalised, imputed data (957 metabolites), natural log
        q = pd.read_csv('/mnt/project-files/biomap/BioMapAI/data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv', index_col=0)
        X = np.log(q.iloc[:, 17:].T.astype(float)); X.columns = q.CHEMICAL_NAME.values
        return X.loc[X.index.intersection(meta.index)]
    if L == 'metabqc_noxeno':  # same, excluding xenobiotics and unnamed (X-) compounds
        q = pd.read_csv('/mnt/project-files/biomap/BioMapAI/data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv', index_col=0)
        q = q[(q.SUPER_PATHWAY != 'Xenobiotics') & q.SUPER_PATHWAY.notna()]
        X = np.log(q.iloc[:, 17:].T.astype(float)); X.columns = q.CHEMICAL_NAME.values
        return X.loc[X.index.intersection(meta.index)]
    X = load_raw(L)
    if L == 'kegg': X = X.div(X.sum(1), axis=0)
    return X


COMP = {'specie': 1e-5, 'kegg': 3e-7}  # ~1st percentile of non-zero relative abundances


def make_pipe(L, model, task='clf'):
    steps = []
    if L in COMP: steps.append(('comp', CompFilter('clr', 0.1, COMP[L])))
    steps.append(('imp', SimpleImputer(strategy='median')))
    if model == 'enet':
        steps.append(('sc', StandardScaler()))
        if L == 'kegg': steps.append(('sel', SelectKBest(f_classif if task == 'clf' else f_regression, k=500)))
        if task == 'clf':
            steps.append(('m', LogisticRegression(penalty='elasticnet', solver='saga', C=0.1, l1_ratio=0.5,
                                                  max_iter=3000, tol=1e-3, class_weight='balanced')))
        else:
            steps.append(('m', ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=5000)))
    elif model == 'rf':
        if task == 'clf':
            steps.append(('m', RandomForestClassifier(300, max_features='sqrt', min_samples_leaf=3, class_weight='balanced',
                                                      n_jobs=1, random_state=0)))
        else:
            steps.append(('m', RandomForestRegressor(300, max_features=0.3, min_samples_leaf=5, n_jobs=1, random_state=0)))
    elif model == 'lgbm':
        kw = dict(n_estimators=200, learning_rate=0.03, num_leaves=7, min_child_samples=10, subsample=0.8,
                  subsample_freq=1, colsample_bytree=0.3, reg_lambda=1.0, n_jobs=1, verbose=-1, random_state=0)
        steps.append(('m', lgb.LGBMClassifier(class_weight='balanced', **kw) if task == 'clf' else lgb.LGBMRegressor(**kw)))
    return Pipeline(steps)


def grid(L, model, task='clf'):
    if model == 'enet':
        if task == 'clf':
            g = {'m__C': [0.01, 0.1, 1.0], 'm__l1_ratio': [0.1, 0.7]}
        else:
            g = {'m__alpha': [0.03, 0.1, 0.3, 1.0], 'm__l1_ratio': [0.1, 0.7]}
        if L == 'kegg':
            g = dict(g); g['m__l1_ratio'] = [0.5]; g['sel__k'] = [100, 500]
        if L in COMP: g = dict(g); g['comp__mode'] = ['clr', 'log']
        if L == 'demo': g = {'m__C': [0.1, 1.0], 'm__l1_ratio': [0.5]} if task == 'clf' else {'m__alpha': [0.1], 'm__l1_ratio': [0.5]}
        return g
    if model == 'rf':
        return {'m__min_samples_leaf': [1, 5], 'm__max_features': ['sqrt', 0.2]} if task == 'clf' else \
               {'m__min_samples_leaf': [5, 15], 'm__max_features': [0.1, 0.3]}
    if model == 'lgbm':
        return {'m__n_estimators': [100, 300], 'm__num_leaves': [4, 15]}
