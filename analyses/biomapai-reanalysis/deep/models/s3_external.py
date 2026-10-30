"""Step 3: external validation. Train on the full BioMapAI cohort, test on the authors' bundled independent cohorts.
(a) Authors' pretrained BioMapAI models applied exactly as in Validation_Independent_cohorts.ipynb (reindex to model
    features, fill 0; metabolome StandardScaler per cohort; their predict() re-standardises per cohort).
(b) Our models (elastic-net LR, RF, LightGBM, tuned by person-grouped CV on BioMapAI) trained only on the features
    shared with the target cohort; microbiome CLR; metabolome log; variant with per-cohort feature z-scoring."""
import sys, json; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from prep import *
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score

IC = '/mnt/project-files/biomap/BioMapAI/data/indenpendent_cohorts/'
os.makedirs(OUT + 'results/s3', exist_ok=True)


# ---------------- external cohorts ----------------
def ext_cohorts():
    C = {}
    sra = pd.read_csv(IC + 'columbia/SraRunTable.txt', index_col=0)
    d = pd.read_csv(IC + 'columbia/metaphlan4_species_full.csv', index_col=0).T
    C['columbia_specie'] = (d, sra.loc[d.index, 'sampletype'].map({'Control': 0, 'Case': 1}), pd.Series(d.index, d.index))
    d = pd.read_table(IC + 'columbia/kegg_gene_normalized.txt', index_col=0).drop(columns='decode').T.astype(float)
    C['columbia_kegg'] = (d, sra.loc[d.index, 'sampletype'].map({'Control': 0, 'Case': 1}), pd.Series(d.index, d.index))
    d = pd.read_csv(IC + 'Raijmakers/metaphlan4_species_full.csv', index_col=0).T
    d = d.loc[[i for i in d.index if 'CFS' in i or 'HV' in i]]
    C['raijmakers_specie'] = (d, pd.Series([0 if 'HV' in i else 1 for i in d.index], d.index), pd.Series(d.index, d.index))
    d = pd.read_table(IC + 'Raijmakers/kegg_gene_normalized.txt', index_col=0).T.astype(float)
    d.index = [i.split('.')[0] for i in d.index]; d = d.loc[[i for i in d.index if 'CFS' in i or 'HV' in i]]
    C['raijmakers_kegg'] = (d, pd.Series([0 if 'HV' in i else 1 for i in d.index], d.index), pd.Series(d.index, d.index))
    # Che 2022 plasma metabolome (raw intensities) - features keyed by InChiKey
    raw = pd.read_csv(IC + 'Int.J.Mol.Sci.2022/r56_metabolomics_yali_v2.csv', index_col=0)
    ph = pd.read_table(IC + 'Int.J.Mol.Sci.2022/r56_pheno_text_binned_v2.txt', index_col=0)
    X = raw.iloc[:, 3:].apply(pd.to_numeric, errors='coerce').fillna(0).T
    X.columns = raw.InChiKey.values
    X = X.T.groupby(level=0).mean().T  # duplicate keys -> mean
    C['che_metabolome'] = (X, ph.loc[X.index, 'Phenotype'].map({'control': 0, 'case': 1}), pd.Series(X.index, X.index))
    # Germain 2022 (jci insight) plasma metabolome, 4 time points per person, keyed by COMP ID
    g = pd.read_excel(IC + 'jci.insight.2022.xlsx', sheet_name='ScaledImpDataZeroDrug&Tobacco', header=None)
    hdr_row = g.index[g.iloc[:, 0].astype(str).str.contains('PATHWAY SORTORDER', na=False)][0]
    fcols = g.iloc[hdr_row, :15].tolist(); body = g.iloc[hdr_row + 1:]
    body = body[pd.to_numeric(body.iloc[:, fcols.index('COMP ID')], errors='coerce').notna()]
    comp = body.iloc[:, fcols.index('COMP ID')].astype(int).values
    ids = g.iloc[0, 16:].values; ph = g.iloc[1, 16:].values; tp = g.iloc[2, 16:].values
    vals = body.iloc[:, 16:].astype(float).values.T
    sid = [f'{a}|{b}' for a, b in zip(ids, tp)]
    X = pd.DataFrame(vals, index=sid, columns=comp)
    y = pd.Series([1 if p == 'CFS' else 0 for p in ph], index=sid)
    C['germain_metabolome'] = (X, y, pd.Series(ids, index=sid))
    C['germain_tp'] = pd.Series(tp, index=sid)
    return C


# ---------------- BioMapAI training data ----------------
def biomap_train(layer):
    if layer in ('specie', 'kegg'):
        X = load_raw(layer)
    else:  # metabolome: less-processed QC-normalised imputed Metabolon data (so a log transform is comparable)
        q = pd.read_csv('/mnt/project-files/biomap/BioMapAI/data/metabolomics/QC-norm Data Common wImp_Jan_06_2022.csv', index_col=0)
        X = q.iloc[:, 17:].T.astype(float)
        X = X.loc[X.index.intersection(meta.index)]
        X.attrs['comp'] = q.COMP_ID.values; X.attrs['inchi'] = q.INCHIKEY.values
    return X


def transform(X, layer, pseudo):
    X = np.asarray(X, float)
    if layer in ('specie', 'kegg'):
        X = X / np.maximum(X.sum(1, keepdims=True), 1e-12)
        L = np.log10(X + pseudo); return L - L.mean(1, keepdims=True)
    Xp = X.copy(); pos = Xp[Xp > 0]; Xp[Xp <= 0] = pos.min() / 2 if pos.size else 1e-6
    return np.log(Xp)


def zcoh(A):
    A = np.asarray(A, float); sd = A.std(0); sd[sd == 0] = 1
    return (A - A.mean(0)) / sd


def fit_model(M, Xtr, ytr, gtr):
    steps = [('sc', StandardScaler())]
    if M == 'enet':
        steps.append(('m', LogisticRegression(penalty='elasticnet', solver='saga', max_iter=3000, tol=1e-3,
                                              class_weight='balanced')))
        gr = {'m__C': [0.01, 0.1, 1.0], 'm__l1_ratio': [0.1, 0.7]}
    elif M == 'rf':
        steps.append(('m', RandomForestClassifier(300, class_weight='balanced', n_jobs=1, random_state=0)))
        gr = {'m__min_samples_leaf': [1, 5], 'm__max_features': ['sqrt', 0.2]}
    else:
        steps.append(('m', lgb.LGBMClassifier(learning_rate=0.03, min_child_samples=10, subsample=0.8, subsample_freq=1,
                                              colsample_bytree=0.3, reg_lambda=1.0, class_weight='balanced', n_jobs=1,
                                              verbose=-1, random_state=0)))
        gr = {'m__n_estimators': [100, 300], 'm__num_leaves': [4, 15]}
    gs = GridSearchCV(Pipeline(steps), gr, cv=StratifiedGroupKFold(5, shuffle=True, random_state=0), scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr, ytr, groups=gtr)
    return gs


def metrics(y, p, persons=None):
    y = np.asarray(y); p = np.asarray(p)
    a, lo, hi = boot_auc(y, p, persons if persons is not None else np.arange(len(y)))
    return dict(auc=a, lo=lo, hi=hi, acc=accuracy_score(y, p > .5), bal_acc=balanced_accuracy_score(y, p > .5), n=len(y))


if __name__ == '__main__':
    C = ext_cohorts(); rows = []
    # (a) authors' pretrained models
    sys.path.insert(0, OUT); from biomapai_dnn import DNN, tf
    fmeta = {L: pd.read_csv(B + f'feature_meta/{L}_feature_meta.csv', index_col=0) for L in ['specie', 'kegg', 'metabolome']}
    for setname, path, sm, ym in [('DeepMECFS', '/mnt/project-files/biomap/pretrained_model_DeepMECFS/', 'DeepMECFS_', 'Y2y_'),
                                  ('codes_AI_model', '/mnt/project-files/biomap/BioMapAI/codes/AI/model/', 'DNN_score_', 'Adjust_score_y_')]:
        for key in ['columbia_specie', 'raijmakers_specie', 'columbia_kegg', 'raijmakers_kegg', 'che_metabolome', 'germain_metabolome']:
            L = key.split('_')[1]; X, y, pers = C[key]
            om = tf.keras.models.load_model(path + sm + L); yl = tf.keras.models.load_model(path + ym + L)
            fm = fmeta[L]
            if key == 'che_metabolome':
                Xa = pd.DataFrame(0.0, index=X.index, columns=fm.index)
                for f in fm.index:
                    if fm.INCHIKEY[f] in X.columns: Xa[f] = X[fm.INCHIKEY[f]].values
                Xa = StandardScaler().fit_transform(Xa)
            elif key == 'germain_metabolome':
                Xa = X.reindex(columns=fm.COMP_ID.values, fill_value=0).values
                Xa = StandardScaler().fit_transform(Xa)
            else:
                Xa = X.reindex(columns=fm.index, fill_value=0).values
            p = DNN.ScoreYModel(om, yl).predict(np.asarray(Xa, 'float32')).ravel()
            r = metrics(y, p, pers.values); r.update(cohort=key, model=f'authors_pretrained_{setname}', variant='as_notebook')
            rows.append(r); print(key, setname, {k: round(v, 3) for k, v in r.items() if k in ('auc', 'acc')}, flush=True)
            if key == 'germain_metabolome':
                pp = pd.DataFrame({'p': p, 'y': y.values, 'person': pers.values}).groupby('person').mean()
                r = metrics(pp.y, pp.p); r.update(cohort='germain_metabolome_personmean', model=f'authors_pretrained_{setname}', variant='as_notebook'); rows.append(r)
            tf.keras.backend.clear_session()
    pd.DataFrame(rows).to_csv(OUT + 'results/s3_external.csv', index=False)

    # (b) our models
    tr = {L: biomap_train(L) for L in ['specie', 'kegg', 'metabolome']}
    for key in ['columbia_specie', 'raijmakers_specie', 'columbia_kegg', 'raijmakers_kegg', 'che_metabolome', 'germain_metabolome']:
        L = key.split('_')[1]; Xe, ye, pe = C[key]; Xt = tr[L]
        if L == 'metabolome':
            ids = Xt.attrs['inchi'] if key == 'che_metabolome' else Xt.attrs['comp']
            Xt = Xt.copy(); Xt.columns = ids; Xt = Xt.loc[:, ~pd.isna(Xt.columns)]; Xt = Xt.T.groupby(level=0).mean().T
            if key == 'che_metabolome':  # exact InChiKey, then fall back to the first (connectivity) block
                ext_cols = {}
                for c in Xt.columns:
                    if c in Xe.columns: ext_cols[c] = c
                    else:
                        cand = [e for e in Xe.columns if isinstance(e, str) and e[:14] == str(c)[:14]]
                        if cand: ext_cols[c] = cand[0]
                shared = list(ext_cols); Xe2 = pd.DataFrame({c: Xe[ext_cols[c]].values for c in shared}, index=Xe.index)
            else:
                shared = [c for c in Xt.columns if c in Xe.columns]; Xe2 = Xe[shared]
            n_exact = len([c for c in fmeta['metabolome'].INCHIKEY if c in Xe.columns]) if key == 'che_metabolome' else None
        else:
            shared = [c for c in Xt.columns if c in Xe.columns]; Xe2 = Xe[shared]
        Xt2 = Xt[shared]
        prev = (Xt2 > 0).mean(0).values >= 0.1 if L != 'metabolome' else np.ones(len(shared), bool)
        pseudo = COMP.get(L, 1e-6)
        At = transform(Xt2.values, L, pseudo)[:, prev]; Ae = transform(Xe2.values, L, pseudo)[:, prev]
        ytr = meta.loc[Xt2.index, 'y'].values; gtr = meta.loc[Xt2.index, 'person'].values
        info = dict(n_shared=len(shared), n_used=int(prev.sum()), n_train_features=Xt.shape[1])
        print(key, info, flush=True)
        variants = ['cohortz'] if L == 'metabolome' else ['raw', 'cohortz']
        for v in variants:
            A1, A2 = (At, Ae) if v == 'raw' else (zcoh(At), zcoh(Ae))
            ens = []
            for M in ['enet', 'rf', 'lgbm']:
                gs = fit_model(M, A1, ytr, gtr)
                p = gs.predict_proba(A2)[:, 1]; ens.append(p)
                r = metrics(ye.values, p, pe.values); r.update(cohort=key, model=M, variant=v, internal_cv_auc=gs.best_score_, **info)
                rows.append(r); print(' ', v, M, 'AUC %.3f acc %.3f (internal %.3f)' % (r['auc'], r['acc'], gs.best_score_), flush=True)
                if key == 'germain_metabolome':
                    for sub, mask in [('personmean', None), ('D1PRE', (C['germain_tp'] == 'D1-PRE').values)]:
                        if mask is None:
                            pp = pd.DataFrame({'p': p, 'y': ye.values, 'person': pe.values}).groupby('person').mean(); r = metrics(pp.y, pp.p)
                        else:
                            r = metrics(ye.values[mask], p[mask])
                        r.update(cohort=f'germain_metabolome_{sub}', model=M, variant=v, **info); rows.append(r)
            p = np.mean(ens, 0); r = metrics(ye.values, p, pe.values); r.update(cohort=key, model='ensemble3', variant=v, **info); rows.append(r)
            print(' ', v, 'ens', 'AUC %.3f acc %.3f' % (r['auc'], r['acc']), flush=True)
            pd.DataFrame(rows).to_csv(OUT + 'results/s3_external.csv', index=False)
