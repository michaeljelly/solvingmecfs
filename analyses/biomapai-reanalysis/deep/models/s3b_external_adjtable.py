"""Step 3 addendum: our models trained on the authors' ADJUSTED metabolome table (metabolome.csv, 730 features, the
same input their pretrained model used) instead of the QC-normalised log data; features matched by COMP_ID (Germain)
or InChiKey (Che); per-cohort z-scoring of log external values; unmatched features dropped (not zero-filled)."""
import sys; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from s3_external import *

C = ext_cohorts(); rows = []
fm = pd.read_csv(B + 'feature_meta/metabolome_feature_meta.csv', index_col=0)
Xt = load_raw('metabolome')
ytr = meta.loc[Xt.index, 'y'].values; gtr = meta.loc[Xt.index, 'person'].values
for key, idcol in [('che_metabolome', 'INCHIKEY'), ('germain_metabolome', 'COMP_ID')]:
    Xe, ye, pe = C[key]
    feats = [f for f in Xt.columns if fm.loc[f, idcol] in set(Xe.columns)]
    Ae = np.log(np.clip(Xe[[fm.loc[f, idcol] for f in feats]].values.astype(float), 1e-9, None))
    Ae = np.where(np.isfinite(Ae), Ae, np.nan); Ae = pd.DataFrame(Ae).fillna(pd.DataFrame(Ae).min()).values
    A1, A2 = zcoh(Xt[feats].values), zcoh(Ae)
    ens = []
    for M in ['enet', 'rf', 'lgbm']:
        gs = fit_model(M, A1, ytr, gtr); p = gs.predict_proba(A2)[:, 1]; ens.append(p)
        r = metrics(ye.values, p, pe.values); r.update(cohort=key, model=M, variant='adjusted_table_cohortz', n_shared=len(feats), internal_cv_auc=gs.best_score_)
        rows.append(r); print(key, M, len(feats), 'AUC %.3f acc %.3f internal %.3f' % (r['auc'], r['acc'], gs.best_score_), flush=True)
    r = metrics(ye.values, np.mean(ens, 0), pe.values); r.update(cohort=key, model='ensemble3', variant='adjusted_table_cohortz', n_shared=len(feats)); rows.append(r)
    print(key, 'ens AUC %.3f acc %.3f' % (r['auc'], r['acc']), flush=True)
pd.DataFrame(rows).to_csv(OUT + 'results/s3b_external_adjtable.csv', index=False)
