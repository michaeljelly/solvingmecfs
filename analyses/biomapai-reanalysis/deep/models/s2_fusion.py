"""Step 2b: late fusion across layers (copes with missing layers) + summary table with person-bootstrap CIs.
Fusion A (mean): average of the logits of whichever layer predictions a sample has.
Fusion B (stack): logistic regression on layer logits (missing -> 0) + missing indicators, trained on the inner
out-of-fold predictions of the outer-training persons only.
Evaluation: visit-level AUC and person-level AUC (mean of a person's visit predictions), predictions averaged over
3 repeats, 95% CI by bootstrapping persons."""
import sys, glob; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from common import *
from sklearn.linear_model import LogisticRegression
from s2_classify import task_samples

TASKS = sys.argv[1].split(',') if len(sys.argv) > 1 else ['all', 'short', 'long']
logit = lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) / (1 - np.clip(p, 1e-4, 1 - 1e-4)))


def evaluate(d, name, task, rows):
    """d: columns sample, person, y, rep, pred (held-out predictions). AUC = mean over repeats, CI by person bootstrap."""
    d = d.dropna(subset=['pred'])
    a, lo, hi = boot_auc_reps(d)
    pa, plo, phi = boot_auc_reps(d, person_level_agg=True)
    bacc = np.mean([((x.pred > .5) == x.y).groupby(x.y).mean().mean() for _, x in d.groupby('rep')])
    rows.append(dict(task=task, model=name, n_visits=d['sample'].nunique(), n_persons=d.person.nunique(), auc_visit=a, lo=lo, hi=hi,
                     auc_person=pa, p_lo=plo, p_hi=phi, bal_acc_visit=bacc))


if __name__ == '__main__':
    rows = []
    for task in TASKS:
        files = sorted(glob.glob(OUT + f'results/s2/pred_{task}_*.csv'))
        P = {}
        for f in files:
            key = f.split(f'pred_{task}_')[1][:-4]  # layer_model
            P[key] = pd.read_csv(f)
            evaluate(P[key][P[key].set == 'test'], key, task, rows)
        m = task_samples(task)
        omics_keys = [k for k in P if not k.startswith('demo') and not k.startswith('metabqc')]
        combos = {}
        for fam in ['enet', 'rf', 'lgbm']:
            ks = [k for k in omics_keys if k.endswith('_' + fam)]
            if len(ks) >= 2: combos['fusion_' + fam] = ks
        if len(omics_keys) >= 2: combos['fusion_allfamilies'] = omics_keys
        demo = [k for k in P if k.startswith('demo')]
        # variant: swap the adjusted metabolome table for the less-processed QC-normalised one
        for fam in ['enet', 'lgbm']:
            if f'metabqc_{fam}' in P:
                combos[f'fusion_{fam}_metabqc'] = [k for k in combos.get('fusion_' + fam, []) if not k.startswith('metabolome')] + [f'metabqc_{fam}']
        if demo and 'fusion_allfamilies' in combos: combos['fusion_allfamilies_plus_demo'] = combos['fusion_allfamilies'] + demo
        for name, ks in list(combos.items()):
            for mode in ['mean', 'stack']:
                out = []
                for rep in range(3):
                    for k in range(5):
                        def mat(setname):
                            cols = {}
                            for key in ks:
                                x = P[key]; x = x[(x.rep == rep) & (x.fold == k) & (x.set == setname)]
                                cols[key] = pd.Series(logit(x.pred.values), index=x['sample'].values)
                            return pd.DataFrame(cols)
                        Te = mat('test'); Tr = mat('trainoof')
                        if mode == 'mean':
                            p = Te.mean(1)
                        else:
                            def feats(M): return np.hstack([M.fillna(0).values, M.isna().values.astype(float)])
                            ytr = m.loc[Tr.index, 'y'].values
                            lr = LogisticRegression(C=0.3, max_iter=2000).fit(feats(Tr), ytr)
                            p = pd.Series(lr.predict_proba(feats(Te))[:, 1], index=Te.index)
                        out.append(pd.DataFrame({'sample': Te.index, 'rep': rep, 'pred': p.values}))
                o = pd.concat(out); o['person'] = m.loc[o['sample'], 'person'].values; o['y'] = m.loc[o['sample'], 'y'].values
                o.to_csv(OUT + f'results/s2/fused_{task}_{name}_{mode}.csv', index=False)
                evaluate(o, f'{name}_{mode}', task, rows)
                print(task, name, mode, round(rows[-1]['auc_visit'], 3), round(rows[-1]['auc_person'], 3), flush=True)
    r = pd.DataFrame(rows)
    r.to_csv(OUT + 'results/s2_summary' + ('' if len(sys.argv) < 2 else '_' + sys.argv[1].replace(',', '_')) + '.csv', index=False)
    print(r[['task', 'model', 'n_visits', 'auc_visit', 'lo', 'hi', 'auc_person', 'p_lo', 'p_hi']].round(3).to_string())
