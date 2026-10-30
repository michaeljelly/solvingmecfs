"""Figures: (1) CSF metabolome global shift + run-order; (2) tilt NE/EPI trajectories; (3) key reproduced group differences."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
B = '/mnt/project-files/mecfs-research/downloads/nih_walitt'; T = f'{B}/results/tidy'; F = f'{B}/figures'
HV, PI = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#52514e', 'axes.labelcolor': '#0b0b0b'})

m = pd.read_csv(f'{T}/csf_metabolomics_wide.csv'); L = np.log2(m.iloc[:, 4:].astype(float))
hv = m.group == 'HV'; fc = L[~hv].mean() - L[hv].mean(); med = L.median(axis=1)
num = m.sample_id.str.extract(r'(\d+)')[0].astype(int)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].hist(fc, bins=40, color='#8a8984', edgecolor='white', linewidth=0.8)
ax[0].axvline(0, color='#0b0b0b', lw=1); ax[0].axvline(fc.median(), color=PI, lw=2)
ax[0].set_xlabel('log2 fold change, PI-ME/CFS vs HV (445 CSF metabolites)'); ax[0].set_ylabel('metabolites')
ax[0].set_title('94%% of CSF metabolites lower in patients\nmedian shift %.2f log2 (%.0f%%)' % (fc.median(), 100 * (2 ** fc.median() - 1)), fontsize=9, loc='left')
for g, c, mk in [('HV', HV, 'o'), ('PI-ME/CFS', PI, 's')]:
    s = m.group == g; ax[1].scatter(num[s], med[s], color=c, s=36, marker=mk, edgecolor='white', linewidth=1, label=g, zorder=3)
ax[1].axhline(0, color='#c3c2b7', lw=1, zorder=1)
ax[1].set_xlabel('Metabolon sample number (NINR-00xxx)'); ax[1].set_ylabel('per-sample median log2 level')
ax[1].set_title('Groups were submitted in blocks; level falls\nwith sample number inside the patient block', fontsize=9, loc='left'); ax[1].legend(frameon=False)
fig.tight_layout(); fig.savefig(f'{F}/fig1_csf_metabolome_global_shift.png', dpi=150)

t = pd.read_csv(f'{T}/tilt_plasma_catecholamines_wide.csv')
mins = ['SUPINE', '4', '8', '12', '16', '20', '24', '28', '32', '36', '40']; x = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
for k, (pre, lab) in enumerate([('NE', 'plasma norepinephrine (nmol/L)'), ('EPI', 'plasma epinephrine (nmol/L)')]):
    for g, c in [('HV', HV), ('PI-ME/CFS', PI)]:
        d = t[t.group == g][[f'{pre}.{mm}' for mm in mins]].values
        for r in d: ax[k].plot(x, r, color=c, lw=0.6, alpha=0.3)
        ax[k].plot(x, np.nanmedian(d, axis=0), color=c, lw=2.5, label=f'{g} median (n={len(d)})')
    ax[k].set_xlabel('minutes of head-up tilt (0 = supine)'); ax[k].set_ylabel(lab); ax[k].legend(frameon=False, fontsize=8)
ax[1].set_yscale('log')
ax[0].set_title('Supine NE lower in patients (MWU p=0.026);\nupright NE similar', fontsize=9, loc='left')
ax[1].set_title('EPI tends higher in patients throughout\n(upright mean perm p=0.027, not FDR-significant)', fontsize=9, loc='left')
fig.tight_layout(); fig.savefig(f'{F}/fig2_tilt_catecholamines.png', dpi=150)

A = pd.read_csv(f'{B}/results/A_scalar_group_tests.csv').sort_values('hedges_g')
A = A[A.p_mwu < 0.1]
fig, ax = plt.subplots(figsize=(8, 5.2))
y = np.arange(len(A)); cols = ['#0b0b0b' if q < 0.05 else '#8a8984' for q in A.q_BH_mwu]
ax.hlines(y, 0, A.hedges_g, color='#c3c2b7', lw=2); ax.scatter(A.hedges_g, y, color=cols, s=40, zorder=3)
ax.set_yticks(y); ax.set_yticklabels([f'{m} (n={a}/{b})' for m, a, b in zip(A.measure, A.n_HV, A.n_PI)], fontsize=7.5)
ax.axvline(0, color='#52514e', lw=1); ax.set_xlabel("Hedges' g, PI-ME/CFS minus HV")
fig.suptitle('Person-level measures with MWU p<0.10\n(black: BH q<0.05 across 42 measures)', fontsize=9, x=0.02, ha='left')
fig.subplots_adjust(left=0.32, top=0.88, bottom=0.1, right=0.97); fig.savefig(f'{F}/fig3_effect_sizes.png', dpi=150)
print('ok')
