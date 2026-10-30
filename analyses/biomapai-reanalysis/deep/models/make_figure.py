"""Figure: case-vs-control AUC per data layer: authors' BioMapAI (visit split vs person split) vs our random forest
and stacked fusion (person split). Points = mean over repeats, bars = 95% person-bootstrap CI."""
import sys; sys.path.insert(0, '/mnt/project-files/biomap/deep/models')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

s1 = pd.read_csv(OUT + 'results/s1_summary.csv')
s2 = pd.read_csv(OUT + 'results/s2_summary_all_short_long.csv'); s2 = s2[s2.task == 'all']
rows = [('Immune (flow)', 'immune', 'immune_rf'), ('Gut species', 'specie', 'specie_rf'), ('Gut KEGG genes', 'kegg', 'kegg_rf'),
        ('Plasma metabolome', 'metabolome', 'metabolome_rf'), ('Routine labs', 'quest', 'quest_rf'),
        ('All layers', 'omics', 'fusion_rf_stack')]
C = {'visit': '#eb6834', 'person': '#2a78d6', 'ours': '#1baf7a'}
fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=150)
fig.patch.set_facecolor('#fcfcfb'); ax.set_facecolor('#fcfcfb')
for i, (lab, l1, l2) in enumerate(rows):
    y = len(rows) - 1 - i
    for split, off, mk in [('visit', 0.22, 's'), ('person', 0.0, 'o')]:
        r = s1[(s1.layer == l1) & (s1.split == split)].iloc[0]
        ax.errorbar(r.auc, y + off, xerr=[[r.auc - r.lo], [r.hi - r.auc]], fmt=mk, ms=6, color=C[split], elinewidth=1.5, capsize=0)
    r = s2[s2.model == l2].iloc[0]
    ax.errorbar(r.auc_visit, y - 0.22, xerr=[[r.auc_visit - r.lo], [r.hi - r.auc_visit]], fmt='D', ms=6, color=C['ours'], elinewidth=1.5, capsize=0)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], color='#0b0b0b')
ax.axvline(0.5, color='#c3c2b7', lw=1, ls='--')
ax.set_xlim(0.45, 0.92); ax.set_xlabel('AUC, ME/CFS vs control (visit level)', color='#52514e')
for s in ['top', 'right']: ax.spines[s].set_visible(False)
for s in ['left', 'bottom']: ax.spines[s].set_color('#c3c2b7')
ax.tick_params(colors='#52514e'); ax.grid(axis='x', color='#e6e5e0', lw=0.6); ax.set_axisbelow(True)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker='s', ls='', color=C['visit'], label="Authors' BioMapAI, visit split (leaky)"),
     Line2D([], [], marker='o', ls='', color=C['person'], label="Authors' BioMapAI, person split"),
     Line2D([], [], marker='D', ls='', color=C['ours'], label='Our random forest / stacked fusion, person split')]
ax.legend(handles=h, loc='upper left', bbox_to_anchor=(0.0, -0.13), fontsize=7.5, frameon=False, ncol=1)
ax.set_title("'All layers': authors = their pre-selected omics.csv (selection leaks); ours = fusion of 5 layers", fontsize=7.5, color='#52514e', loc='left')
plt.tight_layout(); plt.savefig(OUT + 'fig_auc_comparison.png', facecolor=fig.get_facecolor())
print('saved')
