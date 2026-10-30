import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R='results/'
t=pd.read_csv(R+'plasma_targets_G2020_G2022_meta.csv',index_col=0)
ext=pd.read_csv('/mnt/project-files/biomap/deep/coupling/external_replication.csv')
che=ext[(ext.cohort=='Che2022')&ext.metabolite.str.startswith('indolepropionate [Indole')].iloc[0]
rows=[('BioMapAI (US, 147/89)',-0.355,-0.597,-0.113),('Che 2022 (US, 79/91)',che.beta_SD,che.lo,che.hi),
      ('Germain 2022 (US, 60/45)',t.loc['indolepropionate','G2022_est'],t.loc['indolepropionate','G2022_est']-1.96*t.loc['indolepropionate','G2022_se'],t.loc['indolepropionate','G2022_est']+1.96*t.loc['indolepropionate','G2022_se']),
      ('Germain 2020 (US, 26/26)*',t.loc['indolepropionate','G2020_est'],t.loc['indolepropionate','G2020_est']-1.96*t.loc['indolepropionate','G2020_se'],t.loc['indolepropionate','G2020_est']+1.96*t.loc['indolepropionate','G2020_se']),
      ('Pooled (random effects)',t.loc['indolepropionate','meta_RE'],t.loc['indolepropionate','meta_lo'],t.loc['indolepropionate','meta_hi'])]
g=pd.read_csv(R+'gut_nih_prespecified.csv')
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for i,(n,e,lo,hi) in enumerate(rows[::-1]):
    ax[0].plot([lo,hi],[i,i],color='#444'); ax[0].plot(e,i,'s' if 'Pooled' in n else 'o',color='#b2182b' if 'Pooled' in n else '#2166ac',ms=8)
ax[0].set_yticks(range(len(rows))); ax[0].set_yticklabels([r[0] for r in rows[::-1]]); ax[0].axvline(0,color='grey',ls=':')
ax[0].set_xlabel('Patients minus controls (SD units)'); ax[0].set_title('Plasma indole-3-propionate (IPA)')
g=g.sort_values('NIH_est'); cols=['#2166ac' if s else '#999' for s in g.same_dir]
ax[1].barh(range(len(g)),g.NIH_est,color=cols); ax[1].set_yticks(range(len(g))); ax[1].set_yticklabels(g.feature.str.replace('_',' ').str.replace('COMPOSITE ','all '),fontsize=8)
ax[1].axvline(0,color='grey',ls=':'); ax[1].set_xlabel('NIH cohort: patients minus controls (SD, CLR)'); ax[1].set_title('NIH stool (15 ME/CFS, 17 HV)\nblue = same direction as before (13/15)')
plt.figtext(0.01,0.01,'*Germain 2020 labels corrected from drug traces; may share participants with Germain 2022.',fontsize=7)
plt.tight_layout(rect=[0,0.03,1,1]); plt.savefig(R+'fig_ipa_gut_replication.png',dpi=150)
