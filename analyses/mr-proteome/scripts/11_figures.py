# Forest plot: selected MR estimates in DecodeME GWAS-1, DecodeME infection-onset and MVP-only
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
R='/mnt/project-files/mecfs-research/deep2/genetics_mr/results/'
r=pd.read_csv(R+'mr_all_results.csv')
sel=[('PPP:CLEC3B:P05452:OID30756:v1','CLEC3B protein (cis, UKB-PPP)'),('deCODE:17680_12','EPHB1 protein (cis, deCODE)'),('PPP:PRDX6:P30041:OID21404:v1','PRDX6 protein (cis)'),
 ('PPP:RABGAP1L:Q5R372:OID20447:v1','RABGAP1L protein (cis)'),('PPP:MGMT:P16455:OID20588:v1','MGMT protein (cis)'),('EBVload_Kamitaki2026','EBV DNA load in blood'),
 ('Ferritin','Ferritin'),('TSAT','Transferrin saturation'),('Retic_count','Reticulocyte count'),('IRF','Immature reticulocyte fraction'),('VitD_Revez2020','25(OH)D'),
 ('met_IPA_Chen2023_p5e-6','Indolepropionate'),('met_ILA_Chen2023_p5e-6','Indolelactate'),('met_KYN_Chen2023_p5e-6','Kynurenine'),('met_BHB_Chen2023_p5e-6','3-hydroxybutyrate'),('BMI','BMI')]
outs=[('DME_GWAS1','DecodeME GWAS-1','#2a78d6'),('DME_infection','DecodeME infection onset','#eb6834'),('MVP_only','MVP only (independent)','#1baf7a')]
fig,ax=plt.subplots(figsize=(8,9)); y0=np.arange(len(sel))[::-1]
for k,(o,lab,c) in enumerate(outs):
    for i,(e,nm) in enumerate(sel):
        x=r[(r.exposure==e)&(r.outcome==o)]
        if not len(x): continue
        b,se=x.b.iloc[0],x.se.iloc[0]; yy=y0[i]+0.25-0.25*k
        ax.plot([np.exp(b-1.96*se),np.exp(b+1.96*se)],[yy,yy],color=c,lw=2,solid_capstyle='round')
        ax.plot(np.exp(b),yy,'o',color=c,ms=6,mec='white',mew=1.5,label=lab if i==0 else None)
ax.axvline(1,color='#888',lw=1); ax.set_xscale('log'); ax.set_yticks(y0); ax.set_yticklabels([s[1] for s in sel])
ax.set_xlabel('OR for ME/CFS per SD (or unit) higher exposure, 95% CI'); ax.set_xlim(0.3,3); t=[0.3,0.5,0.7,1,1.5,2,3]; ax.set_xticks(t); ax.set_xticklabels([str(v) for v in t]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.grid(axis='x',color='#eee'); ax.legend(loc='lower right',frameon=False,fontsize=9)
ax.set_title('Selected Mendelian randomisation estimates',loc='left',fontsize=12)
plt.tight_layout(); plt.savefig(R+'mr_forest_selected.png',dpi=150)
