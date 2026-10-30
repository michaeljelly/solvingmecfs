from sub_common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
C1,C2,C3,GR,MU,INK='#2a78d6','#eb6834','#1baf7a','#898781','#c9c8c2','#0b0b0b'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':GR,'axes.labelcolor':INK,'xtick.color':'#52514e','ytick.color':'#52514e','figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb','legend.frameon':False})
F=OUT+'fig/'
# Fig 1: symptom clustering vs copula null + PCA map
N=pd.read_csv(OUT+'results/symptom_cluster_copula_null.csv').query("scores=='all12'")
A=pd.read_csv(OUT+'results/symptom_assignments_tp1.csv',index_col=0)
fig,ax=plt.subplots(1,3,figsize=(12,3.6))
ax[0].plot(N.k,N.stab,'-o',color=C1,lw=2,ms=6,label='real data'); ax[0].plot(N.k,N.null_stab_mean,'-o',color=GR,lw=2,ms=6,label='no-cluster null (mean)')
ax[0].fill_between(N.k,0,N.null_stab_95,color=MU,alpha=.35,label='null 95th pct'); ax[0].set_ylim(0,1); ax[0].set_xlabel('k (k-means)'); ax[0].set_ylabel('subsample stability (ARI)'); ax[0].set_title('Symptom clusters: stability',loc='left'); ax[0].legend(fontsize=7,loc='lower left')
ax[1].plot(N.k,N.sil,'-o',color=C1,lw=2,ms=6,label='real'); ax[1].plot(N.k,N.null_sil_mean,'-o',color=GR,lw=2,ms=6,label='null mean'); ax[1].fill_between(N.k,0,N.null_sil_95,color=MU,alpha=.35,label='null 95th pct')
ax[1].set_ylim(0,.25); ax[1].set_xlabel('k'); ax[1].set_ylabel('silhouette'); ax[1].set_title('Symptom clusters: separation',loc='left'); ax[1].legend(fontsize=7,loc='lower left')
cols={0:C3,1:C1,2:C2}; names={0:'milder',1:'moderate, low emotional burden',2:'severe, global'}
for j in range(3):
    s=A[A.sym_k3==j]; ax[2].scatter(s.symPC1,s.symPC2,s=16,color=cols[j],edgecolor='white',lw=.5,label=f'{names[j]} (n={len(s)})')
ax[2].set_xlabel('symptom PC1 (overall severity, 26%)'); ax[2].set_ylabel('symptom PC2 (emotional vs physical, 14%)'); ax[2].set_title('k=3 partition cuts a continuum',loc='left'); ax[2].legend(fontsize=7,loc='upper left')
plt.tight_layout(); plt.savefig(F+'fig1_symptom_clusters.png',dpi=150); plt.close()
# Fig 3: axes
T=pd.read_csv(OUT+'results/axes_tp1_values.csv',index_col=0)
AX=[('axis_lipolysis','Lipolysis (FFA+BHB+acylcarnitine)'),('axis_kyn_inflam','Kynurenine/IFN-g (KT ratio+IFN-g T)'),('axis_dysbiosis','Dysbiosis (R. gnavus - butyrate producers)')]
fig,ax=plt.subplots(1,4,figsize=(15,3.4),gridspec_kw={'width_ratios':[1,1,1,1.25]})
for i,(a,t) in enumerate(AX):
    p=T.loc[T.y==1,a].dropna(); c=T.loc[T.y==0,a].dropna(); b=np.linspace(min(p.min(),c.min()),max(p.max(),c.max()),26)
    ax[i].hist(c,bins=b,color=GR,alpha=.55,density=True,label=f'controls (n={len(c)})'); ax[i].hist(p,bins=b,histtype='step',color=C2,lw=2,density=True,label=f'patients (n={len(p)})')
    ax[i].set_title(t,loc='left',fontsize=8.5); ax[i].set_xlabel('axis score (SD units)'); ax[i].set_yticks([])
ax[0].legend(fontsize=7,loc='upper left')
R=pd.read_csv(OUT+'results/axes_correlations.csv'); comp=['FFA','acylcarnitine','BHB','KT_ratio','IFNg_T','MAIT_CD8','gd_T','butyrate_producers','R_gnavus','NLT.Orthostatic Intolerance','severity']
lab=['FFA','acylcarn.','BHB','Kyn/Trp','IFN-g T','MAIT','gd T','butyrate prod.','R. gnavus','OI (NLT)','severity']
M=np.full((len(comp),len(comp)),np.nan); r=R[R.group=='patients']
for _,x in r.iterrows():
    if x.a in comp and x.b in comp: i,j=comp.index(x.a),comp.index(x.b); M[i,j]=M[j,i]=x.partial_rho
np.fill_diagonal(M,np.nan)
from matplotlib.colors import LinearSegmentedColormap
cm=LinearSegmentedColormap.from_list('d',['#2a78d6','#f0efec','#e34948'])
im=ax[3].imshow(M,cmap=cm,vmin=-.6,vmax=.6); ax[3].set_xticks(range(len(lab))); ax[3].set_xticklabels(lab,rotation=60,ha='right',fontsize=7); ax[3].set_yticks(range(len(lab))); ax[3].set_yticklabels(lab,fontsize=7)
for i in range(len(comp)):
    for j in range(len(comp)):
        if i!=j and not np.isnan(M[i,j]) and abs(M[i,j])>=0.2: ax[3].text(j,i,f'{M[i,j]:.2f}'.replace('0.','.'),ha='center',va='center',fontsize=5.5,color=INK)
ax[3].set_title('Patients: partial Spearman (age/sex/BMI)',loc='left',fontsize=8.5); plt.colorbar(im,ax=ax[3],fraction=.04)
plt.tight_layout(); plt.savefig(F+'fig3_candidate_axes.png',dpi=150); plt.close()
print('ok')
# Fig 2: omics clustering vs null + gut clusters patients vs controls
S=pd.read_csv(OUT+'results/omics_cluster_selection.csv')
fig,ax=plt.subplots(1,3,figsize=(13,3.6))
c=S[(S.space=='concat4')&(S.method=='kmeans')]
ax[0].plot(c.k,c.sil,'-o',color=C1,lw=2,ms=6,label='real (4-layer, k-means)'); ax[0].plot(c.k,c.null_sil,'-o',color=GR,lw=2,ms=6,label='no-cluster null mean'); ax[0].fill_between(c.k,0,c.null_sil95,color=MU,alpha=.35,label='null 95th pct')
ax[0].set_ylim(0,.1); ax[0].set_xlabel('k'); ax[0].set_ylabel('silhouette'); ax[0].set_title('Multi-omics: no cluster structure',loc='left'); ax[0].legend(fontsize=7,loc='upper right')
sp=['concat4','only_immune','only_specie','only_metabolome','only_quest','SNF4']; nm=['4 layers','immune','species','metabolome','routine labs','SNF fusion']
for kk,off,col in [(2,-.2,C1),(3,.2,C2)]:
    r=[S[(S.space==s)&(S.k==kk)&(S.method.isin(['kmeans','spectral']))].iloc[0] for s in sp]
    x=np.arange(len(sp))+off
    ax[1].bar(x,[q.stab-q.null_stab for q in r],width=.36,color=col,label=f'k={kk}')
ax[1].axhline(0,color=GR,lw=1); ax[1].set_xticks(range(len(sp))); ax[1].set_xticklabels(nm,rotation=30,ha='right'); ax[1].set_ylabel('stability minus null stability (ARI)'); ax[1].set_title('Above-null stability: gut k=3, SNF k=3 only',loc='left'); ax[1].legend(fontsize=7)
g=pd.read_csv(OUT+'results/gut_k3_assignments.csv',index_col=0)
pat=np.array([23,58,46]); ctl=np.array([33,42,17])
x=np.arange(3); ax[2].bar(x-.2,100*pat/pat.sum(),.38,color=C2,label=f'patients (n={pat.sum()})'); ax[2].bar(x+.2,100*ctl/ctl.sum(),.38,color=GR,label=f'controls (n={ctl.sum()})')
ax[2].set_xticks(x); ax[2].set_xticklabels(['G1 control-like\n(F. plautii low)','G2 intermediate','G3 dysbiotic\n(R. gnavus, F. plautii up;\nFaecalibacterium down)'],fontsize=7.5); ax[2].set_ylabel('% of group'); ax[2].set_title('Gut clusters: patients vs controls',loc='left'); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(F+'fig2_omics_clusters.png',dpi=150); plt.close()
# Fig 4: temporal persistence + trajectory
Tm=pd.read_csv(OUT+'results/symptom_partition_temporal.csv'); import ast
t=np.array(ast.literal_eval(Tm[(Tm.k==3)&(Tm.tp==2)].table.iloc[0]))
Tr=pd.read_csv(OUT+'results/trajectory_table.csv',index_col=0)
fig,ax=plt.subplots(1,2,figsize=(9,3.6))
from matplotlib.colors import LinearSegmentedColormap
seq=LinearSegmentedColormap.from_list('s',['#f0efec','#2a78d6'])
ax[0].imshow(t,cmap=seq); 
for i in range(3):
    for j in range(3): ax[0].text(j,i,t[i,j],ha='center',va='center',color='white' if t[i,j]>20 else INK)
lb=['milder','moderate','severe']; ax[0].set_xticks(range(3)); ax[0].set_xticklabels(lb); ax[0].set_yticks(range(3)); ax[0].set_yticklabels(lb)
ax[0].set_xlabel('assigned partition at visit 2'); ax[0].set_ylabel('partition at visit 1'); ax[0].set_title('Symptom partition persistence\n(kappa 0.55; mostly = stable severity)',loc='left',fontsize=9)
for s in ['top','right','left','bottom']: ax[0].spines[s].set_visible(False)
for j,col in zip(range(3),[C3,C1,C2]):
    v=Tr[Tr.k3==j].d_sev; ax[1].scatter(np.full(len(v),j)+np.random.default_rng(j).uniform(-.15,.15,len(v)),v,s=14,color=col,edgecolor='white',lw=.5)
    ax[1].plot([j-.25,j+.25],[v.mean()]*2,color=INK,lw=2)
ax[1].axhline(0,color=GR,lw=1); ax[1].set_xticks(range(3)); ax[1].set_xticklabels([f'{l}\n(n={(Tr.k3==j).sum()})' for j,l in enumerate(lb)]); ax[1].set_ylabel('change in severity, first to last visit\n(negative = improved)')
ax[1].set_title('Trajectory by baseline partition',loc='left',fontsize=9)
plt.tight_layout(); plt.savefig(F+'fig4_temporal.png',dpi=150); plt.close()
print('figs done')
