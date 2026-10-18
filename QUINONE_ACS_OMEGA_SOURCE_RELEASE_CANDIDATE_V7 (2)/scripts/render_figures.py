#!/usr/bin/env python3
"""Preview figures from FROZEN, provenance-labeled ACS Omega quinone inputs.
This code performs no new quantum-chemical calculations.
"""
import os, math, csv, json, pathlib, statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import cairosvg
HERE=pathlib.Path(__file__).resolve().parent
DATA=HERE.parent/'data'
FIG=HERE.parent/'figures';FIG.mkdir(exist_ok=True)

def rows(name):
 with open(DATA/name, newline='') as f:return list(csv.DictReader(f))
def save(fig,name):
 fig.savefig(FIG/(name+'.png'),dpi=230, bbox_inches='tight',facecolor='white')
 fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',facecolor='white')
 fig.savefig(FIG/(name+'.svg'),bbox_inches='tight',facecolor='white')
 plt.close(fig)

def source_models():
 from rdkit.Chem import Draw
 c=rows('M6_12_CHEMICAL_IDENTITY_CROSSWALK.csv')
 c={r['kind']:r for r in c if r['site']=='Q0527_R0_SME1'}
 legend=['Parent BAQ','Methylthio (SMe)','Capped-cysteine (AcCysOMe)']
 molecules=[Chem.MolFromSmiles(c[k]['smiles']) for k in ('parent','sme','accysome')]
 for mol in molecules:
  Chem.rdDepictor.Compute2DCoords(mol)
 im=Draw.MolsToGridImage(molecules,molsPerRow=3,subImgSize=(740,520),legends=legend,useSVG=False)
 im.convert('RGB').save(FIG/'FIG01_CHEMICAL_MODELS.png')
 svg=Draw.MolsToGridImage(molecules,molsPerRow=3,subImgSize=(740,520),legends=legend,useSVG=True)
 (FIG/'FIG01_CHEMICAL_MODELS.svg').write_text(svg)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.titlesize':12})
source_models()
# Fig 2: M1 broad exploratory distribution + family-level variation
m1=rows('M1_SHIFTS.csv'); d=np.array([float(r['delta_delta_Ered_eV']) for r in m1]);par={}
for r in m1:par.setdefault(r['parent_id'],[]).append(float(r['delta_delta_Ered_eV']))
mult=[(p,sorted(v)) for p,v in par.items() if len(v)>1]
mult.sort(key=lambda x:max(x[1])-min(x[1]), reverse=True)
fig, (a,b) = plt.subplots(1,2,figsize=(11.2,4.2),gridspec_kw={'width_ratios':[1,1.24]},layout='constrained')
xs=np.arange(len(d)); a.scatter(xs,np.sort(d),s=17)
a.axhline(0,color='0.4',lw=.8);a.set_title('A  Site-specific SMe shifts (n = 38)',loc='left',fontweight='bold');a.set_xlabel('Sites sorted by calculated shift');a.set_ylabel('Parent-relative shift (eV)');a.set_xticks([0,9,19,28,37]);a.set_xlim(-1,38)
for i,(p,vals) in enumerate(mult):
 b.plot([min(vals),max(vals)],[i,i],lw=2.1,alpha=.8)
 b.scatter(vals,[i]*len(vals),s=28,zorder=3)
b.axvline(0,color='0.4',lw=.8)
b.set_yticks(range(len(mult)),[p for p,_ in mult]);b.invert_yaxis();b.set_xlabel('Calculated site shift (eV)');b.set_title('B  Nine multi-site parents',loc='left',fontweight='bold')
fig.suptitle('M1 exploratory site heterogeneity: sulfur-substitution response is not uniform',fontweight='bold',fontsize=12)
save(fig,'FIG02_M1_SITE_HETEROGENEITY')
# Fig 3: all ten original M2R and later 4 corrected, provenance explicitly labeled
m2=rows('M2R_CORE.csv');m3=rows('M3R_SHIFTS.csv');cmap={r['site_id']:r for r in m3};
fig,ax=plt.subplots(figsize=(7,5.8),layout='constrained')
orig=[]
for r in m2:
 x=float(r['sme_shift_eV']); y=float(r['accysome_shift_eV']);site=r['site_id'];orig.append((x,y));
 ax.scatter(x,y,marker='o',s=43,alpha=.78,color='C0')
 if site in cmap:
  corr=float(cmap[site]['m3r_final_shift_eV'])
  ax.annotate('',xy=(x,corr),xytext=(x,y),arrowprops={'arrowstyle':'-|>','lw':1,'color':'0.40','alpha':0.85})
  ax.scatter(x,corr,marker='D',facecolors='none',edgecolors='black',s=70,lw=1.2,zorder=6)
  ax.annotate({'Q0527_R0_SME1':'BAQ','Q0298_R3_SME1':'7-MeBAQ','Q0059_R0_SME3':'BPAQ-C5','Q0528_R0_SME1':'NPQ'}[site],(x,corr),textcoords='offset points',xytext=({'Q0527_R0_SME1':(5,-14),'Q0298_R3_SME1':(5,10),'Q0059_R0_SME3':(6,9),'Q0528_R0_SME1':(6,8)}[site]),fontsize=8)
lim1=-.40;lim2=.29
ax.plot([lim1,lim2],[lim1,lim2],ls='--',color='0.6',lw=.9,label='Perfect transfer (y = x)')
ax.axhline(0,color='0.4',lw=.7);ax.axvline(0,color='0.4',lw=.7)
ax.set(xlim=(lim1,lim2),ylim=(-.44,.32),xlabel='SMe shift (eV)',ylabel='AcCysOMe shift (eV)',title='M2R original ten-site screen and later four-case corrections')
ax.plot([],[], 'o',label='Original M2R (all ten)')
ax.plot([],[], 'D',markerfacecolor='none',markeredgecolor='black',label='Corrected M3R (four)')
ax.legend(fontsize=8,loc='lower right');save(fig,'FIG03_M2R_VINTAGE_CORRECTIONS')
# Fig 4: Operational M5 decomposition
m5=rows('M5_DECOMPOSITION.csv');order=['Q0527_R0_SME1','Q0298_R3_SME1','Q0059_R0_SME3','Q0528_R0_SME1'];lookup={r['site_id']:r for r in m5};labels=['BAQ','7-MeBAQ','BPAQ-C5','NPQ control']
fig,ax=plt.subplots(figsize=(8.2,4.2),layout='constrained');ys=np.arange(4)
geo=np.array([float(lookup[s]['geometry_mediated_component_eV']) for s in order]);side=np.array([float(lookup[s]['sidechain_chemical_component_eV']) for s in order]);
for i in range(4):
 for value,shift,color in [(geo[i],-.12,'C0'),(side[i],.12,'C1')]:
  ax.plot([0,value],[i+shift]*2,lw=3.0,color=color,solid_capstyle='round')
  ax.scatter([value],[i+shift],s=36,color=color,zorder=3)
 ax.annotate(f'Total {geo[i]+side[i]:+.3f} eV',xy=(.21,i),fontsize=8,va='center')
ax.plot([],[],lw=3,color='C0',label='Geometry-dependent term')
ax.plot([],[],lw=3,color='C1',label='Side-chain replacement term')
ax.legend(loc='lower left',frameon=False,fontsize=9,bbox_to_anchor=(0,1.01),ncol=2)
ax.axvline(0,color='0.3',lw=.85);ax.axvline(-.05,color='0.70',lw=.7,ls=':');ax.axvline(.05,color='0.70',lw=.7,ls=':');ax.set_xlim(-.33,.35);ax.set_yticks(ys,labels);ax.invert_yaxis();ax.set_xlabel('Operational component (eV)');ax.set_title('M5: counterfactual geometry and side-chain contributions',loc='left',fontweight='bold',pad=35)
save(fig,'FIG04_M5_COUNTERFACTUAL_PARTITION')
# Fig 5: two method panels (within-method surrogate comparisons only)
m6=rows('M5_M6_COMPARISON.csv'); lookup6={r['site']:r for r in m6}
fig,axes=plt.subplots(1,2,figsize=(10.3,4.3),sharey=True,layout='constrained')
for ax,method,sfield,afield in [(axes[0],'Original method / M3R panel','m5_sme_shift_eV','m3r_accysome_shift_eV'),(axes[1],'M6 fixed-geometry PBE0-D4','m6_sme_shift_eV','m6_accysome_shift_eV')]:
 for i,(site,lab) in enumerate(zip(order,labels)):
  r=lookup6[site];x=float(r[sfield]);y=float(r[afield])
  ax.plot([x,y],[i,i],lw=1.5,color='.55',zorder=1)
  ax.scatter(x,i-.075,s=43,color='C0',label='SMe' if i==0 else None,zorder=3)
  ax.scatter(y,i+.075,s=43,color='C1',marker='D',label='AcCysOMe' if i==0 else None,zorder=3)
 ax.axvline(0,color='.2',lw=.8)
 ax.set_xlabel('Parent-relative electronic shift (eV)')
 ax.set_xlim(-.22,.38);ax.set_title(method,fontweight='bold',fontsize=10)
 ax.grid(axis='x',alpha=.18)
axes[0].set_yticks(ys,labels)
axes[0].invert_yaxis()
axes[1].tick_params(labelleft=False)
axes[0].legend(frameon=False,fontsize=9,loc='upper center',bbox_to_anchor=(1.05,-.12),ncol=2)
fig.suptitle('Surrogate responses across four prespecified quinone cases',fontsize=12,fontweight='bold')
save(fig,'FIG05_M6_FUNCTIONAL_SENSITIVITY')
print('Five proof figures generated. All artwork is draft; captions must state estimator limitations.')
