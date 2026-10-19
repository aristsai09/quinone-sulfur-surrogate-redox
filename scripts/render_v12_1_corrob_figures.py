#!/usr/bin/env python3
"""Rebuild two main figures from *unchanged* archival subset CSVs.
This script does NO quantum calculations. Figures are post hoc presentation of
source-corroborated records; historical full-cohort figures remain distinct.
"""
from pathlib import Path
import csv, math, statistics, json, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'data/v12_1'
OUT=ROOT/'figure_assets'; OUT.mkdir(exist_ok=True)
COLORS={'Q0059':'#4c6b89','Q0527':'#9b5a42','Q0298':'#8579b6','Q0528':'#5e8a70','Q0297':'#948852'}
BLUE='#276184'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.6,'figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
def load(name):
 with (SRC/name).open(newline='',encoding='utf-8') as f: return list(csv.DictReader(f))
m1=load('V11H_RAW_CORROBORATED_M1_34_SITES.csv')
m2=load('V11H_RAW_CORROBORATED_M2R_8_SITES.csv')
assert len(m1)==34 and len({r['parent_id'] for r in m1})==25
assert len(m2)==8 and len({r['parent_id'] for r in m2})==5
assert all(r['site_id'] in {s['adduct_id'] for s in m1} for r in m2)

# Figure 2: same basic grammar as historical picture; visibly primary 34 cases.
vals=sorted([(float(r['delta_delta_Ered_eV']),r['adduct_id'],r['parent_id']) for r in m1])
parents={}
for val,site,par in vals:parents.setdefault(par,[]).append(val)
multi={p:v for p,v in parents.items() if len(v)>1}
assert len(multi)==7
opposite=sum(min(v)<0 and max(v)>0 for v in multi.values());assert opposite==4
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.57),gridspec_kw={'width_ratios':[1.09,1]},constrained_layout=False)
fig.subplots_adjust(left=.065,right=.98,top=.84,bottom=.20,wspace=.22)
fig.suptitle('M1 source-corroborated sulfur-site shifts (34 sites; 25 parents)',fontsize=13,fontweight='bold',y=.965)
for i,(v,site,p) in enumerate(vals):
 ax1.scatter(i,v,s=22,c=BLUE if v<0 else '#b06342',zorder=3,edgecolor='white',linewidth=.3)
ax1.axhline(0,color='#555',lw=.8)
ax1.axhline(.05,color='#929292',lw=.7,ls='--');ax1.axhline(-.05,color='#929292',lw=.7,ls='--')
ax1.set(title='A  Parent-relative SMe electronic shifts',xlabel='Corroborated sites, ordered by shift',ylabel='Calculated shift, eV')
ax1.set_xlim(-1,34);ax1.set_xticks([0,9,19,29,33]);ax1.set_ylim(-.25,.275)
ax1.text(.02,.96,'25 negative  /  9 positive',transform=ax1.transAxes,va='top',fontsize=8)
for j,p in enumerate(sorted(multi,key=lambda k: max(multi[k])-min(multi[k]),reverse=True)):
 y=len(multi)-1-j; v=multi[p]; c=BLUE if min(v)<0 and max(v)>0 else '#7b8390'
 ax2.plot([min(v),max(v)],[y,y],lw=2.1,color=c,alpha=.8,zorder=2)
 ax2.scatter(v,[y]*len(v),c=c,s=29,zorder=4,edgecolors='white',lw=.35)
ax2.axvline(0,color='#555',lw=.8)
ax2.set(yticks=range(len(multi)),yticklabels=list(reversed(sorted(multi,key=lambda k:max(multi[k])-min(multi[k]),reverse=True))),xlabel='Calculated site shift, eV',title='B  Seven multisite parents')
ax2.set_xlim(-.11,.27)
ax2.text(.98,.03,'4/7 parents cross zero',va='bottom',ha='right',transform=ax2.transAxes,fontsize=8)
fig.text(.5,.075,'Retrospective source-provenance restriction; historical 38-site plot retained separately in SI.',ha='center',fontsize=7,color='#444')
for ext in ('png','pdf'):
 fig.savefig(OUT/f'FIG02_M1_CORROBORATED34.{ext}',dpi=300)
plt.close(fig)

# Figure 3: equal-status observation from 8 historical-vintage M2R rows (not mixed into M3R).
siteabbr={
 'Q0059_R0_SME3':'BPAQ C5','Q0527_R0_SME1':'BAQ site 1','Q0527_R0_SME2':'BAQ site 2',
 'Q0298_R3_SME1':'7-MeBAQ site 1','Q0298_R3_SME2':'7-MeBAQ site 2',
 'Q0528_R0_SME1':'NPQ site 1','Q0528_R0_SME2':'NPQ site 2','Q0297_R3_SME2':'DMBAQ site 2'}
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(9.2,7.55),gridspec_kw={'width_ratios':[1.2,1]},constrained_layout=False)
fig.subplots_adjust(left=.095,right=.98,top=.855,bottom=.17,wspace=.43)
fig.suptitle('M2R surrogate transfer: eight source-corroborated matched sites',fontsize=12,fontweight='bold',y=.965)
ax1.set_title('A  Same-site electronic reduction shifts')
lim=(-.27,.28);ax1.set_xlim(lim);ax1.set_ylim((-.48,.16))
xx=[-.32,.36]
ax1.plot(xx,xx,color='#646464',linestyle='--',linewidth=1,label='Identity')
ax1.fill_between(xx,[x-.05 for x in xx],[x+.05 for x in xx],color='#e3ebf0',alpha=.85,label='±0.05 eV tolerance')
ax1.axhline(0,color='#aaa',lw=.65);ax1.axvline(0,color='#aaa',lw=.65)
for r in m2:
 p=r['parent_id'];x=float(r['sme_shift_eV']);y=float(r['accysome_shift_eV'])
 ax1.scatter(x,y,c=COLORS[p],s=62,edgecolor='#222',linewidth=.45,zorder=4)
ax1.set(xlabel='Methylthio shift, eV',ylabel='Capped-cysteine shift, eV')
ax1.legend(loc='lower right',frameon=False)

ordered=sorted(m2,key=lambda r: float(r['abs_difference_eV']))
for i,r in enumerate(ordered):
 diff=float(r['abs_difference_eV']);p=r['parent_id']
 ax2.plot([0,diff],[i,i],color=COLORS[p],lw=1.9,alpha=.82)
 ax2.scatter(diff,i,color=COLORS[p],s=38,zorder=3,edgecolors='#222',linewidth=.3)
ax2.set(yticks=range(8),yticklabels=[siteabbr[r['site_id']] for r in ordered],xlabel='Absolute surrogate discrepancy, eV',title='B  Sitewise transfer mismatch')
ax2.axvline(.05,color='#666',linestyle='--',lw=1)
ax2.text(.053,7.52,'Frozen 0.05 eV criterion',fontsize=7,va='top',color='#444')
ax2.set_xlim(0,.49);ax2.set_ylim(-.7,8.1)
ax2.text(.98,.04,'Median = 0.162917 eV\n8/8 exceed 0.05 eV',transform=ax2.transAxes,ha='right',va='bottom',fontsize=8)
fig.text(.5,.07,'Retrospective provenance restriction; all values retain original M2R vintage (not four-case M3R corrections).',ha='center',fontsize=7,color='#444')
for ext in ('png','pdf'):
 fig.savefig(OUT/f'FIG03_M2R_CORROBORATED8.{ext}',dpi=300)
plt.close(fig)

# Strictly descriptive threshold check, not a new prespecified gate.
thresholds=[.05,.10,.20]
with (SRC/'V12_1_RETROSPECTIVE_THRESHOLD_SENSITIVITY.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['absolute_surrogate_discrepancy_threshold_eV','n_exceed_threshold_out_of_8','n_sites','parent_families_in_panel','analysis_status'])
 for t in thresholds: w.writerow([f'{t:.2f}',sum(float(r['abs_difference_eV'])>=t for r in m2),8,5,'post_hoc_descriptive'])
assert abs(statistics.median([float(r['abs_difference_eV']) for r in m2])-0.162916518)<1e-8
# Secondary mixed-vintage sensitivity; never report as homogeneous.
corrections={r['site_id']:float(r['m3r_final_shift_eV']) for r in load('M3R_SHIFTS.csv')}
mixed=[]
for r in m2:
 site=r['site_id'];acc=corrections.get(site,float(r['accysome_shift_eV']))
 mixed.append({'site_id':site,'parent_id':r['parent_id'],'sme_shift_original_eV':r['sme_shift_eV'],'accysome_shift_post_hoc_mixed_eV':str(acc),'abs_mismatch_post_hoc_mixed_eV':str(abs(acc-float(r['sme_shift_eV']))),'accysome_source_vintage':'M3R_corrected' if site in corrections else 'M2R_historical'})
with (SRC/'V12_1_MIXED_VINTAGE_8SITE_SENSITIVITY.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(mixed[0]));w.writeheader();w.writerows(mixed)
res={'n_m1':len(m1),'n_parents_m1':len(parents),'m1_neg':sum(v[0]<0 for v in vals),'m1_pos':sum(v[0]>0 for v in vals),'m1_multisite_parents':len(multi),'m1_opposite_sign_multisite':opposite,'m2_n':len(m2),'m2_parent_families':len({r['parent_id'] for r in m2}),'m2_median_abs_eV':statistics.median([float(r['abs_difference_eV']) for r in m2]),'threshold_exceed_count':{str(t):sum(float(r['abs_difference_eV'])>=t for r in m2) for t in thresholds},'mixed_vintage_median_abs_eV':statistics.median([float(r['abs_mismatch_post_hoc_mixed_eV']) for r in mixed]),'m3r_corrected_rows':sum(x['accysome_source_vintage']=='M3R_corrected' for x in mixed)}
assert sum(x['accysome_source_vintage']=='M3R_corrected' for x in mixed)==4
(SRC/'V12_1_FIGURE_AND_SENSITIVITY_VERIFICATION.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
