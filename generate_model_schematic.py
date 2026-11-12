"""Generate a vector PDF schematic from the current model defaults."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from model import Parameters, simulate

OUT=Path(__file__).resolve().parent
p=Parameters()
df,meta=simulate(p,treated=False)
r=df.iloc[0]
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(16,11))
fig.subplots_adjust(left=.025,right=.975,bottom=.025,top=.975)
ax.set(xlim=(0,16),ylim=(0,11));ax.axis('off')
ink='#19334a'; muted='#516478'; blue='#176b99'; green='#287961'; red='#b03c45'

def text(x,y,s,size=10,color=ink,weight='normal',ha='center',va='center'):
    return ax.text(x,y,s,fontsize=size,color=color,fontweight=weight,ha=ha,va=va,linespacing=1.45)
def box(x,y,w,h,title,detail='',color=blue,fill='#eef6fb'):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.045,rounding_size=0.12',facecolor=fill,edgecolor=color,linewidth=1.4))
    text(x,y+(.13 if detail else 0),title,11,color,'bold')
    if detail:text(x,y-.21,detail,8.8,color)
def arrow(a,b,color=muted,both=False,style='solid',curve=0):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='<->' if both else '-|>',mutation_scale=13,linewidth=1.5,color=color,linestyle=style,connectionstyle=f'arc3,rad={curve}'))

text(.3,10.55,'sFLT1 degrader model',24,weight='bold',ha='left')
text(.3,10.02,'VEGF-A and PlGF distribution, reversible storage, and ligand fate during sFLT1 removal',12,muted,ha='left')
box(8,8.9,3.1,.8,'Free sFLT1  [S]',f'Production qS  |  natural loss kS',color=ink,fill='#eff2f5')
text(8,8.12,'Shared binding sites: VEGF and PlGF compete for S',10,muted).set_bbox(dict(facecolor='white',edgecolor='none',pad=3))
arrow((5.15,8.9),(6.38,8.9),green)
text(5.08,9.16,'qS input',9,green)
arrow((9.62,8.9),(10.8,8.9))
text(10.9,8.9,'Natural clearance\n+ degrader removal dS',9,ha='left')

for center,ligand,label,color in [(4,'vegf','VEGF-A',blue),(12,'plgf','PlGF',green)]:
    # Light lane boundary keeps independent binder/storage capacities distinct.
    ax.add_patch(FancyBboxPatch((center-3.8,4.16),7.6,3.5,boxstyle='round,pad=0.02,rounding_size=0.14',facecolor='#fafbfc',edgecolor='#dbe3ea',zorder=-1))
    box(center,7.05,2.65,.78,f'sFLT1–{label}  [C]', 'Natural clearance kC  |  degrader dC',color=color)
    box(center,5.2,2.05,.88,f'Free {label}  [L]', 'Available to bind sFLT1',color=color)
    box(center-2.65,5.2,2.05,.88,'Other-bound  [O]','Independent binding-site pool',color=color)
    box(center+2.65,5.2,2.05,.88,('Platelet-stored  [T]' if ligand=='vegf' else 'Placental-stored  [T]'),'Inaccessible while stored',color=color)
    arrow((center,5.69),(center,6.60),color,both=True)
    text(center+.16,6.14,'kon · S · L\nkoff · C',8.6,color,ha='left')
    arrow((center-1.05,5.2),(center-1.58,5.2),color,both=True)
    text(center-1.32,5.98,'Binding / release',8,color)
    arrow((center+1.05,5.2),(center+1.58,5.2),color,both=True)
    text(center+1.32,5.98,'Uptake / release',8,color)
    text(center+2.65,6.49,'Additional release rT → L\n(default zero)',8,color)
    # Replenishment and free production are external sources.
    arrow((center,4.18),(center,4.72),green)
    text(center-.12,4.38,'qL',9,green,ha='right')
    arrow((center+2.65,4.18),(center+2.65,4.72),green)
    text(center+2.48,4.38,'qT',9,green,ha='right')
    arrow((7.0 if center==4 else 9.0,8.46),(center+.35 if center==4 else center-.35,7.49),muted,both=True)
    text(center-2.65,7.1,label,14,color,'bold')

# A single shared fate panel applies identically to both ligand complexes.
box(8,3.42,15.15,.98,'DEGRADER: removes free sFLT1 and sFLT1–ligand complexes',
    'From each complex: sFLT1 is eliminated; 25% of ligand is co-degraded and 75% returns to free ligand (current defaults).',red,'#fff2f1')
text(8,2.62,'dC = d × relative activity on bound sFLT1   |   Other-bound and storage pools have no direct degrader clearance.',10,red)

# Compact parameter summary, generated rather than copied from past results.
rows=[['Starting concentration (pg/mL)','VEGF-A','PlGF'],
      ['Free plasma',f'{r.free_vegf_pg_mL:.2f}',f'{r.free_plgf_pg_mL:.2f}'],
      ['sFLT1-bound / other-bound / stored (each)',f'{r.bound_vegf_pg_mL:.2f}',f'{r.bound_plgf_pg_mL:.2f}'],
      ['Plasma total / all-pools total',f'{r.total_vegf_pg_mL:.2f} / {r.all_pools_vegf_pg_mL:.2f}',f'{r.total_plgf_pg_mL:.2f} / {r.all_pools_plgf_pg_mL:.2f}'],
      ['Kd: sFLT1 / other / effective storage (pM)',f'{p.kd_v_pm:g} / {p.other_kd_pm:g} / {p.platelet_kd_pm:g}',f'{p.kd_p_pm:g} / {p.plgf_other_kd_pm:g} / {p.plgf_placental_kd_pm:g}']]
table=ax.table(cellText=rows,cellLoc='left',colWidths=[.34,.14,.14],bbox=[.015,.068,.62,.152])
table.auto_set_font_size(False);table.set_fontsize(8.7)
for (i,j),cell in table.get_celld().items():
    cell.set_edgecolor('#dbe3ea');cell.set_facecolor('#edf2f6' if i==0 else 'white')
    if i==0:cell.set_text_props(weight='bold',color=ink)
text(10.55,2.18,'NATURAL LOSS & INPUTS',10,weight='bold',ha='left')
text(10.55,1.65,'Every state has first-order natural loss.\nComplex loss destroys ligand; other-site capacity is maintained.\nqS, qL and qT continue throughout treatment.',8.7,ha='left')
text(10.55,.94,'Production is calibrated to initial untreated steady state.\nSingle dose day 0: initial 2/day, half-life 3 days.\nStorage and affinity are effective assumptions.',8.7,ha='left')
text(.3,.3,'Plasma total = L + C + O; all pools = L + C + O + T.  Independent O/T sites for each ligand; only sFLT1 is shared.',9,muted,ha='left')
text(.3,.02,'Exploratory model, not a validated physiological or assay prediction.  Source: model.py defaults; regeneration script: generate_model_schematic.py.',8,muted,ha='left')
fig.savefig(OUT/'model_schematic.pdf',metadata={'Title':'sFLT1 degrader model schematic','Subject':'VEGF and PlGF binding, platelet VEGF / placental PlGF storage and co-degradation'})
fig.savefig(OUT/'model_schematic.png',dpi=150)
plt.close(fig)
print(OUT/'model_schematic.pdf')
