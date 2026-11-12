"""One-at-a-time VEGF/PlGF sensitivity analyses with nine separate figures.

Run: MPLCONFIGDIR=.mplconfig .venv/bin/python sensitivity_analysis.py
"""
from dataclasses import asdict, replace
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from model import Parameters, simulate

OUT=Path(__file__).resolve().parent/'sensitivity_results'
BASE=Parameters()
GROUPS=[
    ('vegf_half_life','Free VEGF half-life',[(f'{h} hour'+('s' if h!=1 else ''),dict(vegf_half_hours=h),BASE.degrader_half_days) for h in [1,4,8,12]],
     'Single dose: initial clearance 2/day, degrader half-life 1 day; 25% co-degradation'),
    ('degrader_half_life','Degrader half-life',[('2 hours',{},2/24),('12 hours',{},12/24)]+[(f'{h} day'+('s' if h!=1 else ''),{},h) for h in [1,3,7]],
     'Single dose: initial clearance 2/day, exponentially decaying; 25% co-degradation'),
    ('codegradation','Ligand co-degradation',[(f'{int(f*100)}%',dict(co_v=f,co_p=f),BASE.degrader_half_days) for f in [.25,.5,.75]],
     'Same fraction for VEGF and PlGF; single dose, initial clearance 2/day, half-life 1 day'),
]
COLORS=['#176b99','#d17b18','#287961','#884da0','#c14465']

def main():
    OUT.mkdir(exist_ok=True)
    summary=[];all_data=[];scenario_meta=[]
    baseline,_=simulate(BASE,treated=False)
    b=baseline.iloc[0]
    baseline_rows=[
        ('Dose / duration','Single dose day 0 / 14 days'),
        ('Initial degrader clearance',f'{BASE.degrader_per_day:g}/day'),
        ('Degrader half-life',f'{BASE.degrader_half_days:g} day'+('s' if BASE.degrader_half_days!=1 else '')),
        ('Co-degradation: VEGF / PlGF',f'{BASE.co_v:.0%} / {BASE.co_p:.0%}'),
        ('Free VEGF / PlGF half-life',f'{BASE.vegf_half_hours*60:.0f} min / {BASE.plgf_half_hours:.2f} h'),
        ('Initial total sFLT1',f'{BASE.sflt1_ng_ml:g} ng/mL'),
        ('Initial free VEGF / PlGF',f'{b.free_vegf_pg_mL:.2f} / {b.free_plgf_pg_mL:.2f} pg/mL'),
        ('All-pools VEGF / PlGF',f'{b.all_pools_vegf_pg_mL:.0f} / {b.all_pools_plgf_pg_mL:.0f} pg/mL'),
        ('sFLT1 / complex half-life',f'{BASE.sflt1_half_hours:g} / {BASE.sv_half_hours:g} h'),
        ('Kd: VEGF / PlGF (all pools)',f'{BASE.kd_v_pm:g} / {BASE.kd_p_pm:g} pM'),
        ('Other-bound / storage loss half-life',f'{BASE.other_complex_half_hours:g} / {BASE.platelet_half_hours:g} h'),
        ('Production','Recalibrate untreated steady state'),
    ]
    pd.DataFrame(baseline_rows,columns=['Parameter','Default']).to_csv(OUT/'baseline_parameters.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
    with PdfPages(OUT/'sensitivity_analysis_all_plots.pdf') as report:
        for group,title,cases,note in GROUPS:
            scenarios=[]
            for label,changes,half in cases:
                p=replace(BASE,**changes)
                df,meta=simulate(p,degrader_half_days=half)
                control,_=simulate(p,treated=False,degrader_half_days=half)
                # True initial untreated steady state in every scenario.
                for ligand in ['vegf','plgf']:
                    key=f'free_{ligand}_pg_mL'
                    np.testing.assert_allclose(control[key],control[key].iloc[0],rtol=1e-6)
                    y=df[key].to_numpy(); peak=int(np.argmax(y))
                    summary.append(dict(analysis=group,scenario=label,ligand=ligand,unit='pg/mL',baseline_pg_mL=y[0],day14_pg_mL=y[-1],day14_percent_change=100*(y[-1]/y[0]-1),peak_pg_mL=y[peak],peak_day=df.day.iloc[peak],production_vegf_pg_mL_day=meta['production_v_pM_day']*p.vegf_mw_kda,production_plgf_pg_mL_day=meta['production_p_pM_day']*p.plgf_mw_kda))
                export=df.copy()
                export.insert(0,'scenario',label);export.insert(0,'analysis',group)
                export['untreated_total_sflt1_ng_mL']=control.total_sflt1_ng_mL.to_numpy()
                export['untreated_free_vegf_pg_mL']=control.free_vegf_pg_mL.to_numpy()
                export['untreated_free_plgf_pg_mL']=control.free_plgf_pg_mL.to_numpy()
                export['degrader_clearance_per_day']=p.degrader_per_day*(np.exp(-np.log(2)*df.day/half) if half else np.ones(len(df)))
                all_data.append(export)
                scenario_meta.append(dict(analysis=group,scenario=label,parameters=asdict(p),derived=meta))
                scenarios.append((label,df,control))
            for ligand,label in [('vegf','Free VEGF-A'),('plgf','Free PlGF'),('sflt1','Total sFLT1')]:
                fig,ax=plt.subplots(figsize=(10,6.4),layout='constrained')
                key='total_sflt1_ng_mL' if ligand=='sflt1' else f'free_{ligand}_pg_mL'
                unit='ng/mL' if ligand=='sflt1' else 'pg/mL' 
                for color,(scenario,df,control) in zip(COLORS,scenarios):
                    ax.plot(df.day,df[key],color=color,lw=2.5,label=scenario)
                control=scenarios[0][2]
                if ligand=='sflt1':
                    for scenario,df,_ in scenarios:
                        y=df[key].to_numpy(); nadir=int(np.argmin(y))
                        summary.append(dict(analysis=group,scenario=scenario,ligand='sflt1',unit='ng/mL',baseline_ng_mL=y[0],day14_ng_mL=y[-1],day14_percent_change=100*(y[-1]/y[0]-1),nadir_ng_mL=y[nadir],nadir_day=df.day.iloc[nadir]))
                ax.plot(control.day,control[key],color='#59616c',ls='--',lw=1.6,label='Untreated steady-state baseline')
                ax.set(xlabel='Time after treatment starts (days)',ylabel=f'{label} ({unit})',xlim=(0,BASE.days))
                fig.suptitle(f'{label}: sensitivity to {title}',fontsize=15,fontweight='bold')
                ax.text(.5,1.025,note,transform=ax.transAxes,ha='center',fontsize=9,color='#516478')
                if ligand=='sflt1' and group!='degrader_half_life':
                    ax.text(.02,.04,'Scenario curves overlap: total sFLT1 turnover is unchanged by these ligand parameters.',transform=ax.transAxes,fontsize=8,color='#516478',bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
                ax.grid(alpha=.2);ax.legend(loc='best',framealpha=.95)
                ax.ticklabel_format(axis='y',style='plain',useOffset=False)
                fig.supxlabel('Production recalibrated to initial untreated steady state for each scenario; all other parameters held at current defaults.\nVertical scales differ across plots. Model predictions; storage pools are exploratory assumptions.',fontsize=8,color='#516478')
                stem=OUT/(f'{group}_total_sflt1' if ligand=='sflt1' else f'{group}_free_{ligand}')
                fig.savefig(stem.with_suffix('.png'),dpi=200)
                fig.savefig(stem.with_suffix('.pdf'))
                report.savefig(fig)
                plt.close(fig)
    pd.DataFrame(summary).to_csv(OUT/'sensitivity_summary.csv',index=False)
    pd.concat(all_data,ignore_index=True).to_csv(OUT/'sensitivity_time_courses.csv',index=False)
    (OUT/'sensitivity_parameters.json').write_text(json.dumps(scenario_meta,indent=2))
    (OUT/'README.md').write_text('''# Sensitivity results

Nine separate plots (PNG and vector PDF), including total sFLT1 for each analysis, and a combined nine-page PDF. Individual plots have no parameter tables. The sensitivity tab retains one baseline-default table at the beginning (also exported as baseline_parameters.csv). The parameter being varied supersedes its baseline value in each scenario.

- Free VEGF half-life: 1, 4, 8, 12 hours. PlGF free half-life stays at its default (16.64 hours).
- Degrader half-life: 2 hours, 12 hours, 1 day, 3 days, 7 days, with single-dose initial effective clearance of 2/day. Effective clearance decays as 2*exp(-ln(2)*t/half-life). No redosing. This assumes activity is proportional to exposure, not a saturable PK/PD model.
- Co-degradation fraction: 25%, 50%, 75%, changed together for both ligands. The remainder returns to free ligand.

These are separate one-at-a-time analyses, not a full factorial combination. All analyses use a single dose on day 0 with initial clearance 2/day and exponential decay. Degrader half-life is 1 day in the VEGF-half-life and co-degradation analyses. Both ligand co-degradation fractions stay at 25% in the other analyses. Free VEGF half-life stays at 34 minutes except in its own analysis. Treatment begins on day 0; simulation and treatment windows end at day 14.

For each scenario, production is recalibrated to the initial untreated kinetic steady state, then kept fixed throughout the simulation. Initial distributions, other binding, platelet VEGF storage and placental PlGF storage retain their current default values. Untreated curves coincide across scenarios and are shown as a single dashed baseline. Production differs across VEGF-half-life scenarios; this is not a fixed-production comparison.

Summary CSV reports ligand day-14 and peak concentrations with peak timing, plus total sFLT1 endpoints and nadirs (ng/mL); decaying degrader responses can peak before day 14. Full treated/untreated time courses and effective clearance are in the time-course CSV. Full defaults and effective inputs are in the parameters JSON. Units for concentrations are pg/mL and time is days. All values are exploratory model predictions, not validated physiological estimates.

Regenerate from the project folder: `MPLCONFIGDIR=.mplconfig .venv/bin/python sensitivity_analysis.py` (requires requirements-dev.txt).
''')
    print(pd.DataFrame(summary)[['analysis','scenario','ligand','day14_pg_mL','day14_percent_change','peak_pg_mL','peak_day']].round(4).to_string(index=False))
    print(f'\nSaved nine plots, combined PDF, CSVs and parameters to {OUT}')

if __name__=='__main__': main()
