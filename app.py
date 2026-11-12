import json
from pathlib import Path
from dataclasses import asdict
import numpy as np
import pandas as pd
import altair as alt
import streamlit as st
from model import Parameters, setup, simulate
from default_parameters import default_parameter_table

st.set_page_config(page_title='sFLT1 degrader model',layout='wide')
st.title('sFLT1 degrader model')
st.caption('Explore free VEGF-A and PlGF following a single dose of an sFLT1 degrader. Illustrative single-compartment model using the supplied project PDF; parameters are assumptions, not fitted patient estimates.')
base=Parameters()
values=asdict(base)

def reset_defaults():
    for key in list(st.session_state):
        del st.session_state[key]

def chart(frame, ylabel, zero=True):
    # Explicit chart has no interval selection, so wheel events cannot zoom.
    data=frame.rename_axis('Day').reset_index().melt('Day',var_name='Series',value_name='Value')
    plot=alt.Chart(data).mark_line(strokeWidth=2,point=alt.OverlayMarkDef(size=8,opacity=0.65)).encode(
        x=alt.X('Day:Q',title='Day'),
        y=alt.Y('Value:Q',title=ylabel,scale=alt.Scale(zero=zero)),
        color=alt.Color('Series:N',legend=alt.Legend(orient='bottom')),
        tooltip=[alt.Tooltip('Day:Q',format='.2f'),'Series:N',alt.Tooltip('Value:Q',format='.4g')],
    ).properties(height=300)
    st.altair_chart(plot,use_container_width=True)

def number(key,label,minimum=0.0):
    values[key]=st.number_input(label,min_value=minimum,value=float(values[key]),format='%.6g',key=key)

model_tab, schematic_tab, sensitivity_tab, defaults_tab = st.tabs(['Interactive model', 'Model schematic', 'Sensitivity analysis', 'Default parameters'])

with defaults_tab:
    st.subheader('Default parameters')
    st.caption('Saved model defaults, independent of sidebar edits.')
    defaults = default_parameter_table()
    st.dataframe(defaults,hide_index=True,width='stretch')
    st.caption('Starting VEGF and PlGF totals include all binding and storage pools. Their plasma totals exclude storage: VEGF 68.83 pg/mL; PlGF 77.91 pg/mL. Kd describes ligand binding to sFLT1. Co-degradation is the fraction of bound ligand destroyed when sFLT1 is degraded; sFLT1 itself is removed.')
    st.download_button('Download default parameters (CSV)',defaults.to_csv(index=False),'default_parameters.csv','text/csv',key='download_default_parameters')

with schematic_tab:
    st.subheader('Model schematic')
    st.caption('Reference diagram showing the model structure and default parameters, including 25% co-degradation. It does not change when you adjust the simulation inputs.')
    assets = Path(__file__).resolve().parent
    st.image(str(assets / 'model_schematic.png'), use_container_width=True)
    st.download_button(
        'Download schematic as PDF',
        data=(assets / 'model_schematic.pdf').read_bytes(),
        file_name='model_schematic.pdf',
        mime='application/pdf',
        key='download_schematic',
    )

with sensitivity_tab:
    st.subheader('Sensitivity analysis')
    st.write('Separate 14-day analyses of free VEGF, free PlGF and total sFLT1, varying one parameter at a time. Production is recalibrated to initial untreated steady state for every scenario.')
    results = Path(__file__).resolve().parent / 'sensitivity_results'
    baseline_defaults = pd.read_csv(results / 'baseline_parameters.csv')
    left = baseline_defaults.iloc[:6].reset_index(drop=True)
    right = baseline_defaults.iloc[6:].reset_index(drop=True).rename(columns={'Parameter':'Parameter (continued)','Default':'Default (continued)'})
    st.markdown('**Baseline defaults for these saved analyses**')
    st.dataframe(pd.concat([left,right],axis=1),hide_index=True,use_container_width=True)
    st.caption('Each sensitivity varies the labeled parameter away from this reference; all other defaults remain fixed. Other-binder and storage affinities are effective assumptions.')
    st.info('These are saved reference results using model defaults, not the current sidebar inputs. All analyses use a single dose on day 0, initial clearance 2/day, and exponential activity decay. The reference degrader half-life is 3 days except when varied. Other settings retain their defaults, including 25% co-degradation when that parameter is not being varied.')
    results = Path(__file__).resolve().parent / 'sensitivity_results'
    for group, title, description in [
        ('vegf_half_life', 'Free VEGF half-life', '1, 4, 8 and 12 hours. Free PlGF half-life remains 16.64 hours.'),
        ('degrader_half_life', 'Degrader half-life', '2 hours, 12 hours, 1 day, 3 days and 7 days. A single dose starts at clearance 2/day; activity decays exponentially, without redosing.'),
        ('codegradation', 'Ligand co-degradation', '25%, 50% and 75%, changed together for VEGF and PlGF. The remaining ligand returns to the free pool.'),
    ]:
        st.markdown(f'### {title}')
        st.caption(description)
        for col, ligand, label in zip(st.columns(2), ['vegf', 'plgf'], ['Free VEGF-A', 'Free PlGF']):
            with col:
                st.markdown(f'**{label}**')
                stem = f'{group}_free_{ligand}'
                st.image(str(results / f'{stem}.png'), width='stretch')
                st.download_button(f'Download {label} plot (PDF)', (results / f'{stem}.pdf').read_bytes(), f'{stem}.pdf', 'application/pdf', key=f'download_{stem}')
        st.markdown('**Total sFLT1**')
        stem = f'{group}_total_sflt1'
        st.image(str(results / f'{stem}.png'), width='stretch')
        st.download_button('Download total sFLT1 plot (PDF)', (results / f'{stem}.pdf').read_bytes(), f'{stem}.pdf', 'application/pdf', key=f'download_{stem}')
        if group != 'degrader_half_life':
            st.caption('Total sFLT1 curves overlap: with identical natural sFLT1/complex clearance and equal degrader activity on free/bound sFLT1, these ligand parameter changes do not alter total sFLT1 turnover.')
    st.caption('Axes use different vertical scales. Storage and other-binding pools are exploratory assumptions. Peaks can occur before day 14; endpoint and peak results are both reported in the summary.')
    st.download_button('Download all nine plots (PDF)', (results / 'sensitivity_analysis_all_plots.pdf').read_bytes(), 'sensitivity_analysis_all_plots.pdf', 'application/pdf')
    st.download_button('Download sensitivity summary (CSV)', (results / 'sensitivity_summary.csv').read_bytes(), 'sensitivity_summary.csv', 'text/csv')
    st.download_button('Download sensitivity time courses (CSV)', (results / 'sensitivity_time_courses.csv').read_bytes(), 'sensitivity_time_courses.csv', 'text/csv')
    st.download_button('Download sensitivity parameters (JSON)', (results / 'sensitivity_parameters.json').read_bytes(), 'sensitivity_parameters.json', 'application/json')

with model_tab:
    with st.sidebar:
        st.header('Simulation')
        st.button('Revert to defaults',on_click=reset_defaults,use_container_width=True)
        st.button('Regenerate model',type='primary',use_container_width=True)
        st.caption('Changes update automatically. Regenerate reruns the model using the current inputs.')
        number('days','Simulation duration (days)',0.01)
        number('start_day','Treatment start (day)')
        number('stop_day','Treatment stop (day)')
        with st.expander('Starting total concentrations',expanded=True):
            number('sflt1_ng_ml','Total sFLT1 (ng/mL)')
            number('vegf_pg_ml','Total plasma VEGF-A, excluding platelets (pg/mL)')
            number('plgf_pg_ml','Total plasma PlGF, excluding placental storage (pg/mL)')
        with st.expander('Natural clearance half-lives',expanded=True):
            number('sflt1_half_hours','Free sFLT1 half-life (hours)',0.001)
            number('vegf_half_hours','Free VEGF-A half-life (hours)',0.001)
            number('plgf_half_hours','Free PlGF half-life (hours)',0.001)
            number('sv_half_hours','sFLT1–VEGF complex half-life (hours)',0.001)
            number('sp_half_hours','sFLT1–PlGF complex half-life (hours)',0.001)
        with st.expander('Degrader',expanded=True):
            number('degrader_half_days','Degrader half-life (days; single dose)',0.001)
            st.caption('One dose at treatment start. Activity decays exponentially; there is no daily redosing. The 3-day default is an adjustable reference assumption.')
            number('degrader_per_day','Initial degrader-mediated free sFLT1 clearance (/day)')
            values['complex_activity']=st.slider('Relative activity on bound sFLT1',0.,1.,1.,0.05,key='complex_activity')
            values['co_v']=st.slider('VEGF co-degradation fraction',0.,1.,float(base.co_v),0.05,key='co_v')
            values['co_p']=st.slider('PlGF co-degradation fraction',0.,1.,float(base.co_p),0.05,key='co_p')
        with st.expander('Other VEGF-binding molecules'):
            st.caption('Default capacity gives equal starting VEGF amounts in sFLT1-bound, other-bound and platelet pools, plus a small free pool. All three affinities default to 20 pM. Capacity zero disables other binding.')
            number('other_capacity_pm','Other binder total site capacity (pM)')
            number('other_kd_pm','Other binder–VEGF Kd (pM)',0.001)
            number('other_kon_pm_day','Other binder association (pM⁻¹ day⁻¹)',0.001)
            number('other_complex_half_hours','Other binder–VEGF complex half-life (hours)',0.001)
        with st.expander('Platelet VEGF storage'):
            st.caption('Stored VEGF is inaccessible to sFLT1 and the degrader. Released VEGF enters the free plasma pool. Storage is expressed per mL of modeled plasma, not per mL of platelet volume.')
            number('platelet_vegf_pg_ml','Starting platelet-stored VEGF (pg/mL equivalent)')
            number('platelet_kd_pm','Effective platelet storage Kd (pM)',0.001)
            number('platelet_kon_pm_day','Effective platelet association (pM⁻¹ day⁻¹)')
            number('platelet_release_per_day','Additional platelet release rate (/day)')
            number('platelet_half_hours','Stored VEGF loss half-life, excluding release (hours)',0.001)
            number('platelet_input_multiplier','Platelet VEGF replenishment multiplier')
            st.caption('Platelet affinity is a phenomenological storage assumption, not a measured receptor Kd. Effective site capacity is inferred from the starting free and stored amounts. Set association to zero to disable reversible exchange. Replenishment balances starting losses; 168 hours is a placeholder loss half-life.')
        with st.expander('Other PlGF-binding molecules'):
            st.caption('Independent effective binding sites for PlGF; default Kd is 150 pM, matching PlGF–sFLT1. These sites are separate from the VEGF other-binder pool.')
            number('plgf_other_capacity_pm','PlGF other-binder total site capacity (pM)')
            number('plgf_other_kd_pm','Other binder–PlGF Kd (pM)',0.001)
            number('plgf_other_kon_pm_day','PlGF other-binder association (pM⁻¹ day⁻¹)',0.001)
            number('plgf_other_complex_half_hours','Other binder–PlGF complex half-life (hours)',0.001)
        with st.expander('Placental PlGF storage'):
            st.caption('Hypothetical effective storage pool using the same equations as VEGF. Default effective Kd is 150 pM. This is a scenario assumption, not a measured placental PlGF reservoir or affinity.')
            number('placental_plgf_pg_ml','Starting placental-stored PlGF (pg/mL equivalent)')
            number('plgf_placental_kd_pm','Effective placental PlGF storage Kd (pM)',0.001)
            number('plgf_placental_kon_pm_day','Effective placental PlGF association (pM⁻¹ day⁻¹)')
            number('plgf_placental_release_per_day','Additional placental PlGF release rate (/day)')
            number('plgf_placental_half_hours','Stored PlGF loss half-life, excluding release (hours)',0.001)
            number('plgf_placental_input_multiplier','Placental PlGF replenishment multiplier')
        with st.expander('Binding kinetics'):
            number('kd_v_pm','VEGF–sFLT1 Kd (pM)',0.001)
            number('kd_p_pm','PlGF–sFLT1 Kd (pM)',0.001)
            number('kon_v_pm_day','VEGF association rate (pM⁻¹ day⁻¹)',0.001)
            number('kon_p_pm_day','PlGF association rate (pM⁻¹ day⁻¹)',0.001)
        with st.expander('Molecular weights'):
            number('sflt1_mw_kda','sFLT1 molecular weight (kDa)',0.001)
            number('vegf_mw_kda','VEGF molecular weight (kDa)',0.001)
            number('plgf_mw_kda','PlGF molecular weight (kDa)',0.001)
        with st.expander('Production inputs',expanded=True):
            mode=st.selectbox('Production calibration',['Recalibrate to current baseline','Hold reference production fixed','Enter production directly'],key='production_mode')
            st.caption('Reference = current default parameters. All modes use your selected starting totals. Fixed or manual inputs can cause untreated baseline drift.')
            manual=[]
            if mode=='Enter production directly':
                for label,q in zip(['sFLT1','VEGF','PlGF'],setup(base)[1]):
                    manual.append(st.number_input(f'{label} input (pM/day)',min_value=0.,value=float(q),format='%.6g',key=f'manual_{label}'))
            number('input_s_multiplier','sFLT1 input multiplier')
            number('input_v_multiplier','VEGF input multiplier')
            number('input_p_multiplier','PlGF input multiplier')

    try:
        p=Parameters(**values)
        production=None if mode.startswith('Recalibrate') else setup(base)[1] if mode.startswith('Hold') else manual
        df,meta=simulate(p,production)
        control,_=simulate(p,production,treated=False)
    except (ValueError,RuntimeError) as exc:
        st.error(str(exc)); st.stop()

    st.write('Treatment clears free sFLT1 and, according to relative activity, ligand-bound sFLT1. Bound ligand is either co-degraded or returned to the free pool. Natural complex clearance removes both partners. Other-bound VEGF/PlGF, platelet-stored VEGF and placental-stored PlGF are not directly targeted by the degrader; they become accessible to sFLT1 only after dissociation or release.')
    cols=st.columns(3)
    for col,key,label in zip(cols,['total_sflt1_ng_mL','free_vegf_pg_mL','free_plgf_pg_mL'],['Total sFLT1 (ng/mL)','Free VEGF-A (pg/mL)','Free PlGF (pg/mL)']):
        initial=df[key].iloc[0]; final=df[key].iloc[-1]
        delta=f'{(final/initial-1)*100:+.2f}% from start' if initial else 'Zero starting concentration'
        col.metric(label,f'{final:.4g}',delta,delta_color='off')
    st.caption('Default all-pools totals are 100 pg/mL for each ligand. VEGF: 6.48 free and 31.17 in each nonfree pool. PlGF: 33.73 free and 22.09 in each nonfree pool. These are model-derived allocations, not measured physiological fractions.')
    st.caption('Metrics show the final simulation day; charts compare treated and untreated trajectories.')
    st.subheader('sFLT1')
    plot=pd.DataFrame({
        'Treated: free':df.free_sflt1_ng_mL.values,
        'Treated: total':df.total_sflt1_ng_mL.values,
        'Untreated: free':control.free_sflt1_ng_mL.values,
        'Untreated: total':control.total_sflt1_ng_mL.values,
    },index=df.day)
    chart(plot,'sFLT1 (ng/mL)')
    for ligand,label in [('vegf','VEGF-A'),('plgf','PlGF')]:
        st.subheader(label)
        for col,kind,title in zip(st.columns(2),['free','total'],['Free plasma','Total plasma']):
            with col:
                st.markdown(f'**{title} {label}**')
                key=f'{kind}_{ligand}_pg_mL'
                plot=pd.DataFrame({'Treated':df[key].values,'Untreated':control[key].values},index=df.day)
                chart(plot,f'{title} {label} (pg/mL)',zero=kind!='free')
                if kind=='free':
                    st.caption('Vertical axis follows the concentration range to show small changes.')
                else:
                    st.caption('Includes free, sFLT1-bound and other-bound ligand; excludes platelet storage for VEGF and placental storage for PlGF.')
    for ligand,label in [('vegf','VEGF'),('plgf','PlGF')]:
        with st.expander(f'{label} distribution across pools',expanded=True):
            storage = 'platelet' if ligand=='vegf' else 'placental'
            storage_label = 'Platelet' if ligand=='vegf' else 'Placental'
            names={f'free_{ligand}_pg_mL':f'Free plasma {label}',f'bound_{ligand}_pg_mL':f'sFLT1-bound {label}',
                   f'other_bound_{ligand}_pg_mL':f'Other-bound {label}',f'{storage}_{ligand}_pg_mL':f'{storage_label}-stored {label}',
                   f'total_{ligand}_pg_mL':'Plasma total (excludes storage)',f'all_pools_{ligand}_pg_mL':'All pools (includes storage)'}
            chart(df.set_index('day')[list(names)].rename(columns=names),f'{label} (pg/mL plasma equivalent)')
            st.caption('Plasma total = free + sFLT1-bound + other-bound ligand. All pools adds platelet VEGF or placental PlGF storage. Neither sum predicts an assay result. Other-binder and storage capacities are independent for VEGF and PlGF; the ligands compete only for sFLT1.')
    st.caption('VEGF release during sample clotting is documented (https://pubmed.ncbi.nlm.nih.gov/9640345/). The placental PlGF storage pool here is an exploratory assumption and is not validated by that VEGF study.')
    with st.expander('Bound fractions and cumulative co-degradation'):
        fractions=pd.DataFrame(index=df.day)
        for ligand in ['vegf','plgf']:
            total=df[f'total_{ligand}_pg_mL'].to_numpy()
            fractions[ligand.upper()+' sFLT1-bound (% of plasma total)']=np.divide(df[f'bound_{ligand}_pg_mL'].to_numpy(),total,out=np.zeros_like(total),where=total>0)*100
            fractions[ligand.upper()+' other-bound (% of plasma total)']=np.divide(df[f'other_bound_{ligand}_pg_mL'].to_numpy(),total,out=np.zeros_like(total),where=total>0)*100
        chart(fractions,'% bound (zero when total is zero)')
        chart(df.set_index('day')[['cumulative_vegf_codegraded_pg_mL','cumulative_plgf_codegraded_pg_mL']].rename(columns={'cumulative_vegf_codegraded_pg_mL':'VEGF co-degraded','cumulative_plgf_codegraded_pg_mL':'PlGF co-degraded'}),'Cumulative loss (pg/mL)')
    with st.expander('Effective production inputs and parameters'):
        st.dataframe(pd.DataFrame({'Protein':['sFLT1','VEGF-A','PlGF'],'Input (pM/day)':[meta['production_s_pM_day'],meta['production_v_pM_day'],meta['production_p_pM_day']],'Input (ng/mL/day for sFLT1; pg/mL/day for ligands)':[meta['production_s_pM_day']*p.sflt1_mw_kda/1000,meta['production_v_pM_day']*p.vegf_mw_kda,meta['production_p_pM_day']*p.plgf_mw_kda]}),hide_index=True)
        st.caption(f"Platelet VEGF replenishment: {meta['platelet_production_pM_day']:.6g} pM/day ({meta['platelet_production_pM_day']*p.vegf_mw_kda:.6g} pg/mL/day equivalent), separate from the free VEGF input above.")
        st.caption(f"Placental PlGF replenishment: {meta['plgf_placental_production_pM_day']:.6g} pM/day ({meta['plgf_placental_production_pM_day']*p.plgf_mw_kda:.6g} pg/mL/day equivalent), separate from free PlGF input.")
        st.json(meta)
    with st.expander('Model equations and assumptions'):
        st.markdown('''Core states are free sFLT1 (S), free VEGF (V), free PlGF (P), and complexes Cᵥ and Cₚ, in pM. Binding competes for one effective sFLT1 site with 1:1 stoichiometry. Molecular weights convert mass concentrations into molar concentrations; dimerization and isoforms are not explicit.

    For ligand L, Rₗ = konₗ S L − koffₗ Cₗ, with koffₗ = Kdₗ konₗ.

    - dS/dt = qS − (kS + d)S − Rᵥ − Rₚ
    - dL/dt = qL − kL L − Rₗ + (1 − fₗ) dC Cₗ
    - dCₗ/dt = Rₗ − (kCₗ + dC)Cₗ

    Here d(t) = d₀ exp[−ln(2)(t − t_start)/half-life] after the single dose at treatment start, within the treatment window; no repeated doses are given. The default half-life is 3 days. dC = d × relative complex activity, and f is the co-degradation fraction. Natural complex clearance destroys both partners. Half-life converts to rate as ln(2)/half-life.

    Initialization solves the **kinetic untreated steady state**, including complex clearance: Cₗ = S L / (Kdₗ + kCₗ/konₗ). Production is then calibrated to balance losses. This avoids drift caused by initializing at binding equilibrium while simultaneously clearing complexes. Changing half-life in recalibration mode also changes production.

    For each ligand L (VEGF or PlGF), two additional states are other-bound ligand O and stored ligand T (platelet VEGF or placental PlGF). Each has independent site capacities and kinetics. Other-binding flux is Rₒ = konₒ (B − O)L − koffₒ O; dO/dt = Rₒ − kₒ O. Clearance destroys bound ligand, with total binding capacity maintained externally.

    Storage exchange is Rₜ = konₜ (Bₜ − T)L − konₜ Kdₜ T. Effective capacity Bₜ = T₀(1 + Kdₜ/L₀) gives zero initial net exchange when starting free ligand is positive. Storage obeys dT/dt = qT − (r + kT)T + Rₜ. The free ligand equation gains −Rₒ − Rₜ + rT. The independently calibrated storage replenishment qT maintains the starting store, scaled by its multiplier. A zero free starting concentration uses Bₜ = T₀ and replenishment also balances initial dissociation. These effective storage models do not imply measured storage receptor affinities.

    Stored ligand is inaccessible to direct sFLT1 binding or degrader clearance until released. VEGF and PlGF compete for sFLT1; the other-binder and storage sites are separate effective pools for each ligand, with no cross-competition. PlGF storage is a hypothetical scenario extension, not an experimentally calibrated pool.

    For each ligand, starting plasma total includes free + sFLT1-bound + other-bound; platelet VEGF or placental PlGF storage is additional. Initialization solves competitive sFLT1 occupancy with both other-binding pools at kinetic steady state. Free-protein production balances plasma clearance minus additional storage release. Negative calibrated production produces an error. Fixed-reference and manual inputs override free-protein production only; each storage replenishment multiplier is separate.

    The PDF supplies illustrative starting concentrations and clearance/affinity assumptions. Association rates of 1 pM⁻¹ day⁻¹ and molecular weights are explicit modeling assumptions. No drug concentration, dose-to-clearance relationship, tissue compartments, receptor signaling, feedback or clinical outcomes are inferred.''')
    st.download_button('Download treated time course',df.to_csv(index=False),'treated.csv','text/csv')
    st.download_button('Download untreated time course',control.to_csv(index=False),'untreated.csv','text/csv')
    st.download_button('Download parameters',json.dumps({'parameters':values,'derived_metadata':meta,'production_mode':mode,'effective_production_pM_day':[meta['production_s_pM_day'],meta['production_v_pM_day'],meta['production_p_pM_day']]},indent=2),'parameters.json','application/json')
