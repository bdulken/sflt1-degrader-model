"""Compact reference table generated from the actual model defaults."""
from dataclasses import asdict
import pandas as pd
from model import Parameters, simulate


def default_parameter_table():
    p=Parameters()
    values=asdict(p)
    covered=set()
    rows=[]

    def row(label,unit='',s=None,v=None,pl=None,g=None):
        def value(key):
            if key is None: return '—'
            covered.add(key)
            return f'{values[key]:.6g}'
        rows.append({'Parameter':label,'Unit':unit,'sFLT1':value(s),'VEGF-A':value(v),'PlGF':value(pl),'Shared default':value(g)})

    row('Simulation duration','days',g='days')
    row('Single dose: administration day','day',g='start_day')
    row('Exposure window ends','day',g='stop_day')
    row('Degrader half-life','days',g='degrader_half_days')
    row('Initial degrader clearance','day⁻¹',g='degrader_per_day')
    row('Relative activity on bound sFLT1','fraction',g='complex_activity')
    row('Ligand co-degradation','fraction',v='co_v',pl='co_p')
    row('Starting total: sFLT1 / plasma ligand','ng/mL sFLT1; pg/mL ligands',s='sflt1_ng_ml',v='vegf_pg_ml',pl='plgf_pg_ml')
    row('Molecular weight','kDa',s='sflt1_mw_kda',v='vegf_mw_kda',pl='plgf_mw_kda')
    row('Free-protein half-life','hours',s='sflt1_half_hours',v='vegf_half_hours',pl='plgf_half_hours')
    row('sFLT1–ligand complex half-life','hours',v='sv_half_hours',pl='sp_half_hours')
    row('sFLT1 binding Kd','pM',v='kd_v_pm',pl='kd_p_pm')
    row('sFLT1 association rate','pM⁻¹ day⁻¹',v='kon_v_pm_day',pl='kon_p_pm_day')
    row('Free-protein input multiplier','×',s='input_s_multiplier',v='input_v_multiplier',pl='input_p_multiplier')
    row('Other-binder site capacity','pM',v='other_capacity_pm',pl='plgf_other_capacity_pm')
    row('Other-binder Kd','pM',v='other_kd_pm',pl='plgf_other_kd_pm')
    row('Other-binder association rate','pM⁻¹ day⁻¹',v='other_kon_pm_day',pl='plgf_other_kon_pm_day')
    row('Other-bound complex half-life','hours',v='other_complex_half_hours',pl='plgf_other_complex_half_hours')
    row('Starting stored ligand','pg/mL plasma equivalent',v='platelet_vegf_pg_ml',pl='placental_plgf_pg_ml')
    row('Effective storage Kd','pM',v='platelet_kd_pm',pl='plgf_placental_kd_pm')
    row('Storage association rate','pM⁻¹ day⁻¹',v='platelet_kon_pm_day',pl='plgf_placental_kon_pm_day')
    row('Additional storage release rate','day⁻¹',v='platelet_release_per_day',pl='plgf_placental_release_per_day')
    row('Stored-ligand loss half-life','hours',v='platelet_half_hours',pl='plgf_placental_half_hours')
    row('Storage replenishment multiplier','×',v='platelet_input_multiplier',pl='plgf_placental_input_multiplier')
    if covered != set(values):
        raise ValueError(f'Default table coverage mismatch: {set(values)-covered}')
    # Derived concentrations and calibrated inputs help interpret the independent parameters.
    df,meta=simulate(p,treated=False)
    initial=df.iloc[0]
    def derived(label,unit,s='—',v='—',pl='—'):
        rows.append({'Parameter':label,'Unit':unit,'sFLT1':s,'VEGF-A':v,'PlGF':pl,'Shared default':'—'})
    derived('Initial free concentration (derived)','ng/mL sFLT1; pg/mL ligands',f'{initial.free_sflt1_ng_mL:.4g}',f'{initial.free_vegf_pg_mL:.4g}',f'{initial.free_plgf_pg_mL:.4g}')
    derived('Initial all-pools ligand total (derived)','pg/mL',v=f'{initial.all_pools_vegf_pg_mL:.4g}',pl=f'{initial.all_pools_plgf_pg_mL:.4g}')
    derived('Free-protein production (calibrated)','ng/mL/day sFLT1; pg/mL/day ligands',f'{meta["production_s_pM_day"]*p.sflt1_mw_kda/1000:.4g}',f'{meta["production_v_pM_day"]*p.vegf_mw_kda:.4g}',f'{meta["production_p_pM_day"]*p.plgf_mw_kda:.4g}')
    derived('Storage replenishment (calibrated)','pg/mL/day plasma equivalent',v=f'{meta["platelet_production_pM_day"]*p.vegf_mw_kda:.4g}',pl=f'{meta["plgf_placental_production_pM_day"]*p.plgf_mw_kda:.4g}')
    return pd.DataFrame(rows)
