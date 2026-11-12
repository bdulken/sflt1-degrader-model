"""Concise reference table generated from model defaults."""
import pandas as pd
from model import Parameters


def default_parameter_table():
    p=Parameters()
    return pd.DataFrame({
        'Parameter':['Starting total concentration','Binding Kd to sFLT1','Ligand co-degradation'],
        'sFLT1':[f'{p.sflt1_ng_ml:g} ng/mL','—','— (target removed)'],
        'VEGF-A':[f'{p.vegf_pg_ml+p.platelet_vegf_pg_ml:.6g} pg/mL',f'{p.kd_v_pm:g} pM',f'{p.co_v:.0%}'],
        'PlGF':[f'{p.plgf_pg_ml+p.placental_plgf_pg_ml:.6g} pg/mL',f'{p.kd_p_pm:g} pM',f'{p.co_p:.0%}'],
    })
