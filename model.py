"""Circulating sFLT1 model with optional alternative-binding and platelet pools for VEGF and PlGF. Concentrations: pM; time: days."""
from dataclasses import dataclass, asdict
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

@dataclass
class Parameters:
    sflt1_ng_ml: float = 10
    vegf_pg_ml: float = 68.82660282096072
    plgf_pg_ml: float = 77.90890665843557
    sflt1_mw_kda: float = 100
    vegf_mw_kda: float = 45
    plgf_mw_kda: float = 30
    sflt1_half_hours: float = 34
    vegf_half_hours: float = 34/60
    plgf_half_hours: float = 24*np.log(2)
    sv_half_hours: float = 34
    sp_half_hours: float = 34
    kd_v_pm: float = 20
    kd_p_pm: float = 150
    kon_v_pm_day: float = 1
    kon_p_pm_day: float = 1
    degrader_half_days: float = 3
    degrader_per_day: float = 2
    complex_activity: float = 1
    co_v: float = 0.25
    co_p: float = 0.25
    start_day: float = 0
    stop_day: float = 14
    days: float = 14
    input_s_multiplier: float = 1
    input_v_multiplier: float = 1
    input_p_multiplier: float = 1

    other_capacity_pm: float = 99.26363022194786
    other_kd_pm: float = 20
    other_kon_pm_day: float = 1
    other_complex_half_hours: float = 34
    platelet_vegf_pg_ml: float = 31.17339717903927
    platelet_kd_pm: float = 20
    platelet_kon_pm_day: float = 1
    platelet_release_per_day: float = 0
    platelet_half_hours: float = 168
    platelet_input_multiplier: float = 1

    plgf_other_capacity_pm: float = 99.30725784046581
    plgf_other_kd_pm: float = 150
    plgf_other_kon_pm_day: float = 1
    plgf_other_complex_half_hours: float = 34
    platelet_plgf_pg_ml: float = 22.091093341564427
    plgf_platelet_kd_pm: float = 150
    plgf_platelet_kon_pm_day: float = 1
    plgf_platelet_release_per_day: float = 0
    plgf_platelet_half_hours: float = 168
    plgf_platelet_input_multiplier: float = 1

    def validate(self):
        d=asdict(self)
        if not all(np.isfinite(v) for v in d.values()):
            raise ValueError('All parameters must be finite.')
        positive=['sflt1_mw_kda','vegf_mw_kda','plgf_mw_kda','sflt1_half_hours','vegf_half_hours','plgf_half_hours','sv_half_hours','sp_half_hours','kd_v_pm','kd_p_pm','kon_v_pm_day','kon_p_pm_day','days','other_kd_pm','other_kon_pm_day','other_complex_half_hours','platelet_half_hours','platelet_kd_pm','plgf_other_kd_pm','plgf_other_kon_pm_day','plgf_other_complex_half_hours','plgf_platelet_kd_pm','plgf_platelet_half_hours','degrader_half_days']
        if any(d[k]<=0 for k in positive) or any(v<0 for v in d.values()):
            raise ValueError('Rates/concentrations must be nonnegative; half-lives, affinities, association rates, molecular weights and duration must be positive.')
        if any(d[k]>1 for k in ['co_v','co_p','complex_activity']):
            raise ValueError('Fractions must be between zero and one.')
        if self.stop_day<self.start_day:
            raise ValueError('Treatment stop must be at or after treatment start.')

def pool_parameters(p):
    """Paired VEGF/PlGF parameters; each ligand has independent storage sites."""
    return dict(
        total=np.array([p.vegf_pg_ml/p.vegf_mw_kda,p.plgf_pg_ml/p.plgf_mw_kda]),
        mw=np.array([p.vegf_mw_kda,p.plgf_mw_kda]),
        kd=np.array([p.kd_v_pm,p.kd_p_pm]),
        kon=np.array([p.kon_v_pm_day,p.kon_p_pm_day]),
        capacity=np.array([p.other_capacity_pm,p.plgf_other_capacity_pm]),
        other_kd=np.array([p.other_kd_pm,p.plgf_other_kd_pm]),
        other_kon=np.array([p.other_kon_pm_day,p.plgf_other_kon_pm_day]),
        other_clear=np.log(2)*24/np.array([p.other_complex_half_hours,p.plgf_other_complex_half_hours]),
        stored=np.array([p.platelet_vegf_pg_ml/p.vegf_mw_kda,p.platelet_plgf_pg_ml/p.plgf_mw_kda]),
        storage_kd=np.array([p.platelet_kd_pm,p.plgf_platelet_kd_pm]),
        storage_kon=np.array([p.platelet_kon_pm_day,p.plgf_platelet_kon_pm_day]),
        release=np.array([p.platelet_release_per_day,p.plgf_platelet_release_per_day]),
        storage_clear=np.log(2)*24/np.array([p.platelet_half_hours,p.plgf_platelet_half_hours]),
        replenish=np.array([p.platelet_input_multiplier,p.plgf_platelet_input_multiplier]),
    )

def setup(p):
    p.validate()
    ks,kv,kp,kcv,kcp=np.log(2)*24/np.array([p.sflt1_half_hours,p.vegf_half_hours,p.plgf_half_hours,p.sv_half_hours,p.sp_half_hours])
    pools=pool_parameters(p)
    st=p.sflt1_ng_ml*1000/p.sflt1_mw_kda
    av=pools['kd']+np.array([kcv,kcp])/pools['kon']
    ao=pools['other_kd']+pools['other_clear']/pools['other_kon']
    def split(s):
        free=np.array([brentq(lambda f:f*(1+s/av[i])+pools['capacity'][i]*f/(ao[i]+f)-total,0,total) if total else 0. for i,total in enumerate(pools['total'])])
        complex_=s*free/av
        other=pools['capacity']*free/(ao+free)
        return free,complex_,other
    s=brentq(lambda s:s+split(s)[1].sum()-st,0,st) if st else 0.
    free,complex_,other=split(s)
    # States: S, V, P, CV, CP, OV, TV, OP, TP.
    y=np.r_[s,free,complex_,other[0],pools['stored'][0],other[1],pools['stored'][1]]
    qlig=np.array([kv,kp])*free+np.array([kcv,kcp])*complex_+pools['other_clear']*other-pools['release']*pools['stored']
    q=np.r_[ks*s+np.dot([kcv,kcp],complex_),qlig]
    return y,q,(ks,kv,kp,kcv,kcp)

_DEFAULT_EXPOSURE = object()

def simulate(p=None, production=None, treated=True, degrader_half_days=_DEFAULT_EXPOSURE):
    p=p or Parameters()
    if degrader_half_days is _DEFAULT_EXPOSURE:
        degrader_half_days=p.degrader_half_days
    if degrader_half_days is not None and (not np.isfinite(degrader_half_days) or degrader_half_days<=0):
        raise ValueError('Degrader half-life must be finite and positive, or None for constant exposure.')
    y0,q,rates=setup(p)
    if production is not None:
        q=np.asarray(production,dtype=float)
        if q.shape!=(3,) or not np.isfinite(q).all() or (q<0).any():
            raise ValueError('Production must contain three finite nonnegative pM/day values.')
    elif (q[1:] < -1e-10).any():
        raise ValueError('Platelet release exceeds the ligand loss supported by the selected starting plasma totals. Reduce platelet release/storage, increase plasma ligand totals, or select manual/fixed production; calibration would otherwise require negative VEGF production or negative PlGF production.')
    q=np.maximum(q,0.)*np.array([p.input_s_multiplier,p.input_v_multiplier,p.input_p_multiplier])
    ks,kv,kp,kcv,kcp=rates
    pools=pool_parameters(p)
    free0=y0[1:3]
    storage_capacity=pools['stored']*(1+np.divide(pools['storage_kd'],free0,out=np.zeros(2),where=free0>0))
    initial_exchange=pools['storage_kon']*((storage_capacity-pools['stored'])*free0-pools['storage_kd']*pools['stored'])
    qt=((pools['storage_clear']+pools['release'])*pools['stored']-initial_exchange)*pools['replenish']
    co=np.array([p.co_v,p.co_p])
    def rhs(t,y,d):
        if degrader_half_days is not None:
            d=d*np.exp(-np.log(2)*max(0.,t-p.start_day)/degrader_half_days)
        s=y[0]; free=y[1:3]; complexes=y[3:5]
        other=y[[5,7]]; stored=y[[6,8]]
        ro=pools['other_kon']*((pools['capacity']-other)*free-pools['other_kd']*other)
        exchange=pools['storage_kon']*((storage_capacity-stored)*free-pools['storage_kd']*stored)
        release=pools['release']*stored
        binding=pools['kon']*(s*free-pools['kd']*complexes)
        dc=d*p.complex_activity
        d_free=q[1:]-np.array([kv,kp])*free-binding+(1-co)*dc*complexes-ro+release-exchange
        d_complex=binding-(np.array([kcv,kcp])+dc)*complexes
        d_other=ro-pools['other_clear']*other
        d_stored=qt-(pools['storage_clear']+pools['release'])*stored+exchange
        return np.r_[q[0]-(ks+d)*s-binding.sum(),d_free,d_complex,d_other[0],d_stored[0],d_other[1],d_stored[1],co*dc*complexes]
    cuts=sorted(set([0.,p.days]+[x for x in [p.start_day,p.stop_day] if 0<x<p.days]))
    grid=np.unique(np.r_[np.linspace(0,p.days,701),cuts])
    state=np.r_[y0,0.,0.]; ts=[]; ys=[]
    for a,b in zip(cuts[:-1],cuts[1:]):
        d=p.degrader_per_day if treated and p.start_day <= (a+b)/2 < p.stop_day else 0.
        sol=solve_ivp(lambda t,y:rhs(t,y,d),(a,b),state,method='BDF',rtol=1e-8,atol=1e-10,dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        t=grid[(grid>=a)&(grid<=b)]
        if ts: t=t[t>a]
        ts.extend(t); ys.extend(sol.sol(t).T); state=sol.y[:,-1]
    raw=np.array(ys)
    if raw.min() < -1e-7: raise RuntimeError('Solver produced negative concentrations.')
    raw=np.maximum(raw,0)
    s=raw[:,0]; complexes=raw[:,3:5]
    columns=dict(day=ts,free_sflt1_ng_mL=s*p.sflt1_mw_kda/1000,total_sflt1_ng_mL=(s+complexes.sum(axis=1))*p.sflt1_mw_kda/1000)
    for i,name in enumerate(['vegf','plgf']):
        mw=pools['mw'][i]; free=raw[:,1+i]; bound=raw[:,3+i]; other=raw[:,5+2*i]; stored=raw[:,6+2*i]
        for label,value in [('free',free),('total',free+bound+other),('bound',bound),('other_bound',other),('platelet',stored),('all_pools',free+bound+other+stored),('cumulative',raw[:,9+i])]:
            key=f'cumulative_{name}_codegraded_pg_mL' if label=='cumulative' else f'{label}_{name}_pg_mL'
            columns[key]=value*mw
    df=pd.DataFrame(columns)
    meta={**asdict(p),'production_s_pM_day':q[0],'production_v_pM_day':q[1],'production_p_pM_day':q[2],
          'platelet_production_pM_day':qt[0],'platelet_effective_capacity_pM':storage_capacity[0],
          'degrader_half_days':degrader_half_days,'plgf_platelet_production_pM_day':qt[1],'plgf_platelet_effective_capacity_pM':storage_capacity[1]}
    return df,meta

def model(kdeg=2,kv=None,co=0.25,**kwargs):
    """Compatibility entry point for vegf_half_life_sensitivity.py."""
    if kv is not None: kwargs['vegf_half_hours']=np.log(2)*24/kv
    return simulate(Parameters(degrader_per_day=kdeg,co_v=co,co_p=co,**kwargs))
