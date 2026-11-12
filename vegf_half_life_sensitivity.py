from pathlib import Path
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from model import model
OUT=Path(__file__).resolve().parent
cases=[('34 minutes',34/60/24,'34min'),('3 hours',3/24,'3h'),('12 hours',12/24,'12h'),('24 hours',1,'24h')]
rows=[]
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,(name,half,suffix) in zip(axs.flat,cases):
    df,params=model(kdeg=2,kv=np.log(2)/half,co=1)
    df.to_csv(OUT/f'vegf_half_life_{suffix}.csv',index=False)
    control,_=model(kdeg=0,kv=np.log(2)/half,co=1)
    ax.plot(df.day,df.free_vegf_pg_mL,label='With degrader',lw=2.5)
    ax.plot(control.day,control.free_vegf_pg_mL,'--',label='No degrader',lw=1.5)
    ax.set(title=f'VEGF half-life: {name}',xlabel='Days',ylabel='Free VEGF-A (pg/mL)')
    ax.grid(alpha=.25);ax.legend()
    fig1,ax1=plt.subplots(figsize=(8,4.8),layout='constrained')
    ax1.plot(df.day,df.free_vegf_pg_mL,label='Free VEGF-A: degrader',lw=2.5)
    ax1.plot(control.day,control.free_vegf_pg_mL,'--',label='Free VEGF-A: no degrader')
    ax1.set(title=f'VEGF half-life {name} | 100% complex co-degradation',xlabel='Days',ylabel='Free VEGF-A (pg/mL)');ax1.grid(alpha=.25);ax1.legend()
    fig1.savefig(OUT/f'vegf_half_life_{suffix}.png',dpi=170);plt.close(fig1)
    for day in [0,1,3,7,14]:
        r=df.iloc[np.argmin(abs(df.day-day))]
        rows.append({'VEGF half-life':name,'day':day,'free VEGF (pg/mL)':r.free_vegf_pg_mL,'total VEGF (pg/mL)':r.total_vegf_pg_mL,'free PlGF (pg/mL)':r.free_plgf_pg_mL,'total sFLT1 (ng/mL)':r.total_sflt1_ng_mL,'VEGF production (pg/mL/day)':params['production_v_pM_day']*45})
fig.suptitle('Free VEGF-A sensitivity to assumed half-life | sFLT1 degrader 2/day',fontsize=14)
fig.savefig(OUT/'vegf_half_life_four_panels.png',dpi=180);plt.close(fig)
pd.DataFrame(rows).to_csv(OUT/'vegf_half_life_comparison.csv',index=False)
print(pd.DataFrame(rows).query('day in [0,7,14]').round(3).to_string(index=False))
