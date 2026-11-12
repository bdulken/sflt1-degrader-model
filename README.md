# Adjustable sFLT1 / VEGF-A / PlGF model

Run locally from this folder:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/streamlit run app.py
```

Open the localhost URL printed by Streamlit. Inputs update automatically. **Regenerate model** reruns current inputs; **Revert to defaults** restores all settings. The **Sensitivity analysis** tab displays the nine saved sensitivity plots with PDF/CSV downloads; these reference results do not track sidebar edits. The **Interactive model** tab contains simulation results; the **Model schematic** tab displays the default reference diagram and offers a PDF download. The schematic does not track edited inputs. Charts show persistent lines/points and do not use scroll zoom. CSV and parameter downloads are available.

## Current defaults

Both ligands start at **100 pg/mL across all pools**, including platelet storage. The three nonfree pools share the nonfree amount equally. Joint kinetic steady-state calibration includes competition between VEGF and PlGF for sFLT1:

| Starting pool (pg/mL) | VEGF-A | PlGF |
| --- | ---: | ---: |
| Free plasma | 6.4798 | 33.7267 |
| sFLT1-bound | 31.1734 | 22.0911 |
| Other-bound | 31.1734 | 22.0911 |
| Platelet-stored | 31.1734 | 22.0911 |
| Plasma total, excluding platelets | 68.8266 | 77.9089 |

sFLT1 starts at 10 ng/mL (100 pM at 100 kDa). VEGF binding affinities default to 20 pM for all three pools; PlGF affinities default to 150 pM. Platelet affinity is an **effective storage assumption**, not a measured receptor Kd. Other-binder capacities are approximately 99.2636 pM for VEGF and 99.3073 pM for PlGF. Equal allocation is established by these defaults; parameter edits do not enforce equal thirds.

Other defaults: free sFLT1 and sFLT1–ligand complexes have 34-hour natural half-lives; free VEGF has a 34-minute half-life; free PlGF clearance is 1/day (16.64-hour half-life). Other-bound complexes default to 34-hour half-lives and stored ligand loss to 168 hours. Association rates are 1 pM^-1 day^-1. Additional platelet release defaults to zero. A single degrader dose is given on day 0; initial clearance is 2/day for free and bound sFLT1 and decays exponentially with a default 3-day half-life, with 25% ligand co-degradation (75% returned to the free pool), during days 0–14. These are scenario assumptions with different evidence bases, not a matched patient dataset.

## Model

Nine concentration states are free sFLT1, free VEGF/PlGF, two sFLT1 complexes, two other-bound pools and two platelet pools. Two additional integral states accumulate degrader-mediated ligand co-degradation. Time is days, concentration is pM. Molecular weights are 100, 45 and 30 kDa for sFLT1, VEGF and PlGF. pM = pg/mL / molecular weight in kDa. One effective binding site and 1:1 ligand stoichiometry are assumed; dimers and isoforms are not explicit.

VEGF and PlGF **compete for sFLT1**. Other-binder and platelet sites are independent effective pools for each ligand, without cross-competition. No named other-binding protein or measured platelet PlGF reservoir is inferred. The PlGF storage extension is hypothetical and independently adjustable.

For either ligand L:

- sFLT1 binding flux: R = kon*S*L − koff*C, with koff = Kd*kon.
- Other binding: RO = konO*(BO−O)*L − koffO*O; dO/dt = RO − kO*O.
- Platelet exchange: RT = konT*(BT−T)*L − konT*KdT*T.
- Platelet storage: dT/dt = qT − (r+kT)*T + RT.
- Free ligand: dL/dt = qL − kL*L − R + (1−f)*dC*C − RO − RT + r*T.
- sFLT1 complex: dC/dt = R − (kC+dC)*C.
- Free sFLT1: dS/dt = qS − (kS+d)*S − RV − RP.

Here d(t) = d0*exp(-ln(2)*(t−start)/degrader_half_days) after the single dose at treatment start, within the treatment window and dC = d times relative complex activity. f is the co-degradation fraction. Natural complex clearance destroys both partners; the ligand-release fraction applies only to degrader-mediated sFLT1 removal. Other-binding capacity is maintained externally as complexes clear. Other-bound and platelet-stored ligands are not directly targeted by the degrader, and cannot bind sFLT1 until released.

Effective platelet site capacity BT = T_initial*(1+KdT/L_initial) is inferred to give zero net initial reversible exchange. Changing Kd while retaining starting storage recalibrates capacity. Setting platelet association to zero disables reversible exchange, retaining additional release/replenishment. Setting the replenishment multiplier to zero removes new stored-ligand input. Stored concentrations are expressed per mL of modeled plasma equivalent, not within platelet volume. For zero starting free ligand, BT=T_initial and replenishment balances initial dissociation as well as loss.

## Calibration and inputs

Initialization solves kinetic untreated steady state including natural complex clearance: C = S*L/(Kd+kC/kon), with both other-bound pools included in the selected plasma totals. Production balances each free and bound clearance pathway minus additional platelet release. A negative required calibrated input is rejected. Half-life converts to clearance using ln(2)/half-life.

- **Recalibrate to current baseline** maintains selected starting totals at untreated steady state when input multipliers equal one. Changing clearance also changes production.
- **Hold reference production fixed** uses the current default scenario's free-protein inputs. Starting states can drift without treatment.
- **Enter production directly** specifies three free-protein inputs in pM/day. Starting states can drift without treatment.

The three free-protein production modes do not override platelet replenishment, which is independently calibrated and multiplied for each ligand. All input multipliers apply throughout the run to treated and untreated scenarios. Both simulations use identical initial states and inputs. Final metrics compare against starting values. Treatment discontinuities are integrated separately to prevent numerical step-over.

`total_vegf_pg_mL` and `total_plgf_pg_mL` are plasma totals excluding platelets. `bound_*` is specifically sFLT1-bound; `other_bound_*` and `platelet_*` are separate. `all_pools_*` includes storage. Bound percentages use plasma totals as denominators. Effective platelet capacities and replenishment inputs are displayed and included in parameter metadata.

## Interpretation and verification

The original project-context PDF contains an earlier conversation, not a validated parameter dataset. The PDF and its extracted transcript are kept locally and excluded from this repository. `model.py` reconstructs the previously missing source; earlier outputs may differ because initialization, pools and defaults have changed. This is an exploratory mechanistic model, not a fitted PK/PD or clinical outcome predictor. Drug PK, tissue exchange, gestational changes, receptor signaling and feedback are omitted.

Neither total is an assay prediction. VEGF release during clotting can increase serum VEGF: [Webb et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9640345/) and [Banks et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9528841/). These VEGF studies do not validate the analogous PlGF pool. Assay recovery and sample processing require separate calibration.

Optional development tools for the schematic/PDF and sensitivity scripts require `.venv/bin/pip install -r requirements-dev.txt`.

Run checks and the existing sensitivity script:

```sh
.venv/bin/python -m unittest discover -s tests
MPLBACKEND=Agg MPLCONFIGDIR=.mplconfig .venv/bin/python vegf_half_life_sensitivity.py
```

The sensitivity script recalibrates production for each VEGF half-life and uses current model defaults. Tests cover untreated steady state, analytic sFLT1 turnover, treatment windows, zero inputs, independent pools, allocation, ligand conservation, and invalid calibration.


## Deploy on Streamlit Community Cloud

1. Push this repository to GitHub on branch `main`.
2. Sign in at https://share.streamlit.io/ with the GitHub account that can access the repository.
3. Select **Create app** and choose the repository, branch `main`, and main file `app.py`.
4. In **Advanced settings**, choose **Python 3.12** and deploy.

The app requires no secrets or external data files. Runtime dependencies are pinned in `requirements.txt`; plotting and PDF development tools are separate. A private GitHub repository can be used with the required GitHub access granted to Streamlit. Review app sharing settings separately from repository visibility. See the [official deployment guide](https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy).

GitHub Actions runs model tests and a dashboard smoke check. The context PDF/transcript, virtual environment, secrets and generated sensitivity outputs are excluded. The schematic and its generation script are included.


## Sensitivity analyses

Run `MPLCONFIGDIR=.mplconfig .venv/bin/python sensitivity_analysis.py` after installing `requirements-dev.txt`. Outputs are saved to `sensitivity_results/`: nine individual PNG/PDF plots (including total sFLT1), one nine-page PDF with baseline-default tables, full time-course/summary CSVs and scenario parameters. Analyses vary one input at a time: free VEGF half-life (1, 4, 8, 12 hours), single-dose degrader half-life (1, 3, 7 days), or joint VEGF/PlGF co-degradation (25%, 50%, 75%). Each scenario recalibrates production to initial untreated steady state. All analyses use a single dose at day 0 with initial clearance 2/day and exponential activity decay, with no redosing. Degrader half-life is 3 days except when varied over 1, 3, and 7 days. The dashboard exposes the degrader half-life. `simulate()` uses the parameter half-life by default; an explicit `degrader_half_days=None` is available only as a constant-exposure counterfactual. A positive override selects another single-dose half-life. The 3-day reference is an adjustable assumption.
