# Sensitivity results

Nine separate plots (PNG and vector PDF), including total sFLT1 for each analysis, and a combined nine-page PDF. Each page starts with a compact baseline-default table (also exported as baseline_parameters.csv). The parameter being varied supersedes its baseline value in each scenario.

- Free VEGF half-life: 1, 4, 8, 12 hours. PlGF free half-life stays at its default (16.64 hours).
- Degrader half-life: 1, 3, 7 days, with single-dose initial effective clearance of 2/day. Effective clearance decays as 2*exp(-ln(2)*t/half-life). No redosing. This assumes activity is proportional to exposure, not a saturable PK/PD model.
- Co-degradation fraction: 25%, 50%, 75%, changed together for both ligands. The remainder returns to free ligand.

These are separate one-at-a-time analyses, not a full factorial combination. All analyses use a single dose on day 0 with initial clearance 2/day and exponential decay. Degrader half-life is 3 days in the VEGF-half-life and co-degradation analyses. Both ligand co-degradation fractions stay at 25% in the other analyses. Free VEGF half-life stays at 34 minutes except in its own analysis. Treatment begins on day 0; simulation and treatment windows end at day 14.

For each scenario, production is recalibrated to the initial untreated kinetic steady state, then kept fixed throughout the simulation. Initial distributions, other binding and platelet storage retain their current default values. Untreated curves coincide across scenarios and are shown as a single dashed baseline. Production differs across VEGF-half-life scenarios; this is not a fixed-production comparison.

Summary CSV reports ligand day-14 and peak concentrations with peak timing, plus total sFLT1 endpoints and nadirs (ng/mL); decaying degrader responses can peak before day 14. Full treated/untreated time courses and effective clearance are in the time-course CSV. Full defaults and effective inputs are in the parameters JSON. Units for concentrations are pg/mL and time is days. All values are exploratory model predictions, not validated physiological estimates.

Regenerate from the project folder: `MPLCONFIGDIR=.mplconfig .venv/bin/python sensitivity_analysis.py` (requires requirements-dev.txt).
