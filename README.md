# Naturalistic Ensembles & Spatial Structure

Codebase and summary outputs for the manuscript *"Do visual ensemble statistics support object and texture classification?"* (Gu, C., Rideaux, R., Wallis, T.S.A., Bex, P.J., & Harrison, W.J.).

## Repository Structure

- `experiments/` — MATLAB/Psychophysics Toolbox stimulus-generation and psychophysics scripts (`objectRecog_v03.m`, `exp_stim_demo.m`, `objectRecog_demos_v01.m`).
- `analysis/exp1/`, `analysis/exp2/`, `analysis/exp3/` — Behavioural aggregation scripts (`prepareData.m`, `meanAccuracy.m`, `objectRecog_groupMeans.m`), participant-level CSVs (`individualData.csv`, `meanData.csv`), and JASP Bayesian ANOVA / one-sample *t*-test files.
- `models/` — Final unified computational observer modelling and working-memory simulation pipeline:
  - `SCModel_unified_allExps.m`: Extracts steerable-filter V1 oriented matched-filter (`v1Diff`), spatial energy (`enDiff`), and energy-weighted / unweighted circular orientation summary predictors (`absMuDiff`, `absSDDiff`) across all 3 experiments (34,743 valid trials, zero condition labels) and evaluates sensory noise grids (`0–0.30π × 0–4 px`).
  - `fit_per_observer_V1.m`: Fits all 6 nested observer models independently to each of the 30 participants via 5-fold out-of-sample cross-validation, computing total BIC, individual BIC winner counts, standardised `β` weights, and emergent `GroupRMSE`.
  - `simulate_ensemble_bias_V1.m`: Simulates all 30 study participants (using their individually fitted `V1_Spatial` sensitivities) on an 8-Gabor working-memory orientation-recall task (`SD = 10°, 20°, 35°`; 18,000 trials) to test whether modest spatial readout uncertainty reproduces classic ensemble biases without summary-statistic computation.
  - `qStats.m`: Fast grouping summary helper used by the modelling scripts.
  - `results/`: Final summary CSVs (`unifiedSummary_v2_V1.csv`, `perObserver_summary.csv`, `simulated_ensemble_bias_summary.csv`).

## External Dependencies

Raw `.mat` session data and large image databases are excluded from version control. To reproduce the full image-level modelling pipeline from scratch, ensure the following are on your MATLAB path:

1. **Psychophysics Toolbox** (for running `experiments/objectRecog_v03.m`).
2. **Image Processing & Circular Statistics Helpers**: `willjharrison/ImageProcessingFoundationsWorkshop` (provides `CircStat2012a`, `helpers/`, and steerable filter utilities).
3. **Stimulus Databases**: THINGS object database (Hebart et al., 2019) and Describable Textures Dataset (DTD; Cimpoi et al., 2014). Update `thingsRoot` and `dtdRoot` at the top of `models/SCModel_unified_allExps.m` to point to your local copies.

## Reproducing the Modelling & Simulation Results

From `models/` in MATLAB:

```matlab
SCModel_unified_allExps;    % Extract predictors & run sensory-noise grid across Exps 1–3
fit_per_observer_V1;        % Fit all 30 observers independently (5-fold CV, BIC, figures)
simulate_ensemble_bias_V1;  % Run N = 30 working-memory ensemble-bias simulation & figure
```
