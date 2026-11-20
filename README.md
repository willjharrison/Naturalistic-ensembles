# Naturalistic Ensemble & Spatial Structure Object Recognition

Codebase for the manuscript *"Do visual ensemble statistics support object and texture classification?"* (Gu, Rideaux, Wallis, Bex, & Harrison).

## Repository Structure
- `experiments/`: MATLAB/Psychtoolbox stimulus generation and experiment scripts (`objectRecog_v03.m`, `objectRecog_demos_v01.m`, `exp_stim_demo.m`).
- `analysis/`: Behavioural analysis scripts (`.m`), summary data (`.csv`), and JASP Bayesian ANOVA/t-test files (`.jasp`) for Experiments 1–3 (`exp1/`, `exp2/`, `exp3/`).
- `models/`:
  - `sc_model/`: Original steered-filter spatial/orientation observer models (`SCModel_v01.m`–`v05`, `kFoldSavedData_v01_parfor.m`, `modelSummary_v01.m`, `simData.m`).
  - `summary_stats_model/`: Spatially invariant circular mean/SD observer model (`sumStatsModel_v01.m`, `qStats.m`, `ensembleModel.csv`).
  - `clean_comparison/`: Refactored single-pass Experiment 1 comparison pipeline (`SCModel_clean_compare.m`, `compare_models.m`, `cleanModelSummary.csv`).

## External Dependencies (Not Bundled)
To avoid duplicating existing toolboxes, the following dependencies must be on your MATLAB path:
- **Psychophysics Toolbox (PTB-3)** — required for running experiment scripts and `BalanceFactors` / `GetSecs`.
- **CircStat2012a** (`circ_mean`, `circ_std`, `circ_dist`, `wrap`) — available in [willjharrison/ImageProcessingFoundationsWorkshop](https://github.com/willjharrison/ImageProcessingFoundationsWorkshop).
- **Custom image-processing & helper functions** (`polarDistFun`, `createLowPassFilter`, etc.) — available in the `helpers/` folder of [willjharrison/ImageProcessingFoundationsWorkshop](https://github.com/willjharrison/ImageProcessingFoundationsWorkshop).
- **MATLAB Toolboxes**: Statistics and Machine Learning Toolbox (`fitglme`, `fitglm`, `grpstats`), Image Processing Toolbox (`imgaussfilt`), and Parallel Computing Toolbox (`parfor`).

## Data Note
Raw `.mat` trial outputs (`Output/`), precomputed Gabor/filter maps (`modIms/`), and source image databases (THINGS, DTD) are excluded via `.gitignore`.
