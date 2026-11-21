%% simulate_ensemble_bias_V1.m
% Directly simulates Utochkin & Brady (2020, JEP:HPP, Experiment 1, Figure 4)
% using our winning purely spatial V1 Oriented Matched-Filter observer
% parameterised by ALL 30 of our real study participants' individually fitted
% V1_Spatial sensitivities from Experiments 1-3.
%
% Version 3 (Fixed Narrow V1 Readout Window + Ring-Anchored Localisation Jitter):
%   - Directly addresses the concern that a wide spatial pooling kernel
%     (sigma_pos = 105 px) could be construed as a neural implementation of
%     ensemble coding.
%   - Primary model uses the EXACT SAME strictly local V1 readout window
%     (sigma_read = 32 px, < 1/3 of the 101.8 px inter-item spacing) for BOTH
%     Set Size 1 and Set Size 4.
%   - At Set Size 4, the ONLY spatial change is trial-by-trial angular
%     localisation / binding jitter of the readout centre along the stimulus
%     ring (theta_read ~ N(theta_cued, sigma_ang^2)), reflecting spatial
%     binding uncertainty among occupied array locations under memory load
%     without ever pooling across the full display on a single trial.
%   - Also runs two full sensitivity sweeps (saved to CSV + 2x3 figure):
%       Sweep A: Deterministic readout width sigma_pos = [32, 50, 70, 90, 110] px
%       Sweep B: Fixed local window (32 px) + ring jitter sigma_ang = [0, 35, 50, 65, 80] deg
%
% Outputs saved in Working/:
%   - fig_simulated_ensemble_bias.png / .pdf       (4-panel primary model: 32 px local + ring jitter)
%   - fig_simulated_sensitivity_sweeps.png / .pdf  (2x3 sweep across sigma_pos & sigma_ang)
%   - simulated_ensemble_bias_summary.csv
%   - simulated_sigma_sweep.csv
%   - simulated_jitter_sweep.csv

clear; close all; clc;
rng(42);

%% 1. Extract All 30 Real Participants' Individual V1_Spatial Sensitivities
fprintf('Loading cached unifiedPredictors_v2_V1.mat for N = 30 observer sensitivities...\n');
load('unifiedPredictors_v2_V1.mat', 'leftV1Mat', 'rightV1Mat', 'validTrial');

expDirs = { ...
    '../Matlab/Output', ...
    '../Matlab/Experiment 2/Output', ...
    '../Matlab/Experiment 3/Output'};

partIDVec = []; respVec = [];
for e = 1:3
    fList = dir(fullfile(expDirs{e}, '*.mat'));
    fList = fList(~ismember({fList.name}, {'.', '..', '.DS_Store'}));
    fList = fList(~contains({fList.name}, {'exclusion'}));
    for f = 1:length(fList)
        ff = load(fullfile(fList(f).folder, fList(f).name), 'task', 'dv', 'subject_name');
        nT = ff.task.nTrials;
        partIDVec = [partIDVec; repmat(e*100 + str2double(ff.subject_name), nT, 1)]; %#ok<AGROW>
        respVec   = [respVec;   ff.dv.chosenSquare(:) - 1]; %#ok<AGROW>
    end
end
partIDVec = partIDVec(validTrial);
respVec   = respVec(validTrial);

cfg = 1; % optimal (OriNoise=0, PosBlur=0)
v1Raw = rightV1Mat(validTrial, cfg) - leftV1Mat(validTrial, cfg);

uParts = unique(partIDVec);
nParts = numel(uParts);
betaV1_obs = zeros(nParts, 1);
for p = 1:nParts
    idxP = (partIDVec == uParts(p));
    Tp = table(respVec(idxP), zscore(v1Raw(idxP)), 'VariableNames', {'resp', 'v1Diff'});
    mV1 = fitglm(Tp, 'resp ~ v1Diff', 'Distribution', 'Binomial');
    betaV1_obs(p) = max(mV1.Coefficients.Estimate(2), 0.10);
end
fprintf('Extracted V1_Spatial beta for N = %d observers (mean = %.3f ± %.3f).\n', ...
    nParts, mean(betaV1_obs), std(betaV1_obs)/sqrt(nParts));

%% 2. Lower-SF Display Geometry & Common Trial Parameters
nItems     = 4;
rangesDeg  = [30, 60, 120];
nRanges    = numel(rangesDeg);
baseOffs   = {[-15, -5, 5, 15], [-30, -10, 10, 30], [-60, -20, 20, 60]};
meanDist   = [10, 20, 40];  % theoretical mean distance to display mean (deg)
nTrials    = 240;           % trials per cell per observer

% Spatial display geometry (4 corners on circle of radius R = 72 px; adjacent = 101.8 px)
imSize     = 256;
ringRadius = 64;
itemSD     = 10;            % broader envelope for lower-SF stimuli
gabFreq    = 0.055;         % lower carrier spatial frequency (cyc/px)
[X, Y]     = meshgrid(1:imSize, 1:imSize);
cx = imSize/2; cy = imSize/2;
cornerAng  = deg2rad([45, 135, 225, 315]);
locX       = cx + ringRadius * cos(cornerAng);
locY       = cy + ringRadius * sin(cornerAng);

E_maps = zeros(imSize, imSize, nItems);
for j = 1:nItems
    E_maps(:,:,j) = exp(-((X - locX(j)).^2 + (Y - locY(j)).^2) / (2 * itemSD^2));
end

phiDeg    = (0:1:359)';
phiRad    = deg2rad(phiDeg);
baseScale = 19.5;
loadScale = 0.60; % mild set-size resource scaling (leaves SS4 range-30 error near ~23 deg!)

%% 3. Baseline Set Size 1 (Fixed Local Readout sigma_read = 32 px, Zero Jitter)
sigmaReadFixed = 32; % < 1/3 of the 101.8 px inter-item spacing!
fprintf('\nRunning Set Size 1 baseline (sigma_read = %d px, sigma_ang = 0 deg)...\n', sigmaReadFixed);
[errSS1, biasSS1] = runRingSpatialSim(betaV1_obs, baseScale, 1.0, cx, cy, ringRadius, cornerAng, ...
    locX, locY, itemSD, sigmaReadFixed, 0, baseOffs, nTrials, phiDeg, phiRad, 1001);
mErr1  = mean(errSS1, 1);   seErr1  = std(errSS1, 0, 1)  / sqrt(nParts);
mBias1 = mean(biasSS1, 1);  seBias1 = std(biasSS1, 0, 1) / sqrt(nParts);

%% 4. Sweep A: Deterministic Spatial Readout Width (sigma_ang = 0 deg)
sigmaSweep = [32, 50, 70, 90, 110];
nSig = numel(sigmaSweep);
swA_err  = zeros(nSig, nRanges); swA_errSE  = zeros(nSig, nRanges);
swA_bias = zeros(nSig, nRanges); swA_biasSE = zeros(nSig, nRanges);
swA_prop = zeros(nSig, nRanges); swA_propSE = zeros(nSig, nRanges);

fprintf('\n=== SWEEP A: Deterministic Spatial Readout Width (itemSD = %d px, sigma_ang = 0 deg) ===\n', itemSD);
fprintf('%-10s | %-18s | %-18s | %-18s\n', 'sigma_pos', 'Err [30,60,120]', 'Bias [30,60,120]', 'Prop [30,60,120]');
fprintf('%s\n', repmat('-', 1, 76));

swA_rows = table();
for k = 1:nSig
    sigK = sigmaSweep(k);
    [eK, bK] = runRingSpatialSim(betaV1_obs, baseScale, loadScale, cx, cy, ringRadius, cornerAng, ...
        locX, locY, itemSD, sigK, 0, baseOffs, nTrials, phiDeg, phiRad, 2000 + k);
    pK = bK ./ meanDist;
    swA_err(k,:)  = mean(eK, 1);  swA_errSE(k,:)  = std(eK, 0, 1) / sqrt(nParts);
    swA_bias(k,:) = mean(bK, 1);  swA_biasSE(k,:) = std(bK, 0, 1) / sqrt(nParts);
    swA_prop(k,:) = mean(pK, 1);  swA_propSE(k,:) = std(pK, 0, 1) / sqrt(nParts);

    fprintf('%4d px    | %4.1f, %4.1f, %4.1f   | %+4.1f, %+4.1f, %+4.1f | %.2f, %.2f, %.2f\n', ...
        sigK, swA_err(k,1), swA_err(k,2), swA_err(k,3), ...
        swA_bias(k,1), swA_bias(k,2), swA_bias(k,3), ...
        swA_prop(k,1), swA_prop(k,2), swA_prop(k,3));

    for rIdx = 1:nRanges
        swA_rows = [swA_rows; table(sigK, 0, rangesDeg(rIdx), ...
            swA_err(k,rIdx), swA_errSE(k,rIdx), swA_bias(k,rIdx), swA_biasSE(k,rIdx), ...
            swA_prop(k,rIdx), swA_propSE(k,rIdx), ...
            'VariableNames', {'SigmaRead_px', 'SigmaAng_deg', 'Range_deg', ...
            'AngDev_SS4', 'SE_AngDev', 'Bias_SS4', 'SE_Bias', 'PropBias_SS4', 'SE_PropBias'})]; %#ok<AGROW>
    end
end
writetable(swA_rows, 'simulated_sigma_sweep.csv');

%% 5. Sweep B: Fixed Narrow Window (sigma_read = 32 px) + Ring-Anchored Localisation Jitter
angJitterSweep = [0, 35, 50, 65, 80]; % angular jitter SD along the 4-item ring (deg)
nJit = numel(angJitterSweep);
swB_err  = zeros(nJit, nRanges); swB_errSE  = zeros(nJit, nRanges);
swB_bias = zeros(nJit, nRanges); swB_biasSE = zeros(nJit, nRanges);
swB_prop = zeros(nJit, nRanges); swB_propSE = zeros(nJit, nRanges);

fprintf('\n=== SWEEP B: Fixed Local Window (sigma_read = %d px) + Ring Localisation Jitter ===\n', sigmaReadFixed);
fprintf('%-10s | %-18s | %-18s | %-18s\n', 'sigma_ang', 'Err [30,60,120]', 'Bias [30,60,120]', 'Prop [30,60,120]');
fprintf('%s\n', repmat('-', 1, 76));

swB_rows = table();
primaryJitIdx = 4; % sigma_ang = 65 deg (featured in primary 4-panel figure)
errSS4_main = []; biasSS4_main = [];

for k = 1:nJit
    angK = angJitterSweep(k);
    [eK, bK] = runRingSpatialSim(betaV1_obs, baseScale, loadScale, cx, cy, ringRadius, cornerAng, ...
        locX, locY, itemSD, sigmaReadFixed, angK, baseOffs, nTrials, phiDeg, phiRad, 3000 + k);
    pK = bK ./ meanDist;
    swB_err(k,:)  = mean(eK, 1);  swB_errSE(k,:)  = std(eK, 0, 1) / sqrt(nParts);
    swB_bias(k,:) = mean(bK, 1);  swB_biasSE(k,:) = std(bK, 0, 1) / sqrt(nParts);
    swB_prop(k,:) = mean(pK, 1);  swB_propSE(k,:) = std(pK, 0, 1) / sqrt(nParts);

    if k == primaryJitIdx
        errSS4_main  = eK;
        biasSS4_main = bK;
    end

    fprintf('%4d deg   | %4.1f, %4.1f, %4.1f   | %+4.1f, %+4.1f, %+4.1f | %.2f, %.2f, %.2f\n', ...
        angK, swB_err(k,1), swB_err(k,2), swB_err(k,3), ...
        swB_bias(k,1), swB_bias(k,2), swB_bias(k,3), ...
        swB_prop(k,1), swB_prop(k,2), swB_prop(k,3));

    for rIdx = 1:nRanges
        swB_rows = [swB_rows; table(sigmaReadFixed, angK, rangesDeg(rIdx), ...
            swB_err(k,rIdx), swB_errSE(k,rIdx), swB_bias(k,rIdx), swB_biasSE(k,rIdx), ...
            swB_prop(k,rIdx), swB_propSE(k,rIdx), ...
            'VariableNames', {'SigmaRead_px', 'SigmaAng_deg', 'Range_deg', ...
            'AngDev_SS4', 'SE_AngDev', 'Bias_SS4', 'SE_Bias', 'PropBias_SS4', 'SE_PropBias'})]; %#ok<AGROW>
    end
end
writetable(swB_rows, 'simulated_jitter_sweep.csv');

%% 6. Print Primary Model Stats & Save Summary CSV
mErr4  = mean(errSS4_main, 1);   seErr4  = std(errSS4_main, 0, 1)  / sqrt(nParts);
mBias4 = mean(biasSS4_main, 1);  seBias4 = std(biasSS4_main, 0, 1) / sqrt(nParts);
propBias4 = biasSS4_main ./ meanDist;
mProp4 = mean(propBias4, 1);     seProp4 = std(propBias4, 0, 1)    / sqrt(nParts);

summaryTbl = table(rangesDeg', mErr1', seErr1', mErr4', seErr4', ...
    mBias1', seBias1', mBias4', seBias4', mProp4', seProp4', ...
    'VariableNames', {'Range_deg', ...
    'AngDev_SS1', 'SE_AngDev_SS1', 'AngDev_SS4', 'SE_AngDev_SS4', ...
    'Bias_SS1', 'SE_Bias_SS1', 'Bias_SS4', 'SE_Bias_SS4', ...
    'PropBias_SS4', 'SE_PropBias_SS4'});
writetable(summaryTbl, 'simulated_ensemble_bias_summary.csv');

[~, pErr4,  ~, stErr4]  = ttest(errSS4_main(:, 3), errSS4_main(:, 1));
[~, pBias4, ~, stBias4] = ttest(biasSS4_main(:, 3), biasSS4_main(:, 1));
[~, pProp4, ~, stProp4] = ttest(propBias4(:, 1), propBias4(:, 3));
fprintf('\nPrimary Fixed-Local + Ring-Jitter Model (sigma_read = %d px, sigma_ang = %d deg, N = %d):\n', ...
    sigmaReadFixed, angJitterSweep(primaryJitIdx), nParts);
fprintf('  Error increases (120° > 30°):             t(%d) = %.2f, p = %.2e (%d/%d observers)\n', ...
    stErr4.df, stErr4.tstat, pErr4, sum(errSS4_main(:,3) > errSS4_main(:,1)), nParts);
fprintf('  Absolute bias increases (120° > 30°):     t(%d) = %.2f, p = %.2e (%d/%d observers)\n', ...
    stBias4.df, stBias4.tstat, pBias4, sum(biasSS4_main(:,3) > biasSS4_main(:,1)), nParts);
fprintf('  Proportional bias decreases (30° > 120°): t(%d) = %.2f, p = %.2e (%d/%d observers)\n', ...
    stProp4.df, stProp4.tstat, pProp4, sum(propBias4(:,1) > propBias4(:,3)), nParts);

%% 7. Figure 1: Primary 4-Panel Display + Fixed-Local Ring-Jitter Replication
colSS1 = [0.00, 0.45, 0.74];
colSS4 = [0.85, 0.33, 0.10];

exOris = mod(35 + baseOffs{2}, 180);
exImg  = zeros(imSize, imSize);
for j = 1:nItems
    th = deg2rad(exOris(j));
    xr = (X - locX(j)) * cos(th) + (Y - locY(j)) * sin(th);
    exImg = exImg + E_maps(:,:,j) .* cos(2*pi*gabFreq * xr);
end

f1 = figure('Color', 'w', 'Visible', 'off', 'Units', 'pixels', 'Position', [60, 140, 1480, 380]);
tiledlayout(1, 4, 'TileSpacing', 'loose', 'Padding', 'compact');

% Panel A: Lower-SF 4-Gabor Display + Fixed 32 px Local Window + Ring-Jittered Readout Centres
nexttile; hold on;
imagesc(exImg, [-1, 1]); colormap(gca, gray); axis image off;
thC = linspace(0, 2*pi, 160);
cueEx = 1;
% Show subtle stimulus ring & 35 example ring-jittered readout centres
plot(cx + ringRadius*cos(thC), cy + ringRadius*sin(thC), ':', 'Color', [1 1 1 0.45], 'LineWidth', 1.0);
rng(42);
jitAngs = cornerAng(cueEx) + deg2rad(randn(35, 1) * angJitterSweep(primaryJitIdx));
scatter(cx + ringRadius*cos(jitAngs), cy + ringRadius*sin(jitAngs), 18, colSS4, 'filled', 'MarkerFaceAlpha', 0.55);
plot(locX(cueEx) + 22*cos(thC), locY(cueEx) + 22*sin(thC), '-', 'Color', [0.95 0.78 0.18], 'LineWidth', 2.0);
plot(locX(cueEx) + sigmaReadFixed*cos(thC), locY(cueEx) + sigmaReadFixed*sin(thC), '--', 'Color', colSS1, 'LineWidth', 2.0);
xlim([1, imSize]); ylim([1, imSize]);
title('A. Display & Local Readout', 'FontSize', 11.0, 'FontWeight', 'bold');
text(6, 236, {'Yellow: cued | Blue: 32px window', 'Red dots: SS4 ring jitter'}, ...
    'Color', 'w', 'FontSize', 7.4, 'BackgroundColor', [0 0 0 0.65]);

% Panel B: Error
nexttile; hold on;
errorbar(rangesDeg, mErr1, seErr1, '-', 'Color', colSS1, 'LineWidth', 2.2, 'CapSize', 5);
errorbar(rangesDeg, mErr4, seErr4, '-', 'Color', colSS4, 'LineWidth', 2.2, 'CapSize', 5);
xlim([10, 140]); ylim([0, 60]);
set(gca, 'XTick', rangesDeg, 'YTick', [0, 20, 40, 60], 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
xlabel('Range (deg)', 'FontSize', 11); ylabel('Angular deviation (deg)', 'FontSize', 11);
title('B. Error (32px + Jitter)', 'FontSize', 10.8, 'FontWeight', 'bold');
legend({'Set size 1', 'Set size 4'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9.5);

% Panel C: Absolute Bias
nexttile; hold on;
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off');
errorbar(rangesDeg, mBias1, seBias1, '-', 'Color', colSS1, 'LineWidth', 2.2, 'CapSize', 5);
errorbar(rangesDeg, mBias4, seBias4, '-', 'Color', colSS4, 'LineWidth', 2.2, 'CapSize', 5);
xlim([10, 140]); ylim([-30, 35]);
set(gca, 'XTick', rangesDeg, 'YTick', [-30, -20, -10, 0, 10, 20, 30], 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
xlabel('Range (deg)', 'FontSize', 11); ylabel('Shift toward mean (circmean, deg)', 'FontSize', 11);
title('C. Absolute Bias', 'FontSize', 10.8, 'FontWeight', 'bold');

% Panel D: Proportional Bias
nexttile; hold on;
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off');
errorbar(rangesDeg, mProp4, seProp4, '-o', 'Color', colSS4, 'MarkerFaceColor', colSS4, ...
    'LineWidth', 2.2, 'MarkerSize', 6, 'CapSize', 5);
xlim([10, 140]); ylim([0, 1.0]);
set(gca, 'XTick', rangesDeg, 'YTick', [0, 0.25, 0.5, 0.75, 1.0], 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
xlabel('Range (deg)', 'FontSize', 11); ylabel('Proportional shift toward mean', 'FontSize', 11);
title('D. Proportional Bias (SS4)', 'FontSize', 10.8, 'FontWeight', 'bold');

exportgraphics(f1, 'fig_simulated_ensemble_bias.png', 'Resolution', 300);
exportgraphics(f1, 'fig_simulated_ensemble_bias.pdf', 'ContentType', 'vector');
set(f1, 'Visible', 'on');

%% 8. Figure 2: 2x3 Sensitivity Sweep Comparison (Deterministic sigma_pos vs. Fixed 32px + Ring Jitter)
cols = parula(nSig + 1);
f2 = figure('Color', 'w', 'Position', [80, 80, 1180, 680]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top row: Sweep A (Deterministic sigma_pos)
nexttile; hold on;
for k = 1:nSig
    errorbar(rangesDeg, swA_err(k,:), swA_errSE(k,:), '-o', 'Color', cols(k,:), 'LineWidth', 1.8, 'MarkerSize', 4);
end
xlim([10, 140]); ylim([0, 60]); set(gca, 'XTick', rangesDeg, 'TickDir', 'out', 'Box', 'off');
ylabel('Angular deviation (deg)'); title('A. Deterministic \sigma_{pos}: Error');
legend(arrayfun(@(x) sprintf('\\sigma_{pos}=%dpx', x), sigmaSweep, 'UniformOutput', false), 'Location', 'northwest', 'Box', 'off', 'FontSize', 8.5);

nexttile; hold on; yline(0, '--', 'Color', [0.6 0.6 0.6]);
for k = 1:nSig
    errorbar(rangesDeg, swA_bias(k,:), swA_biasSE(k,:), '-o', 'Color', cols(k,:), 'LineWidth', 1.8, 'MarkerSize', 4);
end
xlim([10, 140]); ylim([-5, 38]); set(gca, 'XTick', rangesDeg, 'TickDir', 'out', 'Box', 'off');
ylabel('Absolute bias (deg)'); title('B. Deterministic \sigma_{pos}: Bias');

nexttile; hold on; yline(0, '--', 'Color', [0.6 0.6 0.6]);
for k = 1:nSig
    errorbar(rangesDeg, swA_prop(k,:), swA_propSE(k,:), '-o', 'Color', cols(k,:), 'LineWidth', 1.8, 'MarkerSize', 4);
end
xlim([10, 140]); ylim([0, 1.0]); set(gca, 'XTick', rangesDeg, 'TickDir', 'out', 'Box', 'off');
ylabel('Proportional bias'); title('C. Deterministic \sigma_{pos}: Prop. Bias');

% Bottom row: Sweep B (Fixed 32 px Local Window + Ring Localisation Jitter)
nexttile; hold on;
for k = 1:nJit
    errorbar(rangesDeg, swB_err(k,:), swB_errSE(k,:), '-s', 'Color', cols(k,:), 'LineWidth', 1.8, 'MarkerSize', 4);
end
xlim([10, 140]); ylim([0, 60]); set(gca, 'XTick', rangesDeg, 'TickDir', 'out', 'Box', 'off');
xlabel('Range (deg)'); ylabel('Angular deviation (deg)'); title(sprintf('D. Fixed Local (%dpx) + Ring Jitter: Error', sigmaReadFixed));
legend(arrayfun(@(x) sprintf('\\sigma_{ang}=%d°', x), angJitterSweep, 'UniformOutput', false), 'Location', 'northwest', 'Box', 'off', 'FontSize', 8.5);

nexttile; hold on; yline(0, '--', 'Color', [0.6 0.6 0.6]);
for k = 1:nJit
    errorbar(rangesDeg, swB_bias(k,:), swB_biasSE(k,:), '-s', 'Color', cols(k,:), 'LineWidth', 1.8, 'MarkerSize', 4);
end
xlim([10, 140]); ylim([-5, 38]); set(gca, 'XTick', rangesDeg, 'TickDir', 'out', 'Box', 'off');
xlabel('Range (deg)'); ylabel('Absolute bias (deg)'); title(sprintf('E. Fixed Local (%dpx) + Ring Jitter: Bias', sigmaReadFixed));

nexttile; hold on; yline(0, '--', 'Color', [0.6 0.6 0.6]);
for k = 1:nJit
    errorbar(rangesDeg, swB_prop(k,:), swB_propSE(k,:), '-s', 'Color', cols(k,:), 'LineWidth', 1.8, 'MarkerSize', 4);
end
xlim([10, 140]); ylim([0, 1.0]); set(gca, 'XTick', rangesDeg, 'TickDir', 'out', 'Box', 'off');
xlabel('Range (deg)'); ylabel('Proportional bias'); title(sprintf('F. Fixed Local (%dpx) + Ring Jitter: Prop. Bias', sigmaReadFixed));

exportgraphics(f2, 'fig_simulated_sensitivity_sweeps.png', 'Resolution', 300);
exportgraphics(f2, 'fig_simulated_sensitivity_sweeps.pdf', 'ContentType', 'vector');
fprintf('\nSaved fig_simulated_ensemble_bias.png/.pdf, fig_simulated_sensitivity_sweeps.png/.pdf, and sweep CSVs.\n');

%% Local Helpers
function [angDevMat, circBiasMat] = runRingSpatialSim(betaV1_obs, baseScale, loadScale, ...
    cx, cy, ringRadius, cornerAng, locX, locY, itemSD, sigmaRead, sigmaAngDeg, baseOffs, nTrials, phiDeg, phiRad, seedVal)
    rng(seedVal);
    nParts  = numel(betaV1_obs);
    nRanges = numel(baseOffs);
    nItems  = numel(locX);
    angDevMat   = zeros(nParts, nRanges);
    circBiasMat = zeros(nParts, nRanges);

    effVar = sigmaRead^2 + itemSD^2;

    for p = 1:nParts
        gain = betaV1_obs(p) * baseScale * loadScale;
        for rIdx = 1:nRanges
            offs = baseOffs{rIdx};
            errUnit = zeros(nTrials, 1);

            for tr = 1:nTrials
                muTrue  = rand * 360;
                jitter  = (rand(1, nItems) * 6) - 3;
                itemDeg = mod(muTrue + offs + jitter, 360);
                itemDeg = itemDeg(randperm(nItems));
                itemRad = deg2rad(itemDeg);

                cue     = randi(nItems);
                targDeg = itemDeg(cue);
                dTarg   = circDiff360(targDeg, muTrue);

                % Ring-anchored localisation jitter along the occupied stimulus circle
                if sigmaAngDeg > 0
                    rAng = cornerAng(cue) + deg2rad(randn * sigmaAngDeg);
                    rx   = cx + ringRadius * cos(rAng);
                    ry   = cy + ringRadius * sin(rAng);
                else
                    rx = locX(cue);
                    ry = locY(cue);
                end

                % Strictly local V1 readout window at (rx, ry)
                d2 = (rx - locX).^2 + (ry - locY).^2;
                w  = exp(-d2 / (2 * effVar));
                S_phi = (w * cos(itemRad' - phiRad'))' / (norm(w) + eps);

                resp = sampleSoftmax(phiDeg, S_phi, gain);
                e    = circDiff360(resp, targDeg);
                errUnit(tr) = -sign(dTarg) * e;
            end

            [circBiasMat(p, rIdx), angDevMat(p, rIdx)] = circMetrics360(errUnit);
        end
    end
end

function respDeg = sampleSoftmax(phiDeg, S_phi, gain)
    logits = gain * (S_phi - max(S_phi));
    probs  = exp(logits);
    probs  = probs / sum(probs);
    respDeg = phiDeg(find(cumsum(probs) >= rand, 1, 'first'));
end

function d = circDiff360(a, b)
    d = mod(a - b + 180, 360) - 180;
end

function [cMeanDeg, angDevDeg] = circMetrics360(errDeg)
    r = mean(exp(1i * deg2rad(errDeg)));
    cMeanDeg  = rad2deg(angle(r));
    angDevDeg = rad2deg(sqrt(max(2 * (1 - abs(r)), 0)));
end
