%% simulate_ensemble_bias_V1.m
% Directly simulates Utochkin & Brady (2020, JEP:HPP, Experiment 1, Figure 4)
% using our winning purely spatial V1 Oriented Matched-Filter observer
% parameterised by ALL 30 of our real study participants' individually fitted
% V1_Spatial sensitivities from Experiments 1-3.
%
% Exact Utochkin & Brady (2020) Experiment 1 paradigm:
%   - 4 directional items placed at the 4 corners (45, 135, 225, 315 deg,
%     jittered by ±10 deg around the circle).
%   - 3 orientation ranges around a random mean in [0, 360):
%       Range  30 deg: [-15, -5, +5, +15] ± 3 deg jitter  (mean dist = 10 deg)
%       Range  60 deg: [-30, -10, +10, +30] ± 3 deg jitter (mean dist = 20 deg)
%       Range 120 deg: [-60, -20, +20, +60] ± 3 deg jitter (mean dist = 40 deg)
%   - 2 individual-item memory tasks:
%       Set size 1 ("Remember One"): precued single corner -> tight spatial
%         readout kernel (minimal neighbour leakage) & full single-item gain.
%       Set size 4 ("Remember Four"): all 4 corners maintained -> broader
%         spatial readout uncertainty in memory + set-size capacity scaling.
%   - Exact nonparametric metrics from Utochkin & Brady (2020) Figure 4:
%       Panel A (Error): Angular deviation sqrt(2*(1 - R)) in deg (Berens, 2009)
%       Panel B (Bias):  Circular mean of unitized error (+ = toward mean)
%       Panel C:         Proportional shift toward mean (bias / mean distance)
%
% Outputs saved in Working/:
%   - fig_simulated_ensemble_bias.png / .pdf
%   - simulated_ensemble_bias_summary.csv

clear; close all; clc;
rng(2026);

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

%% 2. Utochkin & Brady (2020) Experiment 1 Geometry & Parameters
nItems     = 4;
rangesDeg  = [30, 60, 120];
nRanges    = numel(rangesDeg);
baseOffs   = {[-15, -5, 5, 15], [-30, -10, 10, 30], [-60, -20, 20, 60]};
meanDist   = [10, 20, 40];  % theoretical mean distance to display mean (deg)
nTrials    = 240;           % trials per cell per observer (80 in U&B; 240 for smooth SEs)

% Spatial display geometry (4 corners on circle of radius R = 72 px; adjacent spacing ~102 px)
imSize     = 256;
ringRadius = 72;
itemSD     = 10;            % item spatial envelope SD (px)
[X, Y]     = meshgrid(1:imSize, 1:imSize);
cx = imSize/2; cy = imSize/2;
cornerAng  = deg2rad([45, 135, 225, 315]);
locX       = cx + ringRadius * cos(cornerAng);
locY       = cy + ringRadius * sin(cornerAng);

E_maps = zeros(imSize, imSize, nItems);
for j = 1:nItems
    E_maps(:,:,j) = exp(-((X - locX(j)).^2 + (Y - locY(j)).^2) / (2 * itemSD^2));
end

% Spatial readout kernel widths in working memory (ZERO ensemble computation):
%   Set size 1 (precued corner): tight spatial focus (sigma_pos = 32 px -> ~2% neighbour overlap)
%   Set size 4 (all 4 corners):  broader spatial readout uncertainty (sigma_pos = 105 px)
sigmaPos_SS1 = 32;
sigmaPos_SS4 = 105;

wMat_SS1 = computeOverlapMatrix(X, Y, locX, locY, E_maps, sigmaPos_SS1);
wMat_SS4 = computeOverlapMatrix(X, Y, locX, locY, E_maps, sigmaPos_SS4);

fprintf('Spatial overlap weights (relative to cued corner = 1.0):\n');
fprintf('  Set size 1 (sigma_pos = %3d px): adjacent = %.3f, opposite = %.3f\n', ...
    sigmaPos_SS1, wMat_SS1(1,2), wMat_SS1(1,3));
fprintf('  Set size 4 (sigma_pos = %3d px): adjacent = %.3f, opposite = %.3f\n', ...
    sigmaPos_SS4, wMat_SS4(1,2), wMat_SS4(1,3));

% Candidate response orientations in full 360-deg directional space (like U&B triangles)
phiDeg = (0:1:359)';
phiRad = deg2rad(phiDeg);

%% 3. Simulate All 30 Real Participants on Set Size 1 vs. Set Size 4 x Ranges (30, 60, 120)
% Output matrices: [nParts x nRanges] for Set Size 1 and Set Size 4
angDev_SS1   = zeros(nParts, nRanges);
angDev_SS4   = zeros(nParts, nRanges);
circBias_SS1 = zeros(nParts, nRanges);
circBias_SS4 = zeros(nParts, nRanges);

% Calibrate base decision sharpness to each observer's beta_V1
% (Set size 1 -> ~18 deg baseline angular deviation; Set size 4 -> 1/sqrt(4) = 0.5x resource scaling)
baseScale = 19.5;

fprintf('Simulating N = %d observers across 2 set sizes x 3 ranges (%d trials/cell)...\n', nParts, nTrials);

for p = 1:nParts
    gainSS1 = betaV1_obs(p) * baseScale;
    gainSS4 = gainSS1 * 0.45; % standard working-memory set-size load scaling (1 -> 4 items)

    for rIdx = 1:nRanges
        offs = baseOffs{rIdx};

        errUnit_SS1 = zeros(nTrials, 1);
        errUnit_SS4 = zeros(nTrials, 1);

        for tr = 1:nTrials
            % Generate exact Utochkin & Brady (2020) 4-item display
            muTrue  = rand * 360;
            jitter  = (rand(1, nItems) * 6) - 3; % uniform [-3, +3] deg jitter
            itemDeg = mod(muTrue + offs + jitter, 360);
            % Randomly assign orientations to the 4 corner locations
            itemDeg = itemDeg(randperm(nItems));
            itemRad = deg2rad(itemDeg);

            % Probe one random corner (1..4)
            cue     = randi(nItems);
            targDeg = itemDeg(cue);
            % Signed offset of probed item relative to display mean in [-180, 180)
            dTarg   = circDiff360(targDeg, muTrue);

            % --- Set size 1: tight spatial readout at precued corner ---
            w1 = wMat_SS1(cue, :);
            S1 = (w1 * cos(itemRad' - phiRad'))' / (norm(w1) + eps);
            resp1 = sampleSoftmax(phiDeg, S1, gainSS1);
            e1 = circDiff360(resp1, targDeg);
            % Unitize sign: positive = toward mean, negative = away from mean
            errUnit_SS1(tr) = -sign(dTarg) * e1;

            % --- Set size 4: broader spatial readout uncertainty in memory ---
            w4 = wMat_SS4(cue, :);
            S4 = (w4 * cos(itemRad' - phiRad'))' / (norm(w4) + eps);
            resp4 = sampleSoftmax(phiDeg, S4, gainSS4);
            e4 = circDiff360(resp4, targDeg);
            errUnit_SS4(tr) = -sign(dTarg) * e4;
        end

        % Exact Utochkin & Brady (2020) nonparametric metrics (Berens, 2009)
        [circBias_SS1(p, rIdx), angDev_SS1(p, rIdx)] = circMetrics360(errUnit_SS1);
        [circBias_SS4(p, rIdx), angDev_SS4(p, rIdx)] = circMetrics360(errUnit_SS4);
    end
end

%% 4. Summarise Across N = 30 Simulated Participants & Save CSV
mErr1  = mean(angDev_SS1, 1);   seErr1  = std(angDev_SS1, 0, 1)   / sqrt(nParts);
mErr4  = mean(angDev_SS4, 1);   seErr4  = std(angDev_SS4, 0, 1)   / sqrt(nParts);
mBias1 = mean(circBias_SS1, 1); seBias1 = std(circBias_SS1, 0, 1) / sqrt(nParts);
mBias4 = mean(circBias_SS4, 1); seBias4 = std(circBias_SS4, 0, 1) / sqrt(nParts);

propBias4 = circBias_SS4 ./ meanDist; % proportional shift toward mean (+10, +20, +40 deg)
mProp4  = mean(propBias4, 1);   seProp4 = std(propBias4, 0, 1)    / sqrt(nParts);

summaryTbl = table(rangesDeg', mErr1', seErr1', mErr4', seErr4', ...
    mBias1', seBias1', mBias4', seBias4', mProp4', seProp4', ...
    'VariableNames', {'Range_deg', ...
    'AngDev_SS1', 'SE_AngDev_SS1', 'AngDev_SS4', 'SE_AngDev_SS4', ...
    'Bias_SS1', 'SE_Bias_SS1', 'Bias_SS4', 'SE_Bias_SS4', ...
    'PropBias_SS4', 'SE_PropBias_SS4'});
writetable(summaryTbl, 'simulated_ensemble_bias_summary.csv');

fprintf('\n=== SIMULATED UTOCHKIN & BRADY (2020) FIG 4 RESULTS (V1_Spatial, N = %d) ===\n', nParts);
fprintf('%-8s | %-16s | %-16s | %-16s | %-16s | %-14s\n', ...
    'Range', 'Error SS1 (deg)', 'Error SS4 (deg)', 'Bias SS1 (deg)', 'Bias SS4 (deg)', 'Prop Bias SS4');
fprintf('%s\n', repmat('-', 1, 100));
for rIdx = 1:nRanges
    fprintf('%3d deg  | %5.1f ± %4.1f     | %5.1f ± %4.1f     | %+5.1f ± %4.1f     | %+5.1f ± %4.1f     | %.2f ± %.2f\n', ...
        rangesDeg(rIdx), mErr1(rIdx), seErr1(rIdx), mErr4(rIdx), seErr4(rIdx), ...
        mBias1(rIdx), seBias1(rIdx), mBias4(rIdx), seBias4(rIdx), mProp4(rIdx), seProp4(rIdx));
end

[~, pErr4,  ~, stErr4]  = ttest(angDev_SS4(:, 3), angDev_SS4(:, 1));
[~, pBias4, ~, stBias4] = ttest(circBias_SS4(:, 3), circBias_SS4(:, 1));
[~, pProp4, ~, stProp4] = ttest(propBias4(:, 1), propBias4(:, 3));
fprintf('\nSet size 4 range checks (30° vs 120°, N = %d):\n', nParts);
fprintf('  Error increases (120° > 30°):             t(%d) = %.2f, p = %.2e (%d/%d observers)\n', ...
    stErr4.df, stErr4.tstat, pErr4, sum(angDev_SS4(:,3) > angDev_SS4(:,1)), nParts);
fprintf('  Absolute bias increases (120° > 30°):     t(%d) = %.2f, p = %.2e (%d/%d observers)\n', ...
    stBias4.df, stBias4.tstat, pBias4, sum(circBias_SS4(:,3) > circBias_SS4(:,1)), nParts);
fprintf('  Proportional bias decreases (30° > 120°): t(%d) = %.2f, p = %.2e (%d/%d observers)\n', ...
    stProp4.df, stProp4.tstat, pProp4, sum(propBias4(:,1) > propBias4(:,3)), nParts);

%% 5. Generate Publication-Ready 4-Panel Figure (Example Display + U&B Fig 4 Replication)
colSS1 = [0.00, 0.45, 0.74]; % MATLAB blue (Set size 1)
colSS4 = [0.85, 0.33, 0.10]; % MATLAB orange-red (Set size 4)

% Render a clean example 4-Gabor corner display (60 deg range) for Panel A
gabFreq = 0.12;
exOris  = mod(35 + baseOffs{2}, 180); % 60 deg range around 35 deg
exImg   = zeros(imSize, imSize);
for j = 1:nItems
    th = deg2rad(exOris(j));
    xr = (X - locX(j)) * cos(th) + (Y - locY(j)) * sin(th);
    exImg = exImg + E_maps(:,:,j) .* cos(2*pi*gabFreq * xr);
end

f = figure('Color', 'w', 'Position', [60, 140, 1420, 360]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel A: Example 4-Item Display with Set Size 1 & Set Size 4 Spatial Readout Kernels
nexttile; hold on;
imagesc(exImg, [-1, 1]); colormap(gca, gray); axis image off;
thC = linspace(0, 2*pi, 160);
cueEx = 1; % cued corner
plot(locX(cueEx) + 18*cos(thC), locY(cueEx) + 18*sin(thC), '-', 'Color', [0.95 0.78 0.18], 'LineWidth', 2.0);
plot(locX(cueEx) + sigmaPos_SS1*cos(thC), locY(cueEx) + sigmaPos_SS1*sin(thC), '--', 'Color', colSS1, 'LineWidth', 2.0);
plot(locX(cueEx) + sigmaPos_SS4*cos(thC), locY(cueEx) + sigmaPos_SS4*sin(thC), '--', 'Color', colSS4, 'LineWidth', 2.0);
xlim([1, imSize]); ylim([1, imSize]);
title('A. Display & Readout Windows', 'FontSize', 11.5, 'FontWeight', 'bold');
text(6, 242, 'Yellow: cued | Blue: SS1 (1\sigma) | Red: SS4 (1\sigma)', ...
    'Color', 'w', 'FontSize', 7.8, 'BackgroundColor', [0 0 0 0.60]);

% Panel B: Error (Angular Deviation) — directly matches Utochkin & Brady Fig 4A
nexttile; hold on;
errorbar(rangesDeg, mErr1, seErr1, '-', 'Color', colSS1, 'LineWidth', 2.2, 'CapSize', 5);
errorbar(rangesDeg, mErr4, seErr4, '-', 'Color', colSS4, 'LineWidth', 2.2, 'CapSize', 5);
xlim([10, 140]); ylim([0, 60]);
set(gca, 'XTick', rangesDeg, 'YTick', [0, 20, 40, 60], 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
xlabel('Range (deg)', 'FontSize', 11);
ylabel('Angular deviation (deg)', 'FontSize', 11);
title('B. Error (V1 Spatial Model)', 'FontSize', 11.5, 'FontWeight', 'bold');
legend({'Set size 1', 'Set size 4'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9.5);

% Panel C: Bias (Circular Mean Shift Toward Mean) — directly matches Utochkin & Brady Fig 4B
nexttile; hold on;
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off');
errorbar(rangesDeg, mBias1, seBias1, '-', 'Color', colSS1, 'LineWidth', 2.2, 'CapSize', 5);
errorbar(rangesDeg, mBias4, seBias4, '-', 'Color', colSS4, 'LineWidth', 2.2, 'CapSize', 5);
xlim([10, 140]); ylim([-30, 35]);
set(gca, 'XTick', rangesDeg, 'YTick', [-30, -20, -10, 0, 10, 20, 30], 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
xlabel('Range (deg)', 'FontSize', 11);
ylabel('Shift toward mean (circmean, deg)', 'FontSize', 11);
title('C. Absolute Bias (V1 Spatial Model)', 'FontSize', 11.5, 'FontWeight', 'bold');

% Panel D: Proportional Bias Toward Mean (decreases with range, matching U&B)
nexttile; hold on;
yline(0, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off');
errorbar(rangesDeg, mProp4, seProp4, '-o', 'Color', colSS4, 'MarkerFaceColor', colSS4, ...
    'LineWidth', 2.2, 'MarkerSize', 6, 'CapSize', 5);
xlim([10, 140]); ylim([0, 1.0]);
set(gca, 'XTick', rangesDeg, 'YTick', [0, 0.25, 0.5, 0.75, 1.0], 'TickDir', 'out', 'Box', 'off', 'FontSize', 11);
xlabel('Range (deg)', 'FontSize', 11);
ylabel('Proportional shift toward mean', 'FontSize', 11);
title('D. Proportional Bias (Set Size 4)', 'FontSize', 11.5, 'FontWeight', 'bold');

exportgraphics(f, 'fig_simulated_ensemble_bias.png', 'Resolution', 300);
exportgraphics(f, 'fig_simulated_ensemble_bias.pdf', 'ContentType', 'vector');
fprintf('\nSaved 4-panel fig_simulated_ensemble_bias.png/.pdf and simulated_ensemble_bias_summary.csv.\n');

%% Local Helpers
function wMat = computeOverlapMatrix(X, Y, locX, locY, E_maps, sigP)
    n = numel(locX);
    wMat = zeros(n, n);
    for cue = 1:n
        W_cue = exp(-((X - locX(cue)).^2 + (Y - locY(cue)).^2) / (2 * sigP^2));
        for j = 1:n
            wMat(cue, j) = sum(W_cue .* E_maps(:,:,j), 'all');
        end
        wMat(cue, :) = wMat(cue, :) / wMat(cue, cue);
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
    % Exact nonparametric circular mean & angular deviation (Berens, 2009; Zar, 1999 eq. 26.20)
    r = mean(exp(1i * deg2rad(errDeg)));
    cMeanDeg  = rad2deg(angle(r));
    angDevDeg = rad2deg(sqrt(max(2 * (1 - abs(r)), 0)));
end
