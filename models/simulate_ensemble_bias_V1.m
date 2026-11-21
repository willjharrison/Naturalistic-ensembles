%% simulate_ensemble_bias_V1.m
% Simulates an 8-Gabor working-memory orientation-recall task for ALL 30
% real study participants using EACH observer's own individually fitted
% parameters from Experiments 1-3 (cached unifiedPredictors_v2_V1.mat).
%
% Research question:
%   "If our exact 30 participants (characterised by their own winning
%   V1_Spatial sensitivities from Experiments 1-3) completed a classic
%   working-memory orientation-recall task, would the purely spatial V1
%   model predict a group-level ensemble bias—and its reduction with wider
%   ensemble SD—without ever computing an ensemble mean or SD?"
%
% Outputs saved in Working/:
%   - fig_simulated_ensemble_bias.png / .pdf
%   - simulated_ensemble_bias_summary.csv

clear; close all; clc;
rng(2026);

%% 1. Extract Each of the 30 Observers' Individual Parameters from Exps 1-3
fprintf('Loading cached unifiedPredictors_v2_V1.mat to fit individual observer parameters...\n');
load('unifiedPredictors_v2_V1.mat', 'leftV1Mat', 'rightV1Mat', 'absMuMat', 'absSDMat', 'validTrial');

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
muRaw = absMuMat(validTrial, cfg);
sdRaw = absSDMat(validTrial, cfg);

uParts = unique(partIDVec);
nParts = numel(uParts);
betaV1_obs  = zeros(nParts, 1);
betaEns_obs = zeros(nParts, 1);

for p = 1:nParts
    idxP = (partIDVec == uParts(p));
    Tp = table(respVec(idxP), zscore(v1Raw(idxP)), zscore(muRaw(idxP)), zscore(sdRaw(idxP)), ...
        'VariableNames', {'resp', 'v1Diff', 'absMuDiff', 'absSDDiff'});
    mV1 = fitglm(Tp, 'resp ~ v1Diff', 'Distribution', 'Binomial');
    mEn = fitglm(Tp, 'resp ~ absMuDiff + absSDDiff', 'Distribution', 'Binomial');
    betaV1_obs(p)  = max(mV1.Coefficients.Estimate(2), 0.05); % positive sensitivity floor
    betaEns_obs(p) = max(mEn.Coefficients.Estimate(2), 0.00);
end
fprintf('Extracted individual parameters for N = %d observers (mean beta_V1 = %.3f ± %.3f).\n', ...
    nParts, mean(betaV1_obs), std(betaV1_obs)/sqrt(nParts));

%% 2. 8-Gabor Display Geometry & Spatial Readout Kernels
nGabors      = 8;
imSize       = 256;
ringRadius   = 72;                  % px (inter-item spacing ~ 55 px)
gabSD        = 8;                   % Gabor Gaussian envelope SD (px)
gabFreq      = 0.14;                % cycles/px
ensSDs_deg   = [10, 20, 35];        % Narrow, Medium, Wide ensemble SD (deg)
nTrialsPerSD = 200;                 % If 200 trials/SD/observer = 18,000 simulated trials total

posUncertPx  = [0, 28, 36];         % 0 = veridical; 28 & 36 px = modest spatial WM uncertainty
mainPosIdx   = 3;                   % highlight 36 px in Panels A & B
nPos         = numel(posUncertPx);
nSD          = numel(ensSDs_deg);

phiDeg  = (0:1:179)';
phiRad2 = deg2rad(phiDeg) * 2;

[X, Y] = meshgrid(1:imSize, 1:imSize);
cx = imSize/2; cy = imSize/2;
angLoc = linspace(0, 2*pi, nGabors + 1); angLoc(end) = [];
locX = cx + ringRadius * cos(angLoc);
locY = cy + ringRadius * sin(angLoc);

E_maps = zeros(imSize, imSize, nGabors);
for j = 1:nGabors
    E_maps(:,:,j) = exp(-((X - locX(j)).^2 + (Y - locY(j)).^2) / (2 * gabSD^2));
end

overlapW = zeros(nGabors, nGabors, nPos);
for pIdx = 1:nPos
    sigP = posUncertPx(pIdx);
    for cue = 1:nGabors
        if sigP == 0
            W_cue = E_maps(:,:,cue);
        else
            W_cue = exp(-((X - locX(cue)).^2 + (Y - locY(cue)).^2) / (2 * sigP^2));
        end
        for j = 1:nGabors
            overlapW(cue, j, pIdx) = sum(W_cue .* E_maps(:,:,j), 'all');
        end
        overlapW(cue, :, pIdx) = overlapW(cue, :, pIdx) / overlapW(cue, cue, pIdx);
    end
end

%% 3. Simulate All 30 Observers Independently
edges   = [0, 10, 20, 30, 45, 65];
centers = (edges(1:end-1) + edges(2:end)) / 2;
nBins   = numel(centers);

% Per-observer aggregated metrics (N = 30)
obsMeanBiasV1  = zeros(nParts, nSD, nPos);
obsMeanBiasEns = zeros(nParts, nSD);
obsSlopeV1     = zeros(nParts, nSD, nPos);
obsSlopeEns    = zeros(nParts, nSD);
obsBinBiasV1   = nan(nParts, nSD, nBins);
obsBinBiasEns  = nan(nParts, nSD, nBins);

exDisplay = [];
fprintf('Simulating %d observers x %d SDs x %d trials (%d total trials)...\n', ...
    nParts, nSD, nTrialsPerSD, nParts * nSD * nTrialsPerSD);

for p = 1:nParts
    % Observer's own V1 decision gain & relative ensemble weight
    v1GainP = betaV1_obs(p) * 30;
    relEnsP = betaEns_obs(p) / max(betaV1_obs(p), 0.1);
    % Item precision scales inversely with observer's V1 sensitivity
    sigmaItemP = min(max(6 / max(betaV1_obs(p), 0.15), 8), 25);

    for sIdx = 1:nSD
        ensSD = ensSDs_deg(sIdx);
        absOff  = zeros(nTrialsPerSD, 1);
        bV1_tr  = zeros(nTrialsPerSD, nPos);
        bEns_tr = zeros(nTrialsPerSD, 1);

        for tr = 1:nTrialsPerSD
            muTrue   = rand * 180;
            orisDeg  = mod(muTrue + randn(1, nGabors) * ensSD, 180);
            orisRad2 = deg2rad(orisDeg) * 2;

            meanVec      = mean(exp(1i * orisRad2));
            ensMeanDeg   = mod(rad2deg(angle(meanVec) / 2), 180);
            ensSDEmp_deg = rad2deg(sqrt(max(-2 * log(max(abs(meanVec), 1e-6)), 0)) / 2);

            cue     = randi(nGabors);
            targDeg = orisDeg(cue);
            dTarg   = circDiff180(targDeg, ensMeanDeg);
            absOff(tr) = abs(dTarg);

            if p == 1 && sIdx == 1 && tr == 1
                img = zeros(imSize, imSize);
                for j = 1:nGabors
                    th = deg2rad(orisDeg(j));
                    xr = (X - locX(j)) * cos(th) + (Y - locY(j)) * sin(th);
                    img = img + E_maps(:,:,j) .* cos(2*pi*gabFreq * xr);
                end
                exDisplay = img;
            end

            % Winning V1_Spatial observer (using observer p's own beta_V1)
            for pIdx = 1:nPos
                w = overlapW(cue, :, pIdx);
                S_phi  = (w * cos(orisRad2' - phiRad2'))' / (norm(w) + eps);
                logits = v1GainP * (S_phi - max(S_phi));
                probs  = exp(logits);
                probs  = probs / sum(probs);
                respDeg = phiDeg(find(cumsum(probs) >= rand, 1, 'first'));
                bV1_tr(tr, pIdx) = -sign(dTarg) * circDiff180(respDeg, targDeg);
            end

            % Explicit Ensemble-Prior observer (using observer p's own beta_Ens / beta_V1)
            wEns = min(max(relEnsP * (sigmaItemP^2 / (sigmaItemP^2 + ensSDEmp_deg^2)), 0), 0.8);
            postVec = (1 - wEns) * exp(1i * deg2rad(targDeg)*2) + wEns * exp(1i * deg2rad(ensMeanDeg)*2);
            postMeanDeg = mod(rad2deg(angle(postVec) / 2), 180);
            respEnsDeg  = mod(postMeanDeg + randn * sigmaItemP, 180);
            bEns_tr(tr) = -sign(dTarg) * circDiff180(respEnsDeg, targDeg);
        end

        % Summarise observer p for this SD condition
        m = (absOff >= 5) & (absOff <= 60);
        obsMeanBiasV1(p, sIdx, :) = mean(bV1_tr(m, :), 1);
        obsMeanBiasEns(p, sIdx)   = mean(bEns_tr(m));

        for pIdx = 1:nPos
            c = polyfit(absOff(m), bV1_tr(m, pIdx), 1);
            obsSlopeV1(p, sIdx, pIdx) = c(1);
        end
        cE = polyfit(absOff(m), bEns_tr(m), 1);
        obsSlopeEns(p, sIdx) = cE(1);

        for bIdx = 1:nBins
            bm = (absOff >= edges(bIdx)) & (absOff < edges(bIdx+1));
            if sum(bm) >= 3
                obsBinBiasV1(p, sIdx, bIdx)  = mean(bV1_tr(bm, mainPosIdx));
                obsBinBiasEns(p, sIdx, bIdx) = mean(bEns_tr(bm));
            end
        end
    end
end

%% 4. Print Group-Level Results Across N = 30 Simulated Participants & Save CSV
summaryTbl = table();
fprintf('\n=== SIMULATED GROUP RESULTS ACROSS N = %d PARTICIPANTS (Mean ± SE across observers) ===\n', nParts);
fprintf('%-8s | %-17s | %-17s | %-17s | %-17s\n', ...
    'Ens SD', 'V1 (0 px) Slope', 'V1 (28 px) Slope', 'V1 (36 px) Slope', 'Ensemble Slope');
fprintf('%s\n', repmat('-', 1, 88));

for sIdx = 1:nSD
    sd = ensSDs_deg(sIdx);
    mS_V1  = squeeze(mean(obsSlopeV1(:, sIdx, :), 1));
    seS_V1 = squeeze(std(obsSlopeV1(:, sIdx, :), 0, 1)) / sqrt(nParts);
    mS_En  = mean(obsSlopeEns(:, sIdx));
    seS_En = std(obsSlopeEns(:, sIdx)) / sqrt(nParts);

    mB_V1  = squeeze(mean(obsMeanBiasV1(:, sIdx, :), 1));
    seB_V1 = squeeze(std(obsMeanBiasV1(:, sIdx, :), 0, 1)) / sqrt(nParts);
    mB_En  = mean(obsMeanBiasEns(:, sIdx));
    seB_En = std(obsMeanBiasEns(:, sIdx)) / sqrt(nParts);

    fprintf('%3d deg  | %+5.3f ± %.3f     | %+5.3f ± %.3f     | %+5.3f ± %.3f     | %+5.3f ± %.3f\n', ...
        sd, mS_V1(1), seS_V1(1), mS_V1(2), seS_V1(2), mS_V1(3), seS_V1(3), mS_En, seS_En);

    row = table(sd, mB_V1(1), seB_V1(1), mB_V1(2), seB_V1(2), mB_V1(3), seB_V1(3), mB_En, seB_En, ...
        mS_V1(1), seS_V1(1), mS_V1(2), seS_V1(2), mS_V1(3), seS_V1(3), mS_En, seS_En, ...
        'VariableNames', {'EnsembleSD_deg', ...
        'Bias_V1_0px', 'SE_Bias_0px', 'Bias_V1_28px', 'SE_Bias_28px', 'Bias_V1_36px', 'SE_Bias_36px', ...
        'Bias_Ensemble', 'SE_Bias_Ensemble', ...
        'Slope_V1_0px', 'SE_Slope_0px', 'Slope_V1_28px', 'SE_Slope_28px', 'Slope_V1_36px', 'SE_Slope_36px', ...
        'Slope_Ensemble', 'SE_Slope_Ensemble'});
    summaryTbl = [summaryTbl; row]; %#ok<AGROW>
end

% Paired test across N = 30 observers: Narrow (10 deg) vs Wide (35 deg) slope
[~, pSlope, ~, stSlope] = ttest(obsSlopeV1(:, 1, mainPosIdx), obsSlopeV1(:, 3, mainPosIdx));
nPosBias = sum(mean(obsMeanBiasV1(:, :, mainPosIdx), 2) > 0);
nSteeper = sum(obsSlopeV1(:, 1, mainPosIdx) > obsSlopeV1(:, 3, mainPosIdx));
fprintf('\nObserver-level checks (V1_Spatial, sigma_pos = %d px, N = %d):\n', posUncertPx(mainPosIdx), nParts);
fprintf('  Positive attraction bias toward mean: %d / %d observers\n', nPosBias, nParts);
fprintf('  Steeper bias at Narrow (10°) vs Wide (35°) SD: %d / %d observers (paired t(%d) = %.2f, p = %.2e)\n', ...
    nSteeper, nParts, stSlope.df, stSlope.tstat, pSlope);

writetable(summaryTbl, 'simulated_ensemble_bias_summary.csv');

%% 5. Generate Publication-Ready 3-Panel Figure (Fixed Legend + N=30 SE Error Bars)
f = figure('Color', 'w', 'Position', [80, 120, 1220, 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel A: Example 8-Gabor Display + Spatial Readout Kernel
nexttile; hold on;
imagesc(exDisplay, [-1, 1]); colormap(gca, gray); axis image off;
thC = linspace(0, 2*pi, 100);
plot(locX(1) + 18*cos(thC), locY(1) + 18*sin(thC), '-', 'Color', [0.95 0.75 0.20], 'LineWidth', 2.0);
plot(locX(1) + posUncertPx(mainPosIdx)*cos(thC), locY(1) + posUncertPx(mainPosIdx)*sin(thC), ...
    '--', 'Color', [0.22 0.60 0.90], 'LineWidth', 1.8);
plot(locX(1) + 2*posUncertPx(mainPosIdx)*cos(thC), locY(1) + 2*posUncertPx(mainPosIdx)*sin(thC), ...
    ':', 'Color', [0.22 0.60 0.90], 'LineWidth', 1.3);
title('A. 8-Gabor Display & Spatial Readout', 'FontSize', 11, 'FontWeight', 'bold');
text(8, 242, 'Yellow: cued target | Blue: 1\sigma & 2\sigma memory kernel', ...
    'Color', 'w', 'FontSize', 8, 'BackgroundColor', [0 0 0 0.55]);

% Panel B: Group Mean Attraction Bias (Mean ± 1 SE across N = 30 simulated participants)
nexttile; hold on;
sdCols = [0.15 0.45 0.78; 0.90 0.55 0.15; 0.80 0.25 0.35];
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

hLegB = gobjects(nSD + 1, 1);
for sIdx = 1:nSD
    mB  = nan(nBins, 1); seB = nan(nBins, 1);
    mE  = nan(nBins, 1);
    for bIdx = 1:nBins
        vB = squeeze(obsBinBiasV1(:, sIdx, bIdx));  vB = vB(~isnan(vB));
        vE = squeeze(obsBinBiasEns(:, sIdx, bIdx)); vE = vE(~isnan(vE));
        if numel(vB) >= 5
            mB(bIdx)  = mean(vB);
            seB(bIdx) = std(vB) / sqrt(numel(vB));
            mE(bIdx)  = mean(vE);
        end
    end
    hLegB(sIdx) = errorbar(centers, mB, seB, '-o', 'Color', sdCols(sIdx,:), ...
        'MarkerFaceColor', sdCols(sIdx,:), 'LineWidth', 2.0, 'MarkerSize', 5.5, 'CapSize', 3);
    plot(centers, mE, '--', 'Color', sdCols(sIdx,:), 'LineWidth', 1.4, 'HandleVisibility', 'off');
end
hLegB(nSD+1) = plot(nan, nan, 'k--', 'LineWidth', 1.4);
xlabel('|Target \theta - Ensemble Mean| (deg)', 'FontSize', 10.5);
ylabel('Response Bias Toward Mean (deg)', 'FontSize', 10.5);
legend(hLegB, {'V1 Spatial (SD = 10°)', 'V1 Spatial (SD = 20°)', 'V1 Spatial (SD = 35°)', 'Explicit Ensemble (dashed)'}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', 8.5);
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
title('B. Simulated Group Bias (N = 30)', 'FontSize', 11, 'FontWeight', 'bold');

% Panel C: Group Mean Bias Slope vs Ensemble SD (Fixed Legend + N = 30 SE bars)
nexttile; hold on;
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
hC = gobjects(4, 1);
hC(1) = errorbar(ensSDs_deg, summaryTbl.Slope_V1_0px,  summaryTbl.SE_Slope_0px,  '-o', ...
    'Color', [0.55 0.55 0.55], 'MarkerFaceColor', [0.55 0.55 0.55], 'LineWidth', 1.8, 'CapSize', 3);
hC(2) = errorbar(ensSDs_deg, summaryTbl.Slope_V1_28px, summaryTbl.SE_Slope_28px, '-s', ...
    'Color', [0.35 0.65 0.88], 'MarkerFaceColor', [0.35 0.65 0.88], 'LineWidth', 2.0, 'CapSize', 3);
hC(3) = errorbar(ensSDs_deg, summaryTbl.Slope_V1_36px, summaryTbl.SE_Slope_36px, '-d', ...
    'Color', [0.12 0.38 0.72], 'MarkerFaceColor', [0.12 0.38 0.72], 'LineWidth', 2.2, 'CapSize', 3);
hC(4) = errorbar(ensSDs_deg, summaryTbl.Slope_Ensemble, summaryTbl.SE_Slope_Ensemble, '--^', ...
    'Color', [0.85 0.35 0.22], 'MarkerFaceColor', [0.85 0.35 0.22], 'LineWidth', 2.0, 'CapSize', 3);
xlabel('Ensemble SD (deg)', 'FontSize', 10.5);
ylabel('Bias Slope (deg bias / deg offset)', 'FontSize', 10.5);
set(gca, 'XTick', ensSDs_deg, 'XLim', [5, 40], 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
legend(hC, {'V1 Spatial (\sigma_{pos} = 0 px)', 'V1 Spatial (\sigma_{pos} = 28 px)', ...
            'V1 Spatial (\sigma_{pos} = 36 px)', 'Explicit Ensemble Model'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 8.5);
title('C. Wider Ensemble SD Reduces Bias', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f, 'fig_simulated_ensemble_bias.png', 'Resolution', 300);
exportgraphics(f, 'fig_simulated_ensemble_bias.pdf', 'ContentType', 'vector');
fprintf('\nSaved fig_simulated_ensemble_bias.png/.pdf and simulated_ensemble_bias_summary.csv.\n');

%% Helper
function d = circDiff180(a, b)
    d = mod(a - b + 90, 180) - 90;
end
