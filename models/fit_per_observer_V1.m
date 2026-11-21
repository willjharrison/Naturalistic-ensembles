%% fit_per_observer_V1.m
% Fits each candidate observer model INDEPENDENTLY to every individual
% observer (N = 30 across Experiments 1-3; ~1,150 trials/observer) using
% 5-fold out-of-sample cross-validation on cached unifiedPredictors_v2_V1.mat.
%
% Advantages over a single pooled slope:
%   - Allows each observer (and each experiment) their own sensitivity/gain
%   - Enables individual-level model comparison (BIC winner counts, paired
%     tests on ΔLL across observers, and individual β weight distributions)
%   - Zero condition labels (block, nGabors, cond) in any decision rule
%
% Outputs saved in Working/:
%   - perObserver_summary.csv
%   - perObserver_diagnostics.txt
%   - fig_perObserver_model_vs_behaviour.png / .pdf
%   - fig_perObserver_model_comparison.png   / .pdf

clear; close all; clc;

fprintf('Loading cached unifiedPredictors_v2_V1.mat...\n');
load('unifiedPredictors_v2_V1.mat');

%% 1. Rebuild Trial Metadata Across Experiments 1-3
expDirs = { ...
    '../Matlab/Output', ...
    '../Matlab/Experiment 2/Output', ...
    '../Matlab/Experiment 3/Output'};

expIDVec = []; partIDVec = []; blockVec = [];
condVec = []; nGaborsVec = []; targIntVec = []; respVec = [];

for e = 1:3
    fList = dir(fullfile(expDirs{e}, '*.mat'));
    fList = fList(~ismember({fList.name}, {'.', '..', '.DS_Store'}));
    fList = fList(~contains({fList.name}, {'exclusion'}));
    for f = 1:length(fList)
        ff = load(fullfile(fList(f).folder, fList(f).name), 'task', 'dv', 'subject_name');
        nT = ff.task.nTrials;
        expIDVec   = [expIDVec;   repmat(e, nT, 1)]; %#ok<AGROW>
        partIDVec  = [partIDVec;  repmat(e*100 + str2double(ff.subject_name), nT, 1)]; %#ok<AGROW>
        blockVec   = [blockVec;   repmat(double(contains(fList(f).name, 'align')), nT, 1)]; %#ok<AGROW>
        nGaborsVec = [nGaborsVec; ff.task.nGabors(:)]; %#ok<AGROW>
        targIntVec = [targIntVec; ff.task.targPos(:) - 1]; %#ok<AGROW>
        respVec    = [respVec;    ff.dv.chosenSquare(:) - 1]; %#ok<AGROW>
        if isfield(ff.task, 'sameDiff')
            condVec = [condVec; ff.task.sameDiff(:)]; %#ok<AGROW>
        else
            condVec = [condVec; ff.task.thingOrMat(:)]; %#ok<AGROW>
        end
    end
end

expIDVec   = expIDVec(validTrial);
partIDVec  = partIDVec(validTrial);
blockVec   = blockVec(validTrial);
condVec    = condVec(validTrial);
nGaborsVec = nGaborsVec(validTrial);
targIntVec = targIntVec(validTrial);
respVec    = respVec(validTrial);

% Use optimal noise config cfg = 1 (OriNoise = 0, PosBlur = 0)
cfg = 1;
v1Raw   = rightV1Mat(validTrial, cfg) - leftV1Mat(validTrial, cfg);
enRaw   = rightEnMat(validTrial, cfg) - leftEnMat(validTrial, cfg);
muRaw   = absMuMat(validTrial, cfg);
sdRaw   = absSDMat(validTrial, cfg);
muUnw   = absMuUnwMat(validTrial, cfg);
sdUnw   = absSDUnwMat(validTrial, cfg);

%% 2. Fit 6 Models Independently to Each Observer (5-Fold Out-of-Sample CV)
modelNames = { ...
    'Ensemble_Unweighted', ...
    'Ensemble_Weighted', ...
    'Spatial_EnergyOnly', ...
    'V1_Spatial', ...
    'V1_Spatial_Plus_Energy', ...
    'Full_Combined'};

modelRHS = { ...
    'absMuUnwDiff + absSDUnwDiff', ...
    'absMuDiff + absSDDiff', ...
    'enDiff', ...
    'v1Diff', ...
    'v1Diff + enDiff', ...
    'v1Diff + enDiff + absMuDiff + absSDDiff'};

nModels = numel(modelNames);
uParts  = unique(partIDVec);
nParts  = numel(uParts);
kFolds  = 5;

partExp   = zeros(nParts, 1);
partN     = zeros(nParts, 1);
partLL    = zeros(nParts, nModels);
partAIC   = zeros(nParts, nModels);
partBIC   = zeros(nParts, nModels);
partAcc   = zeros(nParts, nModels);
partBetas = zeros(nParts, 4); % [v1Diff, enDiff, absMuDiff, absSDDiff] from Full_Combined

cvPredAll = zeros(numel(respVec), nModels);

fprintf('Fitting %d models independently to %d observers (5-fold CV)...\n', nModels, nParts);
rng(8042);

for p = 1:nParts
    idxP = find(partIDVec == uParts(p));
    nT   = numel(idxP);
    partExp(p) = expIDVec(idxP(1));
    partN(p)   = nT;

    % Standardise predictors within each observer (preserves zero-label design)
    Tp = table(respVec(idxP), targIntVec(idxP), ...
        zscore(v1Raw(idxP)), zscore(enRaw(idxP)), ...
        zscore(muRaw(idxP)), zscore(sdRaw(idxP)), ...
        zscore(muUnw(idxP)), zscore(sdUnw(idxP)), ...
        'VariableNames', {'resp', 'targInt', 'v1Diff', 'enDiff', ...
        'absMuDiff', 'absSDDiff', 'absMuUnwDiff', 'absSDUnwDiff'});

    perm = randperm(nT);
    Tshuf = Tp(perm, :);

    for m = 1:nModels
        rhs  = modelRHS{m};
        spec = sprintf('resp ~ %s', rhs);
        pShuf = zeros(nT, 1);

        for fold = 1:kFolds
            testMask = ((1:nT) <= nT/kFolds * fold) & ((1:nT) > nT/kFolds * (fold - 1));
            mdl = fitglm(Tshuf(~testMask, :), spec, 'Distribution', 'Binomial');
            pShuf(testMask) = predict(mdl, Tshuf(testMask, :));
        end

        pOrig = zeros(nT, 1);
        pOrig(perm) = pShuf;
        cvPredAll(idxP, m) = pOrig;

        pClip = min(max(pOrig, 1e-6), 1 - 1e-6);
        ll = sum(Tp.resp .* log(pClip) + (1 - Tp.resp) .* log(1 - pClip));
        kParams = numel(strsplit(rhs, '+')) + 1; % predictors + intercept

        partLL(p, m)  = ll;
        partAIC(p, m) = 2 * kParams - 2 * ll;
        partBIC(p, m) = kParams * log(nT) - 2 * ll;
        partAcc(p, m) = mean(round(pOrig) == Tp.resp);
    end

    % Fit full model on all observer trials to extract standardised beta weights
    mFullP = fitglm(Tp, 'resp ~ v1Diff + enDiff + absMuDiff + absSDDiff', 'Distribution', 'Binomial');
    partBetas(p, :) = mFullP.Coefficients.Estimate(2:5)';
end

%% 3. Compute Group RMSE & Print Comprehensive Diagnostics
humAcc = double(respVec == targIntVec);
groupRMSE = zeros(nModels, 1);
for m = 1:nModels
    modAcc = double(targIntVec == 1) .* cvPredAll(:, m) + ...
             double(targIntVec == 0) .* (1 - cvPredAll(:, m));
    indH = qStats(humAcc, [partIDVec, expIDVec, blockVec, condVec, nGaborsVec], 2);
    indM = qStats(modAcc, [partIDVec, expIDVec, blockVec, condVec, nGaborsVec], 2);
    grpH = qStats(indH(:, 6), indH(:, 2:5), 2);
    grpM = qStats(indM(:, 6), indM(:, 2:5), 2);
    groupRMSE(m) = sqrt(mean((grpM(:, 5) - grpH(:, 5)).^2));
end

sumLL  = sum(partLL, 1);
sumAIC = sum(partAIC, 1);
sumBIC = sum(partBIC, 1);
meanAcc = mean(partAcc, 1);

[~, bestBICPerPart] = min(partBIC, [], 2);
[~, bestLLPerPart]  = max(partLL,  [], 2);

summaryTbl = table(modelNames', sumLL', sumLL' - min(sumLL), sumAIC', sumBIC', ...
    meanAcc', groupRMSE, ...
    arrayfun(@(m) sum(bestBICPerPart == m), (1:nModels)'), ...
    arrayfun(@(m) sum(bestLLPerPart  == m), (1:nModels)'), ...
    'VariableNames', {'Model', 'TotalCVLogLik', 'DeltaLL', 'TotalAIC', 'TotalBIC', ...
    'MeanChoiceAcc', 'GroupRMSE', 'BIC_Winners_of_30', 'LL_Winners_of_30'});
writetable(summaryTbl, 'perObserver_summary.csv');

diary('perObserver_diagnostics.txt'); diary on;
fprintf('\n=== PER-OBSERVER UNIFIED MODEL SUMMARY (N = 30 observers, 34,743 trials) ===\n');
disp(summaryTbl);

idxV1   = find(strcmp(modelNames, 'V1_Spatial'));
idxV1E  = find(strcmp(modelNames, 'V1_Spatial_Plus_Energy'));
idxEns  = find(strcmp(modelNames, 'Ensemble_Weighted'));
idxFull = find(strcmp(modelNames, 'Full_Combined'));

fprintf('=== PAIRWISE OBSERVER-LEVEL COMPARISONS ===\n');
fprintf('V1_Spatial beats Ensemble_Weighted on CV LogLik: %d / %d observers (mean ΔLL = %+.1f, median = %+.1f)\n', ...
    sum(partLL(:, idxV1) > partLL(:, idxEns)), nParts, ...
    mean(partLL(:, idxV1) - partLL(:, idxEns)), median(partLL(:, idxV1) - partLL(:, idxEns)));
fprintf('V1_Spatial beats Full_Combined on BIC:           %d / %d observers (total ΔBIC = %+.1f)\n', ...
    sum(partBIC(:, idxV1) < partBIC(:, idxFull)), nParts, ...
    sumBIC(idxV1) - sumBIC(idxFull));
fprintf('V1_Spatial beats V1_Spatial_Plus_Energy on BIC:  %d / %d observers (total ΔBIC = %+.1f)\n', ...
    sum(partBIC(:, idxV1) < partBIC(:, idxV1E)), nParts, ...
    sumBIC(idxV1) - sumBIC(idxV1E));

fprintf('\n=== INDIVIDUAL-OBSERVER BETA WEIGHTS (Full_Combined, Mean ± SE across N=30) ===\n');
bNames = {'v1Diff', 'enDiff', 'absMuDiff', 'absSDDiff'};
for k = 1:4
    fprintf('  %-10s: mean = %+.3f ± %.3f | Exp1 = %+.3f | Exp2 = %+.3f | Exp3 = %+.3f\n', ...
        bNames{k}, mean(partBetas(:,k)), std(partBetas(:,k))/sqrt(nParts), ...
        mean(partBetas(partExp==1, k)), mean(partBetas(partExp==2, k)), mean(partBetas(partExp==3, k)));
end

fprintf('\n=== PER-EXPERIMENT TOTAL BIC & WINNER COUNTS ===\n');
for e = 1:3
    mE = (partExp == e);
    fprintf('Exp %d (N=%d): V1_Spatial BIC=%.1f | Full_Combined BIC=%.1f | Ensemble_Weighted BIC=%.1f | V1 wins BIC vs Full: %d/%d\n', ...
        e, sum(mE), sum(partBIC(mE, idxV1)), sum(partBIC(mE, idxFull)), sum(partBIC(mE, idxEns)), ...
        sum(partBIC(mE, idxV1) < partBIC(mE, idxFull)), sum(mE));
end

fprintf('\n=== ALIGNED BLOCK CONDITION MAIN EFFECTS (Human vs Per-Observer V1_Spatial vs Ensemble) ===\n');
v1Acc  = double(targIntVec == 1) .* cvPredAll(:, idxV1)  + double(targIntVec == 0) .* (1 - cvPredAll(:, idxV1));
ensAcc = double(targIntVec == 1) .* cvPredAll(:, idxEns) + double(targIntVec == 0) .* (1 - cvPredAll(:, idxEns));
condLabels = {{'Similar', 'Dissimilar'}, {'Objects', 'Textures'}, {'Objects', 'Textures'}};
for e = 1:3
    uCond = unique(condVec(expIDVec == e));
    fprintf('Experiment %d:\n', e);
    for cIdx = 1:numel(uCond)
        msk = (expIDVec == e) & (blockVec == 1) & (condVec == uCond(cIdx));
        fprintf('  %-12s | Human: %.3f | V1_Spatial: %.3f | Ensemble: %.3f\n', ...
            condLabels{e}{cIdx}, mean(humAcc(msk)), mean(v1Acc(msk)), mean(ensAcc(msk)));
    end
end
diary off;

%% 4. Figure 1: Per-Observer Model vs Behaviour (3x2 Grid)
col1 = [0.95, 0.75, 0.28]; col2 = [0.95, 0.45, 0.78];
col1D = col1 * 0.72;       col2D = col2 * 0.72;
expTitles = {'Exp 1: Objects', 'Exp 2: Objects vs. Textures', 'Exp 3: Gabor Arrays'};
condNames = {{'Similar', 'Dissimilar'}, {'Objects', 'Textures'}, {'Objects', 'Textures'}};

f1 = figure('Color', 'w', 'Position', [80, 60, 820, 860]);
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for e = 1:3
    uGabs = unique(nGaborsVec(expIDVec == e));
    uCond = unique(condVec(expIDVec == e));
    xPos  = 1:numel(uGabs);

    for bIdx = 1:2
        bVal = 2 - bIdx;
        nexttile; hold on;
        plot([0.6, numel(uGabs)+0.4], [0.5, 0.5], 'k:', 'LineWidth', 1.1);

        hLeg = gobjects(4, 1);
        for cIdx = 1:numel(uCond)
            cVal = uCond(cIdx);
            cCol = col1; cColD = col1D;
            if cIdx == 2, cCol = col2; cColD = col2D; end

            mHum = zeros(numel(uGabs), 1); seHum = zeros(numel(uGabs), 1);
            mSp  = zeros(numel(uGabs), 1); mEn   = zeros(numel(uGabs), 1);

            for gIdx = 1:numel(uGabs)
                mask = (expIDVec == e) & (blockVec == bVal) & ...
                       (condVec == cVal) & (nGaborsVec == uGabs(gIdx));
                uP = unique(partIDVec(mask));
                pH = arrayfun(@(p) mean(humAcc(mask & partIDVec == p)), uP);
                pS = arrayfun(@(p) mean(v1Acc(mask  & partIDVec == p)), uP);
                pE = arrayfun(@(p) mean(ensAcc(mask & partIDVec == p)), uP);
                mHum(gIdx) = mean(pH); seHum(gIdx) = std(pH) / sqrt(numel(pH));
                mSp(gIdx)  = mean(pS); mEn(gIdx)   = mean(pE);
            end

            plot(xPos, mSp, '-',  'Color', cCol,  'LineWidth', 2.2);
            plot(xPos, mEn, '--', 'Color', cColD, 'LineWidth', 1.6);
            errorbar(xPos, mHum, seHum, 'o', 'Color', [0.25 0.25 0.25], ...
                'MarkerFaceColor', cCol, 'MarkerEdgeColor', cCol, 'MarkerSize', 6.5, 'CapSize', 3);

            if cIdx == 1
                hLeg(1) = plot(nan, nan, 'o', 'MarkerFaceColor', col1, 'MarkerEdgeColor', col1, 'MarkerSize', 6.5);
                hLeg(2) = plot(nan, nan, 'o', 'MarkerFaceColor', col2, 'MarkerEdgeColor', col2, 'MarkerSize', 6.5);
                hLeg(3) = plot(nan, nan, 'k-',  'LineWidth', 2.2);
                hLeg(4) = plot(nan, nan, 'k--', 'LineWidth', 1.6);
            end
        end

        xlim([0.6, numel(uGabs)+0.4]); ylim([0.35, 1.0]);
        set(gca, 'XTick', xPos, 'XTickLabel', string(uGabs), 'YTick', [0.5, 0.75, 1.0], ...
            'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
        if e == 1
            if bVal == 1, title('Spatially aligned', 'FontSize', 12, 'FontWeight', 'bold');
            else,         title('Randomised',        'FontSize', 12, 'FontWeight', 'bold'); end
        end
        if bVal == 1
            ylabel({expTitles{e}; 'Proportion correct'}, 'FontSize', 10.5, 'FontWeight', 'bold');
            legend(hLeg, {condNames{e}{1}, condNames{e}{2}, 'V1 Spatial', 'Ensemble Weighted'}, ...
                'Location', 'northwest', 'Box', 'off', 'FontSize', 8.5);
        end
        if e == 3, xlabel('Number of Gabors', 'FontSize', 11); end
    end
end

exportgraphics(f1, 'fig_perObserver_model_vs_behaviour.png', 'Resolution', 300);
exportgraphics(f1, 'fig_perObserver_model_vs_behaviour.pdf', 'ContentType', 'vector');

%% 5. Figure 2: Per-Observer Model Comparison (BIC, Individual Betas, Winner Counts)
f2 = figure('Color', 'w', 'Position', [100, 120, 1140, 390]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel A: Total BIC improvement over Ensemble_Unweighted (higher = better)
nexttile; hold on;
dBIC = max(sumBIC) - sumBIC;
[~, ord] = sort(dBIC, 'ascend');
niceNames = strrep(modelNames(ord), '_', ' ');
barh(1:nModels, dBIC(ord), 'FaceColor', [0.22, 0.49, 0.77], 'EdgeColor', 'none');
set(gca, 'YTick', 1:nModels, 'YTickLabel', niceNames, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
xlabel('\Delta Total BIC vs. Worst Model (higher = better)', 'FontSize', 10);
title('A. Per-Observer BIC Comparison (N = 30)', 'FontSize', 11, 'FontWeight', 'bold');

% Panel B: Individual-observer beta weights (dots = observers, bars = mean ± 95% CI)
nexttile; hold on;
xline(0, 'k:', 'LineWidth', 1);
pLabs = {'V1 oriented match', 'Spatial energy', 'Ensemble mean', 'Ensemble SD'};
expCols = [0.22 0.49 0.77; 0.85 0.37 0.25; 0.25 0.65 0.45];
for k = 1:4
    for e = 1:3
        vals = partBetas(partExp == e, k);
        jitter = (rand(size(vals)) - 0.5) * 0.22;
        scatter(vals, k + jitter, 22, expCols(e,:), 'filled', 'MarkerFaceAlpha', 0.55);
    end
    mB = mean(partBetas(:, k));
    ciB = 1.96 * std(partBetas(:, k)) / sqrt(nParts);
    errorbar(mB, k, 0, 0, ciB, ciB, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6.5, 'LineWidth', 1.8, 'CapSize', 5);
end
set(gca, 'YTick', 1:4, 'YTickLabel', pLabs, 'YLim', [0.4, 4.6], 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
xlabel('Individual-Observer \beta (Mean \pm 95% CI)', 'FontSize', 10);
title('B. Predictor Weights Across 30 Observers', 'FontSize', 11, 'FontWeight', 'bold');

% Panel C: Best Model per Observer (BIC winner count out of 30)
nexttile; hold on;
winCounts = arrayfun(@(m) sum(bestBICPerPart == m), 1:nModels);
barh(1:nModels, winCounts(ord), 'FaceColor', [0.22, 0.49, 0.77], 'EdgeColor', 'none');
set(gca, 'YTick', 1:nModels, 'YTickLabel', niceNames, 'TickDir', 'out', 'Box', 'off', 'FontSize', 10);
xlabel('Number of Observers Best Fit (BIC, out of 30)', 'FontSize', 10);
title('C. Winning Model per Observer (BIC)', 'FontSize', 11, 'FontWeight', 'bold');

exportgraphics(f2, 'fig_perObserver_model_comparison.png', 'Resolution', 300);
exportgraphics(f2, 'fig_perObserver_model_comparison.pdf', 'ContentType', 'vector');

fprintf('\nSaved perObserver_summary.csv, perObserver_diagnostics.txt, and fig_perObserver_*.png/.pdf.\n');
