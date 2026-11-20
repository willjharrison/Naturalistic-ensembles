%% SCModel_clean_compare.m
% Fast, self-contained Experiment 1 model pipeline & comparison.
% Uses existing precomputed stimulus decompositions in:
%   ../Matlab/Output/
%   ../Matlab/Analysis/SCmodel/modIms/
%
% Optimised for speed:
%   - Loads each of the 12,000 modIms/ files ONCE (via parfor) instead of 76x
%   - Applies each spatial blur ONCE per trial (4 blurs instead of 76)
%   - Caches extracted predictors to Working/cleanModelPredictors.mat
%   - Runs 5-fold cross-validated GLMMs in parallel (parfor)
%
% Evaluates 6 candidate models across a 4x4 grid of orientation noise & spatial blur:
%   1. Ensemble_Weighted   : energy-weighted circular mean + SD mismatch (spatially invariant)
%   2. Ensemble_Unweighted : unweighted circular mean + SD mismatch (spatially invariant)
%   3. Spatial_EnergyOnly  : spatial energy dot-product match
%   4. Spatial_OriOnly     : spatially localised orientation error
%   5. Spatial_Combined    : spatial energy + localised orientation
%   6. Full_Combined       : spatial (energy + ori) + weighted ensemble stats
%
% All predictors follow a consistent sign convention:
%   positive value => favours the RIGHT option (targInt = 1, resp = 1).
%
% Outputs saved in Working/: cleanModelPredictors.mat, cleanModelResults.mat, cleanModelSummary.csv.

clear; close all; clc;
tStart = tic;

%% 1. Paths & Grid Configuration
outputDir = '../Matlab/Output';
modImsDir = '../Matlab/Analysis/SCmodel/modIms';

files = dir(fullfile(outputDir, '*.mat'));
files = files(~ismember({files.name}, {'.', '..', '.DS_Store'}));
nFiles = length(files);
assert(nFiles == 20, 'Expected 20 Experiment 1 files in %s, found %d', outputDir, nFiles);

oriNoiseGrid = [0, 0.05, 0.10, 0.20]; % fraction of pi
posBlurGrid  = [0, 2, 4, 8];          % Gaussian blur sigma (pixels)
nOri  = numel(oriNoiseGrid);
nBlur = numel(posBlurGrid);
nReps = 5;

% Build list of configurations: [oriIdx, blurIdx, rep]
cfgList = [];
for o = 1:nOri
    for b = 1:nBlur
        repsHere = nReps;
        if o == 1 && b == 1
            repsHere = 1; % noiseless condition is deterministic
        end
        for r = 1:repsHere
            cfgList = [cfgList; o, b, r]; %#ok<AGROW>
        end
    end
end
nConfigs = size(cfgList, 1);
fprintf('Configured %d (oriNoise x posBlur x rep) evaluations.\n', nConfigs);

%% 2. Load Trial Metadata Once
fprintf('Loading Experiment 1 metadata (%d files)...\n', nFiles);
nTrialsPerFile = 600;
totalTrials = nFiles * nTrialsPerFile;

participant = zeros(totalTrials, 1);
block       = zeros(totalTrials, 1); % 1 = aligned, 0 = randomised
nGabors     = zeros(totalTrials, 1);
sameDiff    = zeros(totalTrials, 1);
targInt     = zeros(totalTrials, 1); % 0 = left, 1 = right
resp        = zeros(totalTrials, 1); % 0 = left, 1 = right
targPosRaw  = zeros(totalTrials, 1);

for f = 1:nFiles
    ff = load(fullfile(files(f).folder, files(f).name), 'task', 'dv', 'subject_name');
    idx = (1:nTrialsPerFile) + nTrialsPerFile * (f - 1);
    whichBlock = double(contains(files(f).name, 'align'));

    participant(idx) = str2double(ff.subject_name);
    block(idx)       = whichBlock;
    nGabors(idx)     = ff.task.nGabors;
    sameDiff(idx)    = ff.task.sameDiff;
    targPosRaw(idx)  = ff.task.targPos;
    targInt(idx)     = ff.task.targPos - 1;
    resp(idx)        = ff.dv.chosenSquare - 1;
end

%% 3. Single-Pass Feature Extraction Across All 12,000 Trials (parfor)
if exist('cleanModelPredictors.mat', 'file')
    fprintf('Found cached cleanModelPredictors.mat — loading...\n');
    load('cleanModelPredictors.mat', 'leftMatchMat', 'rightMatchMat', 'leftOriMat', 'rightOriMat', ...
        'absMuMat', 'absSDMat', 'absMuUnwMat', 'absSDUnwMat');
else
    fprintf('Extracting predictors in a single parallel pass over %d trials...\n', totalTrials);
    leftMatchMat  = zeros(totalTrials, nConfigs);
    rightMatchMat = zeros(totalTrials, nConfigs);
    leftOriMat    = zeros(totalTrials, nConfigs);
    rightOriMat   = zeros(totalTrials, nConfigs);
    absMuMat      = zeros(totalTrials, nConfigs);
    absSDMat      = zeros(totalTrials, nConfigs);
    absMuUnwMat   = zeros(totalTrials, nConfigs);
    absSDUnwMat   = zeros(totalTrials, nConfigs);

    parfor tIdx = 1:totalTrials
        rng(1042 + tIdx);
        fIdx = ceil(tIdx / nTrialsPerFile);
        trialIdx = tIdx - (fIdx - 1) * nTrialsPerFile;
        stimFile = fullfile(modImsDir, sprintf('gabIm_trial_%03d_fileLoop_%02d.mat', trialIdx, fIdx));
        lf = load(stimFile, 'gabEn', 'filtEn', 'distEn', 'imOris', 'distOris', 'gabOris');

        % Phase-invariant orientations in [0, 2*pi)
        imOris0 = lf.imOris;     imOris0(imOris0 < 0)     = imOris0(imOris0 < 0)     + pi; imOris0   = wrap(imOris0 * 2);
        distOris0 = lf.distOris; distOris0(distOris0 < 0) = distOris0(distOris0 < 0) + pi; distOris0 = wrap(distOris0 * 2);
        gabOris0 = lf.gabOris - pi; gabOris0(gabOris0 < 0) = gabOris0(gabOris0 < 0) + pi; gabOris0  = wrap(gabOris0 * 2);

        validGab = ~isnan(gabOris0(:));
        gabOris0V = gabOris0(validGab);
        imOris0V  = imOris0(validGab);
        distOris0V = distOris0(validGab);

        % Precompute the 4 non-negative spatial blur levels ONCE per trial
        filtEnBase = max(lf.filtEn, 0);
        distEnBase = max(lf.distEn, 0);
        gabEnBase  = max(lf.gabEn,  0);

        filtEnB = cell(nBlur, 1);
        distEnB = cell(nBlur, 1);
        gabEnB  = cell(nBlur, 1);
        gabEnBV = cell(nBlur, 1);
        targMatchB = zeros(nBlur, 1);
        distMatchB = zeros(nBlur, 1);

        for bIdx = 1:nBlur
            sig = posBlurGrid(bIdx);
            if sig > 0
                fE = max(imgaussfilt(filtEnBase, sig), 0);
                dE = max(imgaussfilt(distEnBase, sig), 0);
                gE = max(imgaussfilt(gabEnBase,  sig), 0);
            else
                fE = filtEnBase;
                dE = distEnBase;
                gE = gabEnBase;
            end
            filtEnB{bIdx} = fE(:);
            distEnB{bIdx} = dE(:);
            gabEnB{bIdx}  = gE(:);
            gabEnBV{bIdx} = gE(validGab);
            targMatchB(bIdx) = dot(fE(:), gE(:));
            distMatchB(bIdx) = dot(dE(:), gE(:));
        end

        % Allocate row buffers
        rLeftMatch  = zeros(1, nConfigs);
        rRightMatch = zeros(1, nConfigs);
        rLeftOri    = zeros(1, nConfigs);
        rRightOri   = zeros(1, nConfigs);
        rAbsMu      = zeros(1, nConfigs);
        rAbsSD      = zeros(1, nConfigs);
        rAbsMuUnw   = zeros(1, nConfigs);
        rAbsSDUnw   = zeros(1, nConfigs);

        tPos = targPosRaw(tIdx);

        for cfg = 1:nConfigs
            oCfg = cfgList(cfg, 1);
            bCfg = cfgList(cfg, 2);
            oNoise = oriNoiseGrid(oCfg) * pi;

            % Spatial energy match (depends only on blur)
            if tPos == 1
                rLeftMatch(cfg)  = targMatchB(bCfg);
                rRightMatch(cfg) = distMatchB(bCfg);
            else
                rLeftMatch(cfg)  = distMatchB(bCfg);
                rRightMatch(cfg) = targMatchB(bCfg);
            end

            % Apply orientation noise if > 0
            if oNoise > 0
                imO   = wrap(imOris0   + randn(size(imOris0))   * oNoise);
                distO = wrap(distOris0 + randn(size(distOris0)) * oNoise);
                gabOV = wrap(gabOris0V + randn(size(gabOris0V)) * oNoise);
                imOV   = imO(validGab);
                distOV = distO(validGab);
            else
                imO    = imOris0;
                distO  = distOris0;
                gabOV  = gabOris0V;
                imOV   = imOris0V;
                distOV = distOris0V;
            end

            % Spatially localised orientation error
            tOriErr = sqrt(sum(circ_dist(imOV,   gabOV).^2));
            dOriErr = sqrt(sum(circ_dist(distOV, gabOV).^2));
            if tPos == 1
                rLeftOri(cfg)  = tOriErr;
                rRightOri(cfg) = dOriErr;
            else
                rLeftOri(cfg)  = dOriErr;
                rRightOri(cfg) = tOriErr;
            end

            % Energy-weighted ensemble stats
            fEV = filtEnB{bCfg};
            dEV = distEnB{bCfg};
            gEV = gabEnBV{bCfg};
            [~, gSD] = circ_std(gabOV, gEV);
            [~, tSD] = circ_std(imO(:), fEV);
            [~, dSD] = circ_std(distO(:), dEV);
            gMu = circ_mean(gabOV, gEV);
            tMu = circ_mean(imO(:), fEV);
            dMu = circ_mean(distO(:), dEV);

            tMuErr = abs(circ_dist(gMu, tMu));
            dMuErr = abs(circ_dist(gMu, dMu));
            tSDErr = abs(gSD - tSD);
            dSDErr = abs(gSD - dSD);

            if tPos == 1
                rAbsMu(cfg) = tMuErr - dMuErr;
                rAbsSD(cfg) = tSDErr - dSDErr;
            else
                rAbsMu(cfg) = dMuErr - tMuErr;
                rAbsSD(cfg) = dSDErr - tSDErr;
            end

            % Unweighted ensemble stats
            [~, gSDu] = circ_std(gabOV);
            [~, tSDu] = circ_std(imO(:));
            [~, dSDu] = circ_std(distO(:));
            gMuu = circ_mean(gabOV);
            tMuu = circ_mean(imO(:));
            dMuu = circ_mean(distO(:));

            tMuEu = abs(circ_dist(gMuu, tMuu));
            dMuEu = abs(circ_dist(gMuu, dMuu));
            tSDEu = abs(gSDu - tSDu);
            dSDEu = abs(gSDu - dSDu);

            if tPos == 1
                rAbsMuUnw(cfg) = tMuEu - dMuEu;
                rAbsSDUnw(cfg) = tSDEu - dSDEu;
            else
                rAbsMuUnw(cfg) = dMuEu - tMuEu;
                rAbsSDUnw(cfg) = dSDEu - tSDEu;
            end
        end

        leftMatchMat(tIdx, :)  = rLeftMatch;
        rightMatchMat(tIdx, :) = rRightMatch;
        leftOriMat(tIdx, :)    = rLeftOri;
        rightOriMat(tIdx, :)   = rRightOri;
        absMuMat(tIdx, :)      = rAbsMu;
        absSDMat(tIdx, :)      = rAbsSD;
        absMuUnwMat(tIdx, :)   = rAbsMuUnw;
        absSDUnwMat(tIdx, :)   = rAbsSDUnw;
    end

    fprintf('Feature extraction finished in %.1f mins. Saving cache...\n', toc(tStart)/60);
    save('cleanModelPredictors.mat', 'leftMatchMat', 'rightMatchMat', 'leftOriMat', 'rightOriMat', ...
        'absMuMat', 'absSDMat', 'absMuUnwMat', 'absSDUnwMat', 'cfgList', 'oriNoiseGrid', 'posBlurGrid', '-v7.3');
end

%% 4. Define Candidate Models & Run Cross-Validated Fits (parfor)
modelNames = { ...
    'Ensemble_Weighted', ...
    'Ensemble_Unweighted', ...
    'Spatial_EnergyOnly', ...
    'Spatial_OriOnly', ...
    'Spatial_Combined', ...
    'Full_Combined'};

modelRHS = { ...
    'absMuDiff + absSDDiff', ...
    'absMuUnwDiff + absSDUnwDiff', ...
    'matchDiff', ...
    'oriDiff', ...
    'matchDiff + oriDiff', ...
    'matchDiff + oriDiff + absMuDiff + absSDDiff'};

nModels = numel(modelNames);
kFolds = 5;

fprintf('Fitting %d models across %d configurations in parallel...\n', nModels, nConfigs);
cfgMetrics = zeros(nConfigs, nModels, 6); % [choiceLogLik, choiceAcc, idealAlign, idealRand, groupRMSE, indivRMSE]

parfor cfg = 1:nConfigs
    rng(2042 + cfg);

    matchDiff    = zscore(log(max(rightMatchMat(:, cfg), eps)) - log(max(leftMatchMat(:, cfg), eps)));
    oriDiff      = zscore(log(max(leftOriMat(:, cfg), eps))    - log(max(rightOriMat(:, cfg), eps)));
    absMuDiffZ   = zscore(absMuMat(:, cfg));
    absSDDiffZ   = zscore(absSDMat(:, cfg));
    absMuUnwZ    = zscore(absMuUnwMat(:, cfg));
    absSDUnwZ    = zscore(absSDUnwMat(:, cfg));

    T = table(participant, block, sameDiff, nGabors, targInt, resp, ...
        matchDiff, oriDiff, absMuDiffZ, absSDDiffZ, absMuUnwZ, absSDUnwZ, ...
        'VariableNames', {'participant', 'block', 'sameDiff', 'nGabors', 'targInt', 'resp', ...
        'matchDiff', 'oriDiff', 'absMuDiff', 'absSDDiff', 'absMuUnwDiff', 'absSDUnwDiff'});

    assert(~any(any(ismissing(T))), 'NaN values detected in config %d!', cfg);

    permIdx = randperm(totalTrials);
    Tshuf = T(permIdx, :);
    localM = zeros(nModels, 6);

    for m = 1:nModels
        rhs = modelRHS{m};
        humanSpec = sprintf('resp ~ %s + (1 | participant)', rhs);
        idealSpec = sprintf('targInt ~ %s', rhs);

        predHumanProb = zeros(totalTrials, 1);
        predIdealProb = zeros(totalTrials, 1);

        for fold = 1:kFolds
            testMask = ((1:totalTrials) <= totalTrials/kFolds * fold) & ...
                       ((1:totalTrials) >  totalTrials/kFolds * (fold - 1));
            trainTbl = Tshuf(~testMask, :);
            testTbl  = Tshuf(testMask, :);

            mHuman = fitglme(trainTbl, humanSpec, 'Distribution', 'Binomial');
            predHumanProb(testMask) = predict(mHuman, testTbl);

            mIdeal = fitglm(trainTbl, idealSpec, 'Distribution', 'Binomial');
            predIdealProb(testMask) = predict(mIdeal, testTbl);
        end

        pClip = min(max(predHumanProb, 1e-6), 1 - 1e-6);
        choiceLL  = sum(Tshuf.resp .* log(pClip) + (1 - Tshuf.resp) .* log(1 - pClip));
        choiceAcc = mean(round(predHumanProb) == Tshuf.resp);

        idealCorrect = double(round(predIdealProb) == Tshuf.targInt);
        idealAlign = mean(idealCorrect(Tshuf.block == 1));
        idealRand  = mean(idealCorrect(Tshuf.block == 0));

        humanCorrect = double(Tshuf.resp == Tshuf.targInt);
        modelCorrect = double(round(predHumanProb) == Tshuf.targInt);

        indHuman = qStats(humanCorrect, [Tshuf.participant, Tshuf.block, Tshuf.sameDiff, Tshuf.nGabors], 2);
        indModel = qStats(modelCorrect, [Tshuf.participant, Tshuf.block, Tshuf.sameDiff, Tshuf.nGabors], 2);
        indivRMSE = sqrt(mean((indModel(:, 5) - indHuman(:, 5)).^2));

        grpHuman = qStats(indHuman(:, 5), indHuman(:, 2:4), 2);
        grpModel = qStats(indModel(:, 5), indModel(:, 2:4), 2);
        groupRMSE = sqrt(mean((grpModel(:, 4) - grpHuman(:, 4)).^2));

        localM(m, :) = [choiceLL, choiceAcc, idealAlign, idealRand, groupRMSE, indivRMSE];
    end

    cfgMetrics(cfg, :, :) = localM;
end

%% 5. Aggregate Over Repetitions & Save Results
summaryRows = table();
for o = 1:nOri
    for b = 1:nBlur
        mask = (cfgList(:, 1) == o) & (cfgList(:, 2) == b);
        meanM = squeeze(mean(cfgMetrics(mask, :, :), 1)); % [nModels x 6]
        for m = 1:nModels
            row = table(oriNoiseGrid(o), posBlurGrid(b), modelNames(m), ...
                meanM(m, 1), meanM(m, 2), meanM(m, 3), meanM(m, 4), meanM(m, 5), meanM(m, 6), ...
                'VariableNames', {'OriNoisePi', 'PosBlurPx', 'Model', ...
                'ChoiceLogLik', 'ChoiceAcc', 'IdealAcc_Align', 'IdealAcc_Rand', 'GroupRMSE', 'IndivRMSE'});
            summaryRows = [summaryRows; row]; %#ok<AGROW>
        end
    end
end

save('cleanModelResults.mat', 'modelNames', 'modelRHS', 'oriNoiseGrid', 'posBlurGrid', 'cfgList', 'cfgMetrics', 'summaryRows');
writetable(summaryRows, 'cleanModelSummary.csv');

fprintf('\nTotal elapsed time: %.1f mins\n', toc(tStart)/60);
fprintf('\n=== Best Choice Log-Likelihood per Model ===\n');
for m = 1:nModels
    sub = summaryRows(strcmp(summaryRows.Model, modelNames{m}), :);
    [~, bestIdx] = max(sub.ChoiceLogLik);
    disp(sub(bestIdx, :));
end
