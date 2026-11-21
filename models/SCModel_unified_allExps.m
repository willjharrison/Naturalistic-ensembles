%% SCModel_unified_allExps.m
% Unified, image-driven observer model across Experiments 1, 2, and 3.
%
% Version 2 (V1 Oriented Matched-Filter Observer):
%   1. Replaces normalised mean orientation error (/sum(gE)) and cosine
%      similarity with a principled V1 oriented matched-filter observer:
%        S_V1(I, G) = sum(E_I .* E_G .* cos(2*(theta_I - theta_G))) / norm(E_G)
%      and an unoriented spatial energy matched filter:
%        S_En(I, G) = sum(E_I .* E_G) / norm(E_G)
%      Normalising only by norm(E_G) yields exact SDT sqrt(N) evidence
%      accumulation across set sizes (1 -> 256 Gabors) while preserving
%      similar/dissimilar foil overlap and object/texture sparsity differences.
%   2. Pre-indexes THINGSimages and dtd once into an O(1) lookup map and
%      gracefully skips any Exp 2/3 trial whose source image is missing on disk.
%   3. Fits a single unified model across ALL 3 experiments (30 participants)
%      and ALL conditions with NO condition labels in the decision rule.
%   4. Note: Exp 3 currently uses source photos for response options until
%      Exp 3 response Gabor count is confirmed.
%
% Outputs saved in Working/:
%   unifiedPredictors_v2_V1.mat, unifiedResults_v2_V1.mat, unifiedSummary_v2_V1.csv

clear; close all; clc;
tStart = tic;

%% 1. Paths & Configuration
thingsRoot = '/Users/wharrison/Desktop/removed_domain_49A1B3CB-5689-4D03-A900-D7CC29B8FF62/Research/THINGSimages';
dtdRoot    = '/Users/wharrison/Desktop/removed_domain_49A1B3CB-5689-4D03-A900-D7CC29B8FF62/Research/dtd';
modImsDir  = '../Matlab/Analysis/SCmodel/modIms';

assert(exist(thingsRoot, 'dir') == 7, 'THINGSimages folder not found: %s', thingsRoot);
assert(exist(dtdRoot, 'dir')    == 7, 'DTD folder not found: %s', dtdRoot);

expDirs = { ...
    '../Matlab/Output', ...
    '../Matlab/Experiment 2/Output', ...
    '../Matlab/Experiment 3/Output'};

oriNoiseGrid = [0, 0.10, 0.20, 0.30]; % fraction of pi
posBlurGrid  = [0, 2, 4];             % Gaussian blur sigma (px)
nOri  = numel(oriNoiseGrid);
nBlur = numel(posBlurGrid);
nReps = 2;

cfgList = [];
for o = 1:nOri
    for b = 1:nBlur
        repsHere = nReps;
        if o == 1 && b == 1, repsHere = 1; end
        for r = 1:repsHere
            cfgList = [cfgList; o, b, r]; %#ok<AGROW>
        end
    end
end
nConfigs = size(cfgList, 1);

%% 2. Pre-Index THINGSimages & DTD Once
fprintf('Indexing THINGSimages and dtd folders...\n');
tFiles = dir(fullfile(thingsRoot, '**', '*.*'));
tFiles = tFiles(~[tFiles.isdir]);
dFiles = dir(fullfile(dtdRoot, '**', '*.*'));
dFiles = dFiles(~[dFiles.isdir]);

allKeys = cell(numel(tFiles) + numel(dFiles), 1);
allVals = cell(numel(tFiles) + numel(dFiles), 1);
for i = 1:numel(tFiles)
    allKeys{i} = lower(tFiles(i).name);
    allVals{i} = fullfile(tFiles(i).folder, tFiles(i).name);
end
offset = numel(tFiles);
for i = 1:numel(dFiles)
    allKeys{offset + i} = lower(dFiles(i).name);
    allVals{offset + i} = fullfile(dFiles(i).folder, dFiles(i).name);
end
[uniqKeys, ia] = unique(allKeys, 'stable');
fileMap = containers.Map(uniqKeys, allVals(ia));
fprintf('Indexed %d unique image files (%d THINGS, %d DTD).\n', numel(uniqKeys), numel(tFiles), numel(dFiles));

%% 3. Collect All Valid Files Across Experiments 1-3 & Pre-Resolve Image Paths
fileMeta = struct('expID', {}, 'fileLoop', {}, 'path', {}, 'name', {}, 'nTrials', {}, 'p1List', {}, 'p2List', {}, 'validMask', {});
for e = 1:3
    fList = dir(fullfile(expDirs{e}, '*.mat'));
    fList = fList(~ismember({fList.name}, {'.', '..', '.DS_Store'}));
    fList = fList(~contains({fList.name}, {'exclusion'}));
    keptExp = 0; totalExp = 0;

    for f = 1:length(fList)
        fpath = fullfile(fList(f).folder, fList(f).name);
        if e == 1
            tmp = load(fpath, 'task');
            nT = tmp.task.nTrials;
            p1List = cell(nT, 1);
            p2List = cell(nT, 1);
            vMask  = true(nT, 1);
        else
            tmp = load(fpath, 'task', 'targIms');
            nT = tmp.task.nTrials;
            p1List = cell(nT, 1);
            p2List = cell(nT, 1);
            vMask  = false(nT, 1);
            for tr = 1:nT
                k1 = extractFileKey(tmp.targIms{tr, 1});
                k2 = extractFileKey(tmp.targIms{tr, 2});
                if isKey(fileMap, k1) && isKey(fileMap, k2)
                    p1List{tr} = fileMap(k1);
                    p2List{tr} = fileMap(k2);
                    vMask(tr)  = true;
                end
            end
        end
        keptExp  = keptExp  + sum(vMask);
        totalExp = totalExp + nT;

        entry.expID     = e;
        entry.fileLoop  = f;
        entry.path      = fpath;
        entry.name      = fList(f).name;
        entry.nTrials   = nT;
        entry.p1List    = p1List;
        entry.p2List    = p2List;
        entry.validMask = vMask;
        fileMeta(end+1) = entry; %#ok<AGROW>
    end
    fprintf('Experiment %d: %d sessions | %d / %d valid trials (skipping %d).\n', ...
        e, length(fList), keptExp, totalExp, totalExp - keptExp);
end

nSessions = numel(fileMeta);
rawTotal  = sum([fileMeta.nTrials]);

%% 4. Build Trial-Level Metadata Vectors
expIDVec   = zeros(rawTotal, 1);
partIDVec  = zeros(rawTotal, 1);
blockVec   = zeros(rawTotal, 1);
nGaborsVec = zeros(rawTotal, 1);
condVec    = zeros(rawTotal, 1);
targIntVec = zeros(rawTotal, 1);
respVec    = zeros(rawTotal, 1);
validTrial = false(rawTotal, 1);

cursor = 0;
for s = 1:nSessions
    ff = load(fileMeta(s).path, 'task', 'dv', 'subject_name');
    nT = fileMeta(s).nTrials;
    idx = cursor + (1:nT);
    cursor = cursor + nT;

    expIDVec(idx)   = fileMeta(s).expID;
    partIDVec(idx)  = fileMeta(s).expID * 100 + str2double(ff.subject_name);
    blockVec(idx)   = double(contains(fileMeta(s).name, 'align'));
    nGaborsVec(idx) = ff.task.nGabors;
    targIntVec(idx) = ff.task.targPos - 1;
    respVec(idx)    = ff.dv.chosenSquare - 1;
    validTrial(idx) = fileMeta(s).validMask;

    if isfield(ff.task, 'sameDiff')
        condVec(idx) = ff.task.sameDiff;
    elseif isfield(ff.task, 'thingOrMat')
        condVec(idx) = ff.task.thingOrMat;
    end
end

%% 5. Single-Pass Feature Extraction Across All 60 Sessions (parfor)
cacheFile = 'unifiedPredictors_v2_V1.mat';
if exist(cacheFile, 'file')
    fprintf('Loading cached %s...\n', cacheFile);
    load(cacheFile, 'leftV1Mat', 'rightV1Mat', 'leftEnMat', 'rightEnMat', ...
        'absMuMat', 'absSDMat', 'absMuUnwMat', 'absSDUnwMat', 'validTrial');
else
    fprintf('Extracting V1 matched-filter & ensemble predictors across %d sessions (%d valid trials, %d configs)...\n', ...
        nSessions, sum(validTrial), nConfigs);

    sessLeftV1   = cell(nSessions, 1);
    sessRightV1  = cell(nSessions, 1);
    sessLeftEn   = cell(nSessions, 1);
    sessRightEn  = cell(nSessions, 1);
    sessAbsMu    = cell(nSessions, 1);
    sessAbsSD    = cell(nSessions, 1);
    sessAbsMuUnw = cell(nSessions, 1);
    sessAbsSDUnw = cell(nSessions, 1);

    parfor s = 1:nSessions
        rng(5042 + s);
        meta = fileMeta(s);
        nT = meta.nTrials;

        lV1 = nan(nT, nConfigs); rV1 = nan(nT, nConfigs);
        lEn = nan(nT, nConfigs); rEn = nan(nT, nConfigs);
        aMu = nan(nT, nConfigs); aSD = nan(nT, nConfigs);
        aMuU = nan(nT, nConfigs); aSDU = nan(nT, nConfigs);

        if meta.expID == 1
            ff = load(meta.path, 'task');
            for tr = 1:nT
                stimFile = fullfile(modImsDir, sprintf('gabIm_trial_%03d_fileLoop_%02d.mat', tr, meta.fileLoop));
                lf = load(stimFile, 'gabEn', 'filtEn', 'distEn', 'imOris', 'distOris', 'gabOris');
                [lV1(tr,:), rV1(tr,:), lEn(tr,:), rEn(tr,:), aMu(tr,:), aSD(tr,:), aMuU(tr,:), aSDU(tr,:)] = ...
                    computeV1Predictors(lf.filtEn, lf.distEn, lf.gabEn, lf.imOris, lf.distOris, lf.gabOris, ...
                    ff.task.targPos(tr), oriNoiseGrid, posBlurGrid, cfgList);
            end
        else
            ff = load(meta.path, 'task', 'stimuli', 'filt', 'mask', ...
                'allOris', 'nSFs', 'nOris', 'rows', 'cols', 'masterX', 'masterY', 'encodingGamma');
            whichBlock = contains(meta.name, 'align');
            imSize   = ff.stimuli.imSize;
            encGamma = ff.encodingGamma;
            p1List   = meta.p1List;
            p2List   = meta.p2List;
            vMask    = meta.validMask;

            for tr = 1:nT
                if ~vMask(tr), continue; end
                im1Raw = double(imread(p1List{tr})) / 255;
                if ndims(im1Raw) == 3, im1Raw = rgb2gray(im1Raw); end
                im2Raw = double(imread(p2List{tr})) / 255;
                if ndims(im2Raw) == 3, im2Raw = rgb2gray(im2Raw); end

                im1 = min(max(imresize((im1Raw.^encGamma)*2 - 1, [imSize imSize]), -1), 1);
                im2 = min(max(imresize((im2Raw.^encGamma)*2 - 1, [imSize imSize]), -1), 1);
                targIm = ((im1/2 + .5).^(1/encGamma))*2 - 1;
                distIm = ((im2/2 + .5).^(1/encGamma))*2 - 1;

                [imOris, filtEn]   = steerFilterIm(targIm, ff.filt, ff.mask, ff.allOris, ff.nSFs, ff.nOris);
                [distOris, distEn] = steerFilterIm(distIm, ff.filt, ff.mask, ff.allOris, ff.nSFs, ff.nOris);

                gabIm = buildGaborIm(imOris, filtEn, ff.task, ff.stimuli, tr, whichBlock, ...
                    ff.rows, ff.cols, ff.masterX, ff.masterY);
                [gabOris, gabEn] = steerFilterIm(gabIm, ff.filt, ff.mask, ff.allOris, ff.nSFs, ff.nOris);

                [lV1(tr,:), rV1(tr,:), lEn(tr,:), rEn(tr,:), aMu(tr,:), aSD(tr,:), aMuU(tr,:), aSDU(tr,:)] = ...
                    computeV1Predictors(filtEn, distEn, gabEn, imOris, distOris, gabOris, ...
                    ff.task.targPos(tr), oriNoiseGrid, posBlurGrid, cfgList);
            end
        end

        sessLeftV1{s}   = lV1;  sessRightV1{s}  = rV1;
        sessLeftEn{s}   = lEn;  sessRightEn{s}  = rEn;
        sessAbsMu{s}    = aMu;  sessAbsSD{s}    = aSD;
        sessAbsMuUnw{s} = aMuU; sessAbsSDUnw{s} = aSDU;
    end

    leftV1Mat   = vertcat(sessLeftV1{:});
    rightV1Mat  = vertcat(sessRightV1{:});
    leftEnMat   = vertcat(sessLeftEn{:});
    rightEnMat  = vertcat(sessRightEn{:});
    absMuMat    = vertcat(sessAbsMu{:});
    absSDMat    = vertcat(sessAbsSD{:});
    absMuUnwMat = vertcat(sessAbsMuUnw{:});
    absSDUnwMat = vertcat(sessAbsSDUnw{:});

    fprintf('Feature extraction finished in %.1f mins. Saving %s...\n', toc(tStart)/60, cacheFile);
    save(cacheFile, 'leftV1Mat', 'rightV1Mat', 'leftEnMat', 'rightEnMat', ...
        'absMuMat', 'absSDMat', 'absMuUnwMat', 'absSDUnwMat', 'validTrial', ...
        'cfgList', 'oriNoiseGrid', 'posBlurGrid', '-v7.3');
end

%% 6. Filter to Complete Valid Trials & Fit Unified Models (5-Fold CV in parfor)
expIDVec   = expIDVec(validTrial);
partIDVec  = partIDVec(validTrial);
blockVec   = blockVec(validTrial);
condVec    = condVec(validTrial);
nGaborsVec = nGaborsVec(validTrial);
targIntVec = targIntVec(validTrial);
respVec    = respVec(validTrial);

leftV1Mat   = leftV1Mat(validTrial, :);
rightV1Mat  = rightV1Mat(validTrial, :);
leftEnMat   = leftEnMat(validTrial, :);
rightEnMat  = rightEnMat(validTrial, :);
absMuMat    = absMuMat(validTrial, :);
absSDMat    = absSDMat(validTrial, :);
absMuUnwMat = absMuUnwMat(validTrial, :);
absSDUnwMat = absSDUnwMat(validTrial, :);

totalTrials = numel(respVec);
fprintf('Fitting 6 unified models across %d valid trials (%d noise configs)...\n', totalTrials, nConfigs);

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
kFolds = 5;
cfgMetrics = zeros(nConfigs, nModels, 7); % [choiceLL, aic, bic, choiceAcc, idealAlign, idealRand, groupRMSE]

parfor cfg = 1:nConfigs
    rng(6042 + cfg);

    v1Diff     = zscore(rightV1Mat(:, cfg) - leftV1Mat(:, cfg));
    enDiff     = zscore(rightEnMat(:, cfg) - leftEnMat(:, cfg));
    absMuDiffZ = zscore(absMuMat(:, cfg));
    absSDDiffZ = zscore(absSDMat(:, cfg));
    absMuUnwZ  = zscore(absMuUnwMat(:, cfg));
    absSDUnwZ  = zscore(absSDUnwMat(:, cfg));

    T = table(expIDVec, partIDVec, blockVec, condVec, nGaborsVec, targIntVec, respVec, ...
        v1Diff, enDiff, absMuDiffZ, absSDDiffZ, absMuUnwZ, absSDUnwZ, ...
        'VariableNames', {'expID', 'partID', 'block', 'cond', 'nGabors', 'targInt', 'resp', ...
        'v1Diff', 'enDiff', 'absMuDiff', 'absSDDiff', 'absMuUnwDiff', 'absSDUnwDiff'});

    assert(~any(any(ismissing(T))), 'NaNs detected in config %d', cfg);

    permIdx = randperm(totalTrials);
    Tshuf = T(permIdx, :);
    localM = zeros(nModels, 7);

    for m = 1:nModels
        rhs = modelRHS{m};
        humanSpec = sprintf('resp ~ %s + (1 | partID)', rhs);
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
        nParams   = numel(strsplit(rhs, '+')) + 2;
        aic       = 2 * nParams - 2 * choiceLL;
        bic       = nParams * log(totalTrials) - 2 * choiceLL;
        choiceAcc = mean(round(predHumanProb) == Tshuf.resp);

        idealCorrect = double(round(predIdealProb) == Tshuf.targInt);
        idealAlign = mean(idealCorrect(Tshuf.block == 1));
        idealRand  = mean(idealCorrect(Tshuf.block == 0));

        humanCorrect = double(Tshuf.resp == Tshuf.targInt);
        modelCorrect = double(Tshuf.targInt == 1) .* predHumanProb + ...
                       double(Tshuf.targInt == 0) .* (1 - predHumanProb);
        indHum = qStats(humanCorrect, [Tshuf.partID, Tshuf.expID, Tshuf.block, Tshuf.cond, Tshuf.nGabors], 2);
        indMod = qStats(modelCorrect, [Tshuf.partID, Tshuf.expID, Tshuf.block, Tshuf.cond, Tshuf.nGabors], 2);
        grpHum = qStats(indHum(:, 6), indHum(:, 2:5), 2);
        grpMod = qStats(indMod(:, 6), indMod(:, 2:5), 2);
        groupRMSE = sqrt(mean((grpMod(:, 5) - grpHum(:, 5)).^2));

        localM(m, :) = [choiceLL, aic, bic, choiceAcc, idealAlign, idealRand, groupRMSE];
    end
    cfgMetrics(cfg, :, :) = localM;
end

%% 7. Save & Print Summary
summaryRows = table();
for o = 1:nOri
    for b = 1:nBlur
        mask = (cfgList(:, 1) == o) & (cfgList(:, 2) == b);
        meanM = squeeze(mean(cfgMetrics(mask, :, :), 1));
        for m = 1:nModels
            row = table(oriNoiseGrid(o), posBlurGrid(b), modelNames(m), ...
                meanM(m,1), meanM(m,2), meanM(m,3), meanM(m,4), meanM(m,5), meanM(m,6), meanM(m,7), ...
                'VariableNames', {'OriNoisePi', 'PosBlurPx', 'Model', ...
                'ChoiceLogLik', 'AIC', 'BIC', 'ChoiceAcc', 'IdealAcc_Align', 'IdealAcc_Rand', 'GroupRMSE'});
            summaryRows = [summaryRows; row]; %#ok<AGROW>
        end
    end
end

save('unifiedResults_v2_V1.mat', 'modelNames', 'modelRHS', 'oriNoiseGrid', 'posBlurGrid', 'cfgList', 'cfgMetrics', 'summaryRows', 'validTrial');
writetable(summaryRows, 'unifiedSummary_v2_V1.csv');

fprintf('\nTotal elapsed time: %.1f mins\n', toc(tStart)/60);
fprintf('\n=== Best Unified V1 Fit per Model (All 3 Experiments) ===\n');
for m = 1:nModels
    sub = summaryRows(strcmp(summaryRows.Model, modelNames{m}), :);
    [~, bestIdx] = max(sub.ChoiceLogLik);
    disp(sub(bestIdx, :));
end

%% model fits
plot_unified_allExps

%% Local Helper Functions
function key = extractFileKey(relPath)
    raw = strrep(char(relPath), '\', '/');
    [~, nm, ext] = fileparts(raw);
    key = lower([nm ext]);
end

function [oris, en] = steerFilterIm(im, filt, mask, allOris, nSFs, nOris)
    sz = size(im, 1);
    imFFT = fft2(im);
    filtMag = zeros(sz, sz, nOris, nSFs);
    oriSin  = zeros(sz, sz, nOris, nSFs);
    oriCos  = zeros(sz, sz, nOris, nSFs);
    oris    = zeros(sz, sz, nSFs);
    en      = zeros(sz, sz, nSFs);
    for sf = 1:nSFs
        for ori = 1:nOris
            thisIm = ifft2(imFFT .* filt(:,:,ori,sf), 'symmetric');
            filtMag(:,:,ori,sf) = fftshift(thisIm);
            thisOri = allOris(ori);
            oriSin(:,:,ori,sf) = filtMag(:,:,ori,sf) .* sin(deg2rad(thisOri));
            oriCos(:,:,ori,sf) = filtMag(:,:,ori,sf) .* cos(deg2rad(thisOri));
        end
        oris(:,:,sf) = atan2(sum(oriSin(:,:,:,sf), 3), sum(oriCos(:,:,:,sf), 3));
        en(:,:,sf)   = abs(complex(filtMag(:,:,1,sf), filtMag(:,:,2,sf))) .* mask;
    end
end

function gabIm = buildGaborIm(imOris, filtEn, task, stimuli, trial, whichBlock, rows, cols, masterX, masterY)
    imSize = stimuli.imSize;
    outIm = zeros(imSize * 2);
    for sf = 1:size(imOris, 3)
        nPlaces = task.nGabors(trial);
        thisOris = imOris(:,:,sf) + pi;
        sd = stimuli.filtSize(sf);
        thisSFMag = filtEn(:,:,sf);
        x = masterX; y = masterY;
        for loc = 1:nPlaces
            [~, thisPlace] = max(thisSFMag(:));
            gabCoords = [rows(thisPlace), cols(thisPlace)];
            thisSFMag = thisSFMag .* (1 - Gaussian2D(sd * stimuli.gaborSeparation, gabCoords, imSize));
            if ~whichBlock
                [xN, yN] = pol2cart(rand*2*pi, rand*imSize/2);
                gabCoords = [xN, yN] + imSize/2;
            end
            gabCoords = min(max(gabCoords, 1), imSize*2);
            th = thisOris(thisPlace);
            x = (x - gabCoords(1)) / sd;
            y = (y - gabCoords(2)) / sd;
            outIm = outIm + cos(th)*(-2*x.*exp(-(x.^2+y.^2))) + sin(th)*(-2*y.*exp(-(x.^2+y.^2)));
            x = x*sd + gabCoords(1);
            y = y*sd + gabCoords(2);
        end
    end
    gabIm = maxContrast(outIm);
    gabIm = gabIm(1:imSize, 1:imSize);
end

function [rLV1, rRV1, rLEn, rREn, rAM, rAS, rAMu, rASu] = computeV1Predictors( ...
    filtEn0, distEn0, gabEn0, imOris0, distOris0, gabOris0, tPos, oriNoiseGrid, posBlurGrid, cfgList)

    nBlur = numel(posBlurGrid);
    nConfigs = size(cfgList, 1);

    % Map orientations to doubled-angle [0, 2*pi) so cos(o1 - o2) = cos(2*(theta1 - theta2))
    imOris0(imOris0 < 0)     = imOris0(imOris0 < 0)     + pi; imOris0   = wrap(imOris0 * 2);
    distOris0(distOris0 < 0) = distOris0(distOris0 < 0) + pi; distOris0 = wrap(distOris0 * 2);
    gabOris0 = gabOris0 - pi;
    gabOris0(gabOris0 < 0)   = gabOris0(gabOris0 < 0)   + pi; gabOris0  = wrap(gabOris0 * 2);

    validGab   = ~isnan(gabOris0(:));
    gabOris0V  = gabOris0(validGab);
    imOris0V   = imOris0(validGab);
    distOris0V = distOris0(validGab);

    filtEnBase = max(filtEn0, 0);
    distEnBase = max(distEn0, 0);
    gabEnBase  = max(gabEn0,  0);

    filtEnB  = cell(nBlur, 1); distEnB  = cell(nBlur, 1);
    filtEnBV = cell(nBlur, 1); distEnBV = cell(nBlur, 1);
    gabEnBV  = cell(nBlur, 1); gabNormB = zeros(nBlur, 1);
    targEnB  = zeros(nBlur, 1); distEnMatchB = zeros(nBlur, 1);

    for bIdx = 1:nBlur
        sig = posBlurGrid(bIdx);
        if sig > 0
            fE = max(imgaussfilt(filtEnBase, sig), 0);
            dE = max(imgaussfilt(distEnBase, sig), 0);
            gE = max(imgaussfilt(gabEnBase,  sig), 0);
        else
            fE = filtEnBase; dE = distEnBase; gE = gabEnBase;
        end
        filtEnB{bIdx}  = fE(:);
        distEnB{bIdx}  = dE(:);
        filtEnBV{bIdx} = fE(validGab);
        distEnBV{bIdx} = dE(validGab);
        gEV            = gE(validGab);
        gabEnBV{bIdx}  = gEV;
        gNorm          = norm(gEV) + eps;
        gabNormB(bIdx) = gNorm;

        % Unoriented spatial energy matched filter (SDT sqrt(N) SNR scaling)
        targEnB(bIdx)      = sum(filtEnBV{bIdx} .* gEV) / gNorm;
        distEnMatchB(bIdx) = sum(distEnBV{bIdx} .* gEV) / gNorm;
    end

    rLV1 = zeros(1, nConfigs); rRV1 = zeros(1, nConfigs);
    rLEn = zeros(1, nConfigs); rREn = zeros(1, nConfigs);
    rAM  = zeros(1, nConfigs); rAS  = zeros(1, nConfigs);
    rAMu = zeros(1, nConfigs); rASu = zeros(1, nConfigs);

    for cfg = 1:nConfigs
        oCfg = cfgList(cfg, 1);
        bCfg = cfgList(cfg, 2);
        oNoise = oriNoiseGrid(oCfg) * pi;

        if tPos == 1
            rLEn(cfg) = targEnB(bCfg); rREn(cfg) = distEnMatchB(bCfg);
        else
            rLEn(cfg) = distEnMatchB(bCfg); rREn(cfg) = targEnB(bCfg);
        end

        if oNoise > 0
            imO    = wrap(imOris0   + randn(size(imOris0))   * oNoise);
            distO  = wrap(distOris0 + randn(size(distOris0)) * oNoise);
            gabOV  = wrap(gabOris0V + randn(size(gabOris0V)) * oNoise);
            imOV   = imO(validGab);
            distOV = distO(validGab);
        else
            imO = imOris0; distO = distOris0; gabOV = gabOris0V;
            imOV = imOris0V; distOV = distOris0V;
        end

        % V1 Oriented Matched-Filter SNR: sum(E_I * E_G * cos(2*dTheta)) / ||E_G||
        gEV   = gabEnBV{bCfg};
        gNorm = gabNormB(bCfg);
        tV1 = sum(filtEnBV{bCfg} .* gEV .* cos(imOV   - gabOV)) / gNorm;
        dV1 = sum(distEnBV{bCfg} .* gEV .* cos(distOV - gabOV)) / gNorm;
        if tPos == 1
            rLV1(cfg) = tV1; rRV1(cfg) = dV1;
        else
            rLV1(cfg) = dV1; rRV1(cfg) = tV1;
        end

        % Energy-weighted ensemble summary statistics
        fEV = filtEnB{bCfg}; dEV = distEnB{bCfg};
        [~, gSD] = circ_std(gabOV, gEV);
        [~, tSD] = circ_std(imO(:), fEV);
        [~, dSD] = circ_std(distO(:), dEV);
        gMu = circ_mean(gabOV, gEV);
        tMu = circ_mean(imO(:), fEV);
        dMu = circ_mean(distO(:), dEV);

        tMuErr = abs(circ_dist(gMu, tMu)); dMuErr = abs(circ_dist(gMu, dMu));
        tSDErr = abs(gSD - tSD);           dSDErr = abs(gSD - dSD);
        if tPos == 1
            rAM(cfg) = tMuErr - dMuErr; rAS(cfg) = tSDErr - dSDErr;
        else
            rAM(cfg) = dMuErr - tMuErr; rAS(cfg) = dSDErr - tSDErr;
        end

        % Unweighted ensemble summary statistics
        [~, gSDu] = circ_std(gabOV);
        [~, tSDu] = circ_std(imO(:));
        [~, dSDu] = circ_std(distO(:));
        gMuu = circ_mean(gabOV); tMuu = circ_mean(imO(:)); dMuu = circ_mean(distO(:));
        tMuEu = abs(circ_dist(gMuu, tMuu)); dMuEu = abs(circ_dist(gMuu, dMuu));
        tSDEu = abs(gSDu - tSDu);           dSDEu = abs(gSDu - dSDu);
        if tPos == 1
            rAMu(cfg) = tMuEu - dMuEu; rASu(cfg) = tSDEu - dSDEu;
        else
            rAMu(cfg) = dMuEu - tMuEu; rASu(cfg) = dSDEu - tSDEu;
        end
    end
end
