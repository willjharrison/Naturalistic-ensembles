%% compare_models.m
% Compares candidate observer model variations using 5-fold cross-validation
% against existing simulation data (../Matlab/Analysis/SCmodel/allFits_cells.mat)
% and empirical group accuracy (../Matlab/Analysis/SCmodel/meanAccuracy.mat).
%
% Outputs are saved inside Working/ to keep existing folders untouched.

clear; close all; clc;

scModelDir = '../Matlab/Analysis/SCmodel';
load(fullfile(scModelDir, 'meanAccuracy.mat'), 'meanAccuracy');
Y = meanAccuracy.meanAccuracy;

fprintf('Loading allFits_cells.mat...\n');
tic;
load(fullfile(scModelDir, 'allFits_cells.mat'), 'allData');
fprintf('Loaded in %.1f s\n', toc);

nFiles = numel(allData);
datSize = size(allData);
k = 5;

%% Define model variations
modelNames = { ...
    'SpatialOnly_Match', ...
    'Spatial_Orientation', ...
    'Spatial_Match_Plus_Ori', ...
    'SummaryStats_Weighted', ...
    'SummaryStats_Unweighted', ...
    'Full_Combined'};

modelSpecs = { ...
    'targInt ~ leftRightMatchDiff', ...
    'targInt ~ leftRightOriDiff', ...
    'targInt ~ leftRightOriDiff + leftRightMatchDiff', ...
    'targInt ~ leftRightMuDiff + leftRightSDDiff', ...
    'targInt ~ leftRightMuUnwDiff + leftRightSDUnwDiff', ...
    'targInt ~ leftRightOriDiff + leftRightMatchDiff + leftRightMuDiff + leftRightSDDiff'};

nModels = numel(modelNames);
summaryRows = table();
allGF = cell(nModels, 1);

for m = 1:nModels
    modName = modelNames{m};
    modSpec = modelSpecs{m};
    fprintf('\n=== Running Model %d/%d: %s ===\n', m, nModels, modName);

    gf = zeros(k, 5, nFiles); % [noiseID, repID, logLik, adjR2, trialAcc]

    parfor fileLoop = 1:nFiles
        rng(1042 + fileLoop);
        loadedDat = allData{fileLoop};
        [noiseID, repID] = ind2sub(datSize, fileLoop);

        % Derived difference predictors (positive = favours right option, or lower error for right)
        loadedDat.leftRightMatchDiff = loadedDat.leftMatchTrans - loadedDat.rightMatchTrans;
        loadedDat.leftRightMuDiff    = abs(loadedDat.leftMu) - abs(loadedDat.rightMu);
        loadedDat.leftRightSDDiff    = abs(loadedDat.leftSDDiff) - abs(loadedDat.rightSDDiff);
        loadedDat.leftRightMuUnwDiff = abs(loadedDat.leftMuUnweighted) - abs(loadedDat.rightMuUnweighted);
        loadedDat.leftRightSDUnwDiff = abs(loadedDat.leftSDUnweighted) - abs(loadedDat.rightSDUnweighted);

        nTrials = height(loadedDat);
        shuffledData = loadedDat(randperm(nTrials), :);
        localGF = zeros(k, 5);

        for foldLoop = 1:k
            testMask = ((1:nTrials) <= nTrials/k * foldLoop) & ((1:nTrials) > nTrials/k * (foldLoop - 1));
            testData = shuffledData(testMask, :);
            trainData = shuffledData(~testMask, :);

            dataModel = fitglme(trainData, modSpec, 'Distribution', 'Binomial');
            modProb = predict(dataModel, testData);

            trialCorrect = double(round(modProb) == testData.targInt);
            indSum = qStats(trialCorrect, [testData.participant, testData.block, testData.sameDiff, testData.nGabors]);
            groupSum = qStats(indSum(:, 5), indSum(:, 2:4));

            X = groupSum(:, 4);
            thisFit = fitglm(X, Y);

            localGF(foldLoop, :) = [noiseID, repID, thisFit.LogLikelihood, thisFit.Rsquared.Adjusted, mean(trialCorrect)];
        end

        gf(:, :, fileLoop) = localGF;
    end

    allGF{m} = gf;

    % Aggregate across folds & reps per noise level
    gfRes = reshape(permute(gf, [3 1 2]), numel(gf)/5, 5);
    noiseSummary = qStats(gfRes(:, 3:5), gfRes(:, 1), 2); % [noiseID, meanLogLik, meanAdjR2, meanTrialAcc]
    [bestLL, bestIdx] = max(noiseSummary(:, 2));

    row = table({modName}, {modSpec}, noiseSummary(bestIdx, 1), bestLL, noiseSummary(bestIdx, 3), noiseSummary(bestIdx, 4), ...
        'VariableNames', {'Model', 'Formula', 'BestNoiseID', 'BestLogLik', 'AdjR2_atBest', 'TrialAcc_atBest'});
    summaryRows = [summaryRows; row];
    disp(row);
end

%% Save & display summary
save('modelComparisonResults.mat', 'modelNames', 'modelSpecs', 'allGF', 'summaryRows');
writetable(summaryRows, 'modelComparisonSummary.csv');
disp(summaryRows);
