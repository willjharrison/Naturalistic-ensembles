clear
close all
clc

% f = dir('noiseModelOutput_fine/*.mat');



%% set up variables
noiseID = [];
repID = [];
foldID = [];
partID = [];
block = [];
sameDiff = [];
nGabors = [];
targInt = [];
modProb = [];

datNames = {'noiseID', 'repID', 'foldID', 'partID', 'block', 'sameDiff', ...
    'nGabors', 'targInt', 'modProb'};
fittedData = table(noiseID, repID, foldID, partID, block, sameDiff, ...
    nGabors, targInt, modProb);

%% load behavioural data
load('meanAccuracy.mat');
Y = (meanAccuracy.meanAccuracy);

%% load model sims
tic;
load('allFits_cells.mat');
timeToOpen=toc

nFiles = numel(allData);
datSize = size(allData);
%% run it
counter = 1;
k = 5;
tic;
parfor fileLoop = 1:nFiles
    
    %% open data
%     loadedDat = load(fullfile(f(fileLoop).folder, f(fileLoop).name));
    loadedDat = allData{fileLoop};
    
    % grab noise id and repetition id
%     noiseID = str2double(f(fileLoop).name(19:21));
%     repID = str2double(f(fileLoop).name(31:32));
    
    [noiseID, repID] = ind2sub(datSize,fileLoop);
    
    loadedDat.leftRightMatchDiff = loadedDat.leftMatchTrans - loadedDat.rightMatchTrans;
    %% use k-fold cross validation
    
    
    
    % shuffle data before fitting so we remove sequential effects from
    % different folds
    shuffledData = loadedDat(randperm(size(loadedDat,1)),:);
    
    modSpec = ['targInt ~ leftRightOriDiff + leftRightMatchDiff'];
    
    for foldLoop = 1:k
        
        testData = shuffledData((1:size(loadedDat,1))<=size(loadedDat,1)/k*foldLoop & (1:size(loadedDat,1)) > size(loadedDat,1)/k*(foldLoop -1),:);
        trainData = shuffledData(~((1:size(loadedDat,1))<=size(loadedDat,1)/k*foldLoop) | ~((1:size(loadedDat,1)) > size(loadedDat,1)/k*(foldLoop -1)),:);
        
        dataModel = fitglme(trainData, modSpec,'Distribution', 'Binomial');
        
        
        newData = splitvars(table([repmat(noiseID,size(testData,1),1) ...
            repmat(repID, size(testData,1),1) ...
            repmat(foldLoop, size(testData,1),1) ...
            table2array(testData(:,1:5)) predict(dataModel, testData)]));
        newData.Properties.VariableNames = datNames;
        
        fittedData = vertcat(fittedData, newData);
        
        % summarise data and fit to behaviour
        indSum = qStats(round(newData.modProb)==newData.targInt, table2array(newData(:,[4 5 6 7])));
        groupSum = qStats(indSum(:,5), indSum(:,2:4));
        
        X = groupSum(:,4);
        
        
        thisFit = fitglm(X, Y);
        gf(foldLoop,:,fileLoop) = [noiseID repID thisFit.LogLikelihood thisFit.Coefficients.Estimate' thisFit.Rsquared.Adjusted];
        
        modFitData(:,:,foldLoop,fileLoop) = [repmat(noiseID,length(groupSum),1) ...
            repmat(repID,length(groupSum),1) ...
            groupSum(:,1:3) thisFit.predict];
        
        
    end
    
%     if exist('l')
%         fprintf(repmat('\b',1,l))
%     end
%     l = fprintf('File: %s out of %s\rElapsed: %.1f mins\rProgress: %2.1f%%\r', string(fileLoop), string(nFiles), toc/60,(fileLoop)/nFiles*100);
%     pause(.1);
    
end
totalTime = toc/60;
fprintf('Total processing time was: %.1f mins\n', totalTime)

%% reshape to make analysis easier
gfRes = reshape(permute(gf,[3 1 2]), numel(gf)/size(gf,2), size(gf,2));
modSummary = qStats(gfRes(:,3:6), gfRes(:,1),2);
[~,idx] = max(modSummary(:,2));
whichNoise = modSummary(idx,1);


modFitRes = reshape(permute(modFitData,[4 3 1 2]), numel(modFitData)/6, 6);
noiseIdx = modFitRes(:,1) == whichNoise;

a = qStats(modFitRes(noiseIdx,6), modFitRes(noiseIdx,[1 3:5]))
% a(:,1) = round(a(:,1))

%%
modResults = qStats(round(fittedData.modProb)==fittedData.targInt, table2array(fittedData(:,[1 5:7])));
save('modResults.mat','modResults','fittedData','modFitData','gf')