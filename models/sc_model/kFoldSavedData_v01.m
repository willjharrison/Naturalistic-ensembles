clear
close all
clc

f = dir('noiseModelOutput/*.mat');

nFiles = length(f);

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

load('meanAccuracy.mat');

%% run it
counter = 1;
tic;
for fileLoop = 1:nFiles
    
    %% open data
    load(fullfile(f(fileLoop).folder, f(fileLoop).name))
    
    % grab noise id and repetition id
    noiseID = str2double(f(fileLoop).name(19:21));
    repID = str2double(f(fileLoop).name(31:32));
    
    allData.leftRightMatchDiff = allData.leftMatchTrans - allData.rightMatchTrans;
    %% use k-fold cross validation
    
    k = 5;
    
    % shuffle data before fitting so we remove sequential effects from
    % different folds
    shuffledData = allData(randperm(size(allData,1)),:);
    
    modSpec = ['targInt ~ leftRightOriDiff + leftRightMatchDiff'];
    
    for foldLoop = 1:k
        
        testData = shuffledData((1:size(allData,1))<=size(allData,1)/k*foldLoop & (1:size(allData,1)) > size(allData,1)/k*(foldLoop -1),:);
        trainData = shuffledData(~((1:size(allData,1))<=size(allData,1)/k*foldLoop) | ~((1:size(allData,1)) > size(allData,1)/k*(foldLoop -1)),:);
        
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
        
        if counter == 1
            Y = (meanAccuracy.meanAccuracy);
        end
        modFit.f{counter} = fitglm(X, Y);
        modFit.gf(counter,:) = [noiseID repID modFit.f{counter}.LogLikelihood];
        if counter == 1
            modFit.allData = [repmat(noiseID,length(groupSum),1) ...
                repmat(repID,length(groupSum),1) groupSum(:,1:3) modFit.f{counter}.predict];
        else
            modFit.allData = [modFit.allData; repmat(noiseID,length(groupSum),1) ...
                repmat(repID,length(groupSum),1) ...
                groupSum(:,1:3) modFit.f{counter}.predict];
        end
        counter = counter + 1;
    end
    
    if exist('l')
        fprintf(repmat('\b',1,l))
    end
    l = fprintf('File: %s out of %s\rElapsed: %.1f mins\rProgress: %2.1f%%\r', string(fileLoop), string(nFiles), toc/60,(fileLoop)/nFiles*100);
    pause(.1);
    
end
% if exist('l')
%         fprintf(repmat('\b',1,l))
%     end
%     l = fprintf('File: %s \rElapsed: %.1f mins\rProgress: %2.1f%%\r', string(fileLoop), toc/60,(fileLoop)/nFiles*100);
%%
modResults = qStats(round(fittedData.modProb)==fittedData.targInt, table2array(fittedData(:,[1 5:7])));
save('modResults.mat','modResults','fittedData','modFit')