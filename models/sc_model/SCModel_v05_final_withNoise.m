clear
close all
clc

%% find all images
encodingGamma = 2;
load('../../../../SteeredFilterPlacement/THINGS/Variables/image_concept_index.mat')
load('../../../../SteeredFilterPlacement/THINGS/Variables/image_paths.mat')
load('../../../../SteeredFilterPlacement/THINGS/Variables/unique_id.mat')

%% find files

files = dir('../../Output/');
files = files(~ismember({files.name}, {'.','..','.DS_Store'}));
nFiles = length(files);

simScore = @(imA, imB) dot(imA(:), imB(:))/(norm(imA(:))*norm(imB(:)));

participant = [];
block = [];
nGabors = [];
sameDiff = [];
targInt = [];
resp = [];
leftMatch = zeros(12000,1);
rightMatch = leftMatch;
leftOri = leftMatch;
rightOri = leftMatch;
leftOriMean = leftMatch;
rightOriMean = rightMatch;
leftMu = leftMatch;
rightMu = leftMatch;
leftSD = leftMatch;
rightSD = leftMatch;
leftMuUnweighted = leftMatch;
rightMuUnweighted = leftMatch;
leftSDUnweighted = leftMatch;
rightSDUnweighted = leftMatch;

targCat = [];
distCat = [];
offset = .001; % allows to compute log of energy, which otherwise includes 0s.
tic;
for fileLoop = 1:nFiles
    
    load([files(fileLoop).folder '/' files(fileLoop).name],'task','dv','subject_name','stimuli')
    thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';
    
    if fileLoop == 1
        [~, radDist] = polarDistFun(stimuli.imSize);
        
        mask = ...
            createLowPassFilter(radDist,stimuli.imSize/2 - 6, 6);
    end
    
    
    % get the block data - is this aligned or randomised?
    whichBlock = contains(files(fileLoop).name,'align');
    
    % store experiment variables
    participant = [participant; repmat(str2double(subject_name),600,1)];
    block = [block; repmat(whichBlock,600,1)];
    nGabors = [nGabors; task.nGabors];
    sameDiff = [sameDiff; task.sameDiff];
    targInt = [targInt; task.targPos - 1];
    targCat = [targCat; stimuli.allCats(:,1)];
    distCat = [distCat; stimuli.allCats(:,2)];
    resp = [resp; dv.chosenSquare-1];
    
    
    %% start trial loop
    for trial = 1:task.nTrials
        
        %% load image files
        stimFileName = sprintf(['modIms/gabIm_trial_%03d_fileLoop_%02d.mat'], trial, fileLoop);
        load(stimFileName, 'gabIm', 'targIm', 'distIm','gabEn','filtEn','distEn','imOris','distOris','theta','gabOris');
        %
        % first make measures phase-invariant
        imOris = imOris;
        imOris(imOris<0) = imOris(imOris<0) + pi;
        imOris = wrap(imOris*2);
        
        distOris = distOris;
        distOris(distOris<0) = distOris(distOris<0) + pi;
        distOris = wrap(distOris*2);
        
        gabOris = gabOris - pi; % correct for image generation (see above)
        gabOris(gabOris<0) = gabOris(gabOris<0) + pi;
        gabOris = wrap(gabOris*2);
        %
        %% gab En similarity
        
        if task.targPos(trial) == 1
            leftMatch(trial + 600*(fileLoop-1),1) = dot(filtEn(:), gabEn(:));% max(xcorr2(filtEn, gabEn),[],[1 2]);
            rightMatch(trial + 600*(fileLoop-1),1) = dot(distEn(:), gabEn(:)); %max(xcorr2(distEn, gabEn),[],[1 2]);
            
            
        else
            
            leftMatch(trial + 600*(fileLoop-1),1) = dot(distEn(:), gabEn(:)); %max(xcorr2(distEn, gabEn),[],[1 2]);
            rightMatch(trial + 600*(fileLoop-1),1) = dot(filtEn(:), gabEn(:)); % max(xcorr2(filtEn, gabEn),[],[1 2]);
        end
        
        
        %% orientation similarity
        if task.targPos(trial) == 1
            
            leftOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(imOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            rightOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(distOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            
            leftOriMean(trial + 600*(fileLoop-1),1) = circ_mean(circ_dist(imOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))));
            rightOriMean(trial + 600*(fileLoop-1),1) = circ_mean(circ_dist(distOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))));
            
        else
            
            leftOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(distOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            rightOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(imOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            
            leftOriMean(trial + 600*(fileLoop-1),1) = circ_mean(circ_dist(distOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))));
            rightOriMean(trial + 600*(fileLoop-1),1) = circ_mean(circ_dist(imOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))));
            
        end
        
        %% summary stat similarity weighted by energy
        
        [~,gabSD(trial + 600*(fileLoop-1),1)] = circ_std(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:))));
        [~,targSD] = circ_std(imOris(:), filtEn(:));
        [~,distSD] = circ_std(distOris(:), distEn(:));
        
        if task.targPos(trial) == 1
            
            leftMu(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:)))),circ_mean(imOris(:), filtEn(:)));
            rightMu(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:)))),circ_mean(distOris(:), distEn(:)));
            
            leftSD(trial + 600*(fileLoop-1),1) = targSD;
            rightSD(trial + 600*(fileLoop-1),1) = distSD;
            
        else
            
            leftMu(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:)))),circ_mean(distOris(:), distEn(:)));
            rightMu(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:)))),circ_mean(imOris(:), filtEn(:)));
            
            leftSD(trial + 600*(fileLoop-1),1) = distSD;
            rightSD(trial + 600*(fileLoop-1),1) = targSD;
            
        end
        
        
        %% summary stat similarity unweighted
        
        if task.targPos(trial) == 1
            
            leftMuUnweighted(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:)))),circ_mean(imOris(:)));
            rightMuUnweighted(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:)))),circ_mean(distOris(:)));
            
            leftSDUnweighted(trial + 600*(fileLoop-1),1) = (circ_std(gabOris(~isnan(gabOris(:))))-circ_std(imOris(:)));
            rightSDUnweighted(trial + 600*(fileLoop-1),1) = (circ_std(gabOris(~isnan(gabOris(:))))-circ_std(distOris(:)));
            
        else
            
            leftMuUnweighted(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:)))),circ_mean(distOris(:)));
            rightMuUnweighted(trial + 600*(fileLoop-1),1) = circ_dist(circ_mean(gabOris(~isnan(gabOris(:)))),circ_mean(imOris(:)));
            
            leftSDUnweighted(trial + 600*(fileLoop-1),1) = (circ_std(gabOris(~isnan(gabOris(:))))-circ_std(distOris(:)));
            rightSDUnweighted(trial + 600*(fileLoop-1),1) = (circ_std(gabOris(~isnan(gabOris(:))))-circ_std(imOris(:)));
            
        end
        
        if mod(trial,300) == 0
            
            if exist('l')
                fprintf(repmat('\b',1,l))
            end
            l = fprintf('File: %s \rTrials: %s\rElapsed: %.1f mins\rProgress: %2.1f%%\r', files(fileLoop).name(17:24), string(trial), toc/60,(fileLoop-1)/nFiles*100);
            pause(.1);
        end
        
        
        %%
    end
    
end % end fileLoop

%% summarise predictors
clc
nNoiseLevels = 10;
noiseReps = 20;
oriNoise = linspace(0,1,nNoiseLevels)*4*std([leftOri; rightOri]);
matchNoise = linspace(0,1,nNoiseLevels)*4*std([leftMatch; rightMatch]);
k = 5;
rng(1042); % for reproducibility
fittedData = zeros(12000*nNoiseLevels,8);
obsAv = [];
for repLoop = 1%:noiseReps
    for noiseLoop = 1:nNoiseLevels
        
        % add late gaussian noise
        oriNoiseLevel = oriNoise(noiseLoop);
        matchNoiseLevel = matchNoise(noiseLoop);
        leftOriNoise = leftOri + randn(size(leftOri))*oriNoiseLevel;
        rightOriNoise = rightOri + randn(size(rightOri))*oriNoiseLevel;
        leftRightOriDiff = leftOri - rightOri + randn(size(leftOri))*oriNoiseLevel;
        leftMatchNoise = leftMatch + randn(size(leftMatch))*matchNoiseLevel;
        rightMatchNoise = rightMatch + randn(size(rightMatch))*matchNoiseLevel;
        allData = table(participant, block, sameDiff, nGabors, targInt, ...
            leftOriNoise, rightOriNoise, leftRightOriDiff, leftMatch, rightMatch, leftMatchNoise, rightMatchNoise);
        
        %% use k-fold cross validation
        
        % shuffle data before fitting so we remove sequential effects from
        % different folds
        shuffledData = allData(randperm(size(allData,1)),:);
        
        %     modSpec = ['targInt ~ leftOriNoise*rightOriNoise+leftMatchNoise+rightMatchNoise'];
        modSpec = ['targInt ~ leftOriNoise+rightOriNoise+leftMatchNoise+rightMatchNoise'];
        
        for foldLoop = 1:k
            
            testData = shuffledData((1:size(allData,1))<=size(allData,1)/k*foldLoop & (1:size(allData,1)) > size(allData,1)/k*(foldLoop -1),:);
            trainData = shuffledData(~((1:size(allData,1))<=size(allData,1)/k*foldLoop) | ~((1:size(allData,1)) > size(allData,1)/k*(foldLoop -1)),:);
            
            dataModel = fitglme(trainData, modSpec,'Distribution', 'Binomial');
            %         fittedData = [fittedData; [repmat(foldLoop, size(testData,1),1) repmat(noiseLoop, size(testData,1),1) table2array(testData(:,1:5)) predict(dataModel, testData)]];
            fittedData((1:size(allData,1)/k)+size(allData,1)/k*(foldLoop-1) + 12000*(noiseLoop-1),:) = [repmat(foldLoop, size(testData,1),1) repmat(noiseLoop, size(testData,1),1) table2array(testData(:,1:5)) predict(dataModel, testData)];
        end
        %     dataModel
        %     size(fittedData)
        
    end
    
    % analyse k-fold data
    
    % aggregate across folds
    foldAv = qStats(round(fittedData(:,8))==fittedData(:,7), fittedData(:,2:6),2);
    
    % average over observers
    obsAv = [obsAv; qStats(foldAv(:,6), foldAv(:,[1 3:5]))];
    
end

obsAv = qStats(obsAv(:,5), obsAv(:,1:4));

%% observer accuracy for Tom
% a = qStats(allData.resp == allData.targInt, [allData.participant, allData.block, allData.sameDiff, allData.nGabors],2);
% a = qStats(a(:,5), a(:,2:4));
%
% block = a(:,1);
% sameDiff = a(:,2);
% nGabors = a(:,3);
% meanAccuracy = a(:,4);
% sem = a(:,5);
%
% meanAccuracy = table(block,sameDiff,nGabors,meanAccuracy,sem);
% writetable(meanAccuracy, 'meanAccuracy.csv')


%%
totalTime = toc
save('allFits_v05.mat')