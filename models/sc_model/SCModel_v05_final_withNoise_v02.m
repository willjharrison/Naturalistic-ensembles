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

%% observer accuracy for Tom
a = qStats(allData.resp == allData.targInt, [allData.participant, allData.block, allData.sameDiff, allData.nGabors],2);
a = qStats(a(:,5), a(:,2:4));

block = a(:,1);
sameDiff = a(:,2);
nGabors = a(:,3);
meanAccuracy = a(:,4);
sem = a(:,5);

meanAccuracy = table(block,sameDiff,nGabors,meanAccuracy,sem);
% save('meanAccuracy.mat','meanAccuracy')
% writetable(meanAccuracy, 'meanAccuracy.csv')

%% summarise predictors
% clc


% add late gaussian noise
% oriNoiseLevel = oriNoise(noiseLoop);
% matchNoiseLevel = matchNoise(noiseLoop);
% leftOriNoise = leftOri + randn(size(leftOri))*oriNoiseLevel;
% rightOriNoise = rightOri + randn(size(rightOri))*oriNoiseLevel;
leftRightOriDiff = leftOri - rightOri; % + randn(size(leftOri))*oriNoiseLevel;
% leftMatchNoise = leftMatch + randn(size(leftMatch))*matchNoiseLevel;
% rightMatchNoise = rightMatch + randn(size(rightMatch))*matchNoiseLevel;
allData = table(participant, block, sameDiff, nGabors, targInt, resp, ...
    leftOri, rightOri, leftRightOriDiff, leftMatch, rightMatch); %, leftMatchNoise, rightMatchNoise);

% fit noise
noiseLoops = 10;
k = 5;
rng(1042); % for reproducibility
costFun = @(x) simData(allData, k, noiseLoops, x);

oriNoiseLevels = [0 logspace(1,1.5,5)]';
% [oriNoiseLevels, matchNoiseLevels] = BalanceFactors(1,0,oriNoiseLevels,oriNoiseLevels);
[oriNoiseLevels, matchNoiseLevels] = BalanceFactors(1,0,oriNoiseLevels,oriNoiseLevels);

nNoiseLevels = length(oriNoiseLevels);
negLL = [];%zeros(length(oriNoiseLevels),1);
tic;
parfor noiseLoop = 1:nNoiseLevels
       
    negLL = [negLL; costFun([oriNoiseLevels(noiseLoop,1), matchNoiseLevels(noiseLoop,1)])];
%     negLL = [negLL; costFun([oriNoiseLevels(noiseLoop,1)])];
    
end
timeToRunParFor = toc
figure;
if length(unique(matchNoiseLevels))>1
    subplot(1,2,1)
    plot(oriNoiseLevels,negLL)
    hold
    subplot(1,2,2)
    plot(matchNoiseLevels,negLL,'*')
else
    plot(oriNoiseLevels,negLL)
end
[~,minLL] = min(negLL);

noiseParams = [oriNoiseLevels(minLL) matchNoiseLevels(minLL)]
% noiseParams = [oriNoiseLevels(minLL)]

save('modelPredicts.mat')

% obsAv(:,4) = 1./(1 + exp(-(g.Coefficients.Estimate(1) + g.Coefficients.Estimate(2).*obsAv(:,4))))

%% summarise model data
[~, g, fittedData, dataModel] = simData(allData, k, 1, noiseParams);
% yHat = g.Coefficients.Estimate(1) + g.Coefficients.Estimate(2).*fittedData(:,7);

% yHat = 1./(1 + exp(-(g.Coefficients.Estimate(1) + g.Coefficients.Estimate(2).*fittedData(:,7))));

foldAv = qStats(round(g.predict)==fittedData(:,6), fittedData(:,2:5),2);

% average over observers
obsAv = [qStats(foldAv(:,5), foldAv(:,[2:4]))]


%%

oriNoiseLevels = [0 logspace(0,2,20)]';
% [oriNoiseLevels, matchNoiseLevels] = BalanceFactors(1,0,oriNoiseLevels,oriNoiseLevels);
[oriNoiseLevels, matchNoiseLevels] = BalanceFactors(1,0,oriNoiseLevels,oriNoiseLevels);
nNoiseLevels = length(oriNoiseLevels);
fd = [];
tic;
for noiseLoop = 1:nNoiseLevels

    thisPs = [oriNoiseLevels(noiseLoop,1), matchNoiseLevels(noiseLoop,1)];
    [n, g, fittedData, dataModel] = simData(allData, k, 1, thisPs);

%     fd = [fd; repmat(noiseLoop, size(fittedData,1),1) fittedData];

    % average over folds/observers
    foldAv = qStats(round(g.predict)==fittedData(:,6), fittedData(:,2:5),2);
    fd = [fd; repmat(noiseLoop, 20,1) qStats(foldAv(:,5), foldAv(:,[2:4]))];

end

% 


% 
% %% use k-fold cross validation
% 
% % shuffle data before fitting so we remove sequential effects from
% % different folds
% shuffledData = allData(randperm(size(allData,1)),:);
% 
% %     modSpec = ['targInt ~ leftOriNoise*rightOriNoise+leftMatchNoise+rightMatchNoise'];
% modSpec = ['targInt ~ leftRightOriDiff'];
% fittedData = zeros(12000,7);
% for foldLoop = 1:k
%     
%     testData = shuffledData((1:size(allData,1))<=size(allData,1)/k*foldLoop & (1:size(allData,1)) > size(allData,1)/k*(foldLoop -1),:);
%     trainData = shuffledData(~((1:size(allData,1))<=size(allData,1)/k*foldLoop) | ~((1:size(allData,1)) > size(allData,1)/k*(foldLoop -1)),:);
%     
%     dataModel = fitglme(trainData, modSpec,'Distribution', 'Binomial');
%     %         fittedData = [fittedData; [repmat(foldLoop, size(testData,1),1) repmat(noiseLoop, size(testData,1),1) table2array(testData(:,1:5)) predict(dataModel, testData)]];
%     fittedData((1:size(allData,1)/k)+size(allData,1)/k*(foldLoop-1),:) = [repmat(foldLoop, size(testData,1),1) table2array(testData(:,1:5)) predict(dataModel, testData)];
%     
% end
% %     dataModel
% %     size(fittedData)
% 
% 
% 
% % analyse k-fold data
% 
% % aggregate across folds
% foldAv = qStats(round(fittedData(:,7))==fittedData(:,6), fittedData(:,2:5),2);
% 
% % average over observers
% obsAv = [qStats(foldAv(:,5), foldAv(:,[2:4]))]
% 
% 
% 
% 
% %%
% totalTime = toc
% save('allFits_v05.mat')