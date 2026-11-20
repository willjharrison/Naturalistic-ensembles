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
%         gabOris = gabOris - pi;
        
%         % smooth gabEn to simulate minor spatial imprecision
%         spatialSD = 2;
%         gabEn = imgaussfilt(gabEn,spatialSD);
%         filtEn = imgaussfilt(filtEn,spatialSD);
%         distEn = imgaussfilt(distEn,spatialSD);
%         
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
%         % add noise to the orientation estimates
%         oriNoiseSD = deg2rad(45);
%         imOris = imOris + randn(size(imOris))*oriNoiseSD;
%         distOris = distOris + randn(size(distOris))*oriNoiseSD;
%         gabOris = gabOris + randn(size(gabOris))*oriNoiseSD;
        
        %% gab En similarity 
        
        if task.targPos(trial) == 1
            leftMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(filtEn, gabEn),[],[1 2]);
            rightMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(distEn, gabEn),[],[1 2]);
            
            
        else
            
            leftMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(distEn, gabEn),[],[1 2]);
            rightMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(filtEn, gabEn),[],[1 2]);
        end
        
         
        %% orientation similarity
        if task.targPos(trial) == 1
 
            leftOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(imOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            rightOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(distOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
        %             leftOri = [leftOri; circ_mean(circ_dist(imOris(:),gabOris(~isnan(gabOris(:)))))];
%             rightOri = [rightOri; circ_mean(circ_dist(distOris(:),gabOris(~isnan(gabOris(:)))))];

        else
            
            leftOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(distOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            rightOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(imOris(~isnan(gabOris(:))),gabOris(~isnan(gabOris(:)))).^2))];
            
            %             leftOri = [leftOri; circ_mean(circ_dist(distOris(:),gabOris(~isnan(gabOris(:)))))];
%             rightOri = [rightOri; circ_mean(circ_dist(imOris(:),gabOris(~isnan(gabOris(:)))))];

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
%         clc
%         circ_std(imOris(:), filtEn(:))
%         circ_std(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:))))
%         circ_std(gabOris(~isnan(gabOris(:))), gabEn(~isnan(gabOris(:))))
%         circ_std(distOris(:), distEn(:))

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
%%%%%%
%% using log values doesn't work
%%%%%%
%         %% gab En
%         
%         spatialPoolSD = 2;
%         gabEn = imgaussfilt(gabEn,spatialPoolSD);
%         filtEn = imgaussfilt(filtEn,spatialPoolSD);
%         distEn = imgaussfilt(distEn,spatialPoolSD);
%         
%         if task.targPos(trial) == 1
%             leftMatch2 = [leftMatch2; dot(log((filtEn(:) + offset)),log((gabEn(:) + offset)))];
%             rightMatch2 = [rightMatch2; dot(log((distEn(:) + offset)),log((gabEn(:) + offset)))];
%             
%             
%         else
%             
%             leftMatch2 = [leftMatch2; dot((log((distEn(:) + offset))),log((gabEn(:) + offset)))];
%             rightMatch2 = [rightMatch2; dot((log((filtEn(:) + offset))),log((gabEn(:) + offset)))];
%         end
%         
%          %% with divisive normalisation
%         
%         if task.targPos(trial) == 1
%             leftMatch = [leftMatch; dot(log((filtEn(:) + offset)./sum((filtEn(:) + offset))),log((gabEn(:) + offset)./sum((gabEn(:) + offset))))];
%             rightMatch = [rightMatch; dot(log((distEn(:) + offset)./sum((distEn(:) + offset))),log((gabEn(:) + offset)./sum((gabEn(:) + offset))))]; 
%             
%         else
%             
%             leftMatch = [leftMatch; dot(log((distEn(:) + offset)./sum((distEn(:) + offset))),log((gabEn(:) + offset)./sum((gabEn(:) + offset))))];
%             rightMatch = [rightMatch; dot(log((filtEn(:) + offset)./sum((filtEn(:) + offset))),log((gabEn(:) + offset)./sum((gabEn(:) + offset))))];
%         end
%         
%         if any([isinf(leftMatch) isnan(leftMatch)],[1 2])
%             any([isinf(leftMatch) isnan(leftMatch)],[1 2])
%             return
%         end
        
        %%
        
        
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

noiseLevel = 0;
noisyLeft = leftMatch + randn(size(leftMatch))*std(leftMatch)*noiseLevel;
noisyRight = rightMatch + randn(size(rightMatch))*std(rightMatch)*noiseLevel;
logLeftMatch = (log(abs(leftMatch)).*sign(leftMatch));
logRightMatch = (log(abs(rightMatch)).*sign(rightMatch));
leftRightRatio = (leftMatch./rightMatch);
leftRightDiff = (leftMatch) - (rightMatch);

leftOriTrans = log(leftOri);%.*(1/2);
rightOriTrans = log(rightOri);%.*(1/2);
leftRightOriDiff = (leftOri - rightOri);
leftRightOriRatio = leftOri./rightOri;

noiseSD = 0;
leftMatchTrans = zscore(log(abs(leftMatch)).*sign(leftMatch));
rightMatchTrans = zscore(log(abs(rightMatch)).*sign(rightMatch));
% leftRightDiff = (leftMatchTrans) - (rightMatchTrans);
leftOriTrans = zscore(log(leftOri)) + randn(size(leftOri))*noiseSD;
rightOriTrans = zscore(log(rightOri)) + randn(size(leftOri))*noiseSD;
% leftRightDiff = (leftMatchTrans.^2 + leftOriTrans.^2) > (rightMatchTrans.^2 + rightOriTrans.^2);
leftRightDiff = (leftMatchTrans.^2) > (rightMatchTrans.^2);
leftRightOriDiff = leftOriTrans - rightOriTrans; %sqrt((leftOriTrans.*leftMatchTrans).^2) > sqrt((rightOriTrans.*rightMatchTrans).^2);

leftSDDiff = sqrt(gabSD.^2 + leftSD.^2);
rightSDDiff = sqrt(gabSD.^2 + rightSD.^2);
accuracy = resp == targInt;
allData = table(participant,block,sameDiff,nGabors,targInt,resp,accuracy,... % experiment factors
    targCat, distCat, ... % random effects of image category
    leftMatch,rightMatch,leftMatchTrans, rightMatchTrans,leftRightDiff,leftRightRatio,... % spatial image statistics
    leftOri, rightOri, leftRightOriDiff, leftOriTrans, rightOriTrans, leftMu, rightMu, gabSD, leftSD, rightSD, leftSDDiff, rightSDDiff, ...% orientation image statistics with weights
    leftMuUnweighted, rightMuUnweighted, leftSDUnweighted, rightSDUnweighted); % orientation image statistics without weights

%%
% 
% statDistUnweighted = sqrt(leftMuUnweighted.^2 + leftSDUnweighted.^2) > sqrt(rightMuUnweighted.^2 + rightSDUnweighted.^2);
% qStats(statDistUnweighted==allData.targInt,  table2array(allData(:,2:4)))
% 
% statDist = sqrt(leftMu.^2 + leftSD.^2) > sqrt(rightMu.^2 + rightSD.^2);
% a = qStats(statDist==allData.targInt,  table2array(allData(:,1:4)));
% qStats(a(:,5), a(:,2:4))

%% simple fit

% allData = allData(allData.block==1,:);
clc
% modSpec = ['targInt ~ leftOri*leftSD + rightMatch - 1']; % <- this works!
% modSpec = ['targInt ~ leftOri*leftSD - 1'];% <- this works really well with hardcore fit!
modSpec = ['targInt ~ leftOri'];
dataModel = fitglme(allData, modSpec,'Distribution', 'Binomial') %,'Link','probit')
allData.modelResp = dataModel.predict;

% regress model onto human
modPerformance = [qStats(round(allData.modelResp) == allData.targInt, table2array(allData(:,1:4)),2)];
humanPerformance = [qStats(allData.accuracy, table2array(allData(:,1:4)),2)];

x1 = modPerformance(:,5);
y = humanPerformance(:,5);
g1 = humanPerformance(:,1);
dat = table(x1,y,g1);

pred = 'x1^1';
f1 = fitglme(dat, ['y ~ ' pred ' + (' pred ' | g1)'])
% f1 = fitglme(dat, ['y ~ ' pred ' + (x1 | g1)'])
fittedModel = f1.predict; %x1*-.7652 + (x1.^2)*.782 + .7050;% f1.predict;

% summarise
results = qStats(fittedModel, humanPerformance(:,2:4))
[results(results(:,1) == 1 & results(:,2) == 0,4) results(results(:,1) == 1 & results(:,2) ==1,4) results(results(:,1) == 1 & results(:,2) == 0,4) - results(results(:,1) == 1 & results(:,2) ==1,4)]


%%%%%%
%% try left vs right stats in separate model, then as combined predictors
% modSpec = ['targInt ~ leftOriTrans*leftMatch - 1'];
% leftModel = fitglme(allData, modSpec,'Distribution', 'Binomial') %,'Link','probit')
% leftModelResp = leftModel.predict;
% 
% modSpec = ['targInt ~ rightOriTrans*rightMatch - 1'];
% rightModel = fitglme(allData, modSpec,'Distribution', 'Binomial') %,'Link','probit')
% rightModelResp = rightModel.predict;
% 
% x1 = leftModelResp;
% x2 = rightModelResp;
% y = allData.resp;
% g1 = allData.participant;
% dat = table(x1,x2,y,g1);
% fFit = fitglme(dat, 'y ~ x1*x2 + (1 | g1)', 'Distribution', 'Binomial') 
% 
% results = qStats(round(fFit.predict)==allData.targInt, table2array(allData(:,1:4)));
% results = qStats(results(:,5), results(:,2:4))
% [results(results(:,1) == 1 & results(:,2) == 0,4) results(results(:,1) == 1 & results(:,2) ==1,4) results(results(:,1) == 1 & results(:,2) == 0,4) - results(results(:,1) == 1 & results(:,2) ==1,4)]




%% cross fold



%% hardcore fit
% clc

% regress model onto human
x1 = allData.modelResp;
y = allData.resp;
g1 = allData.participant;
g2 = allData.targCat; 
g3 = allData.distCat;
dat = table(x1,y,g1,g2,g3);

pred = 'x1^1';
% f1 = fitglme(dat, ['y ~ ' pred ' + (' pred ' | g1)'])
f1 = fitglme(dat, ['y ~ ' pred ' + (x1 | g1)'])
fittedModel = round(f1.predict)==allData.targInt; %x1*-.7652 + (x1.^2)*.782 + .7050;% f1.predict;

% summarise
results = qStats(fittedModel, table2array(allData(:,1:4)));
results = qStats(results(:,5), results(:,2:4))
[results(results(:,1) == 1 & results(:,2) == 0,4) results(results(:,1) == 1 & results(:,2) ==1,4) results(results(:,1) == 1 & results(:,2) == 0,4) - results(results(:,1) == 1 & results(:,2) ==1,4)]

%%
totalTime = toc
save('allFits_v04.mat')