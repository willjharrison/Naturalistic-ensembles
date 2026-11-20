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
leftMatch2 = leftMatch;
rightMatch2 = leftMatch;
leftOri = leftMatch;
rightOri = leftMatch;
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
        theta = theta - pi;
        
        %% gab En with spatial pooling
        
%         spatialPoolSD = .01;
%         noiseSD = 1;
%         gabEn = imgaussfilt(gabEn + noiseSD*std(gabEn(:))*randn(size(gabEn)),spatialPoolSD);
%         filtEn = imgaussfilt(filtEn + noiseSD*std(filtEn(:))*randn(size(gabEn)),spatialPoolSD);
%         distEn = imgaussfilt(distEn + noiseSD*std(distEn(:))*randn(size(gabEn)),spatialPoolSD);
        
%         if task.targPos(trial) == 1
%             leftMatch(trial + 600*(fileLoop-1),1) = [dot(((filtEn(:))),((gabEn(:))))];
%             rightMatch(trial + 600*(fileLoop-1),1) = [dot(((distEn(:))),((gabEn(:))))];
%             
%             
%         else
%             
%             leftMatch(trial + 600*(fileLoop-1),1) = [dot((((distEn(:)))),((gabEn(:))))];
%             rightMatch(trial + 600*(fileLoop-1),1) = [dot((((filtEn(:)))),((gabEn(:))))];
%         end
        
        if task.targPos(trial) == 1
            leftMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(filtEn, gabEn),[],[1 2]);
            rightMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(distEn, gabEn),[],[1 2]);
            
            
        else
            
            leftMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(distEn, gabEn),[],[1 2]);
            rightMatch(trial + 600*(fileLoop-1),1) = max(xcorr2(filtEn, gabEn),[],[1 2]);
        end
        
%         
        %% stim correlation 

        if task.targPos(trial) == 1
            leftMatch2(trial + 600*(fileLoop-1),1) = [dot(((targIm(:))),((gabIm(:))))];
            rightMatch2(trial + 600*(fileLoop-1),1) = [dot(((distIm(:))),((gabIm(:))))];
            
            
        else
            
            leftMatch2(trial + 600*(fileLoop-1),1) = [dot((((distIm(:)))),((gabIm(:))))];
            rightMatch2(trial + 600*(fileLoop-1),1) = [dot((((targIm(:)))),((gabIm(:))))];
        end
        
         %% with divisive normalisation
        
%         if task.targPos(trial) == 1
%             leftMatch2(trial + 600*(fileLoop-1),1) = [dot(((filtEn(:))./sum((filtEn(:)))),((gabEn(:))./sum((gabEn(:)))))];
%             rightMatch2(trial + 600*(fileLoop-1),1) = [dot(((distEn(:))./sum((distEn(:)))),((gabEn(:))./sum((gabEn(:)))))]; 
%             
%             
%         else
%             
%             leftMatch2(trial + 600*(fileLoop-1),1) = [dot(((distEn(:))./sum((distEn(:)))),((gabEn(:))./sum((gabEn(:)))))];
%             rightMatch2(trial + 600*(fileLoop-1),1) = [dot(((filtEn(:))./sum((filtEn(:)))),((gabEn(:))./sum((gabEn(:)))))];
%             
% 
%         end
        
        if any([isinf(leftMatch) isnan(leftMatch)],[1 2])
            any([isinf(leftMatch) isnan(leftMatch)],[1 2])
            return
        end
      
        %% find orientation similarity
        if task.targPos(trial) == 1
 
            leftOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(imOris(:),gabOris(:)).^2))];
            rightOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(distOris(:),gabOris(:)).^2))];
        %             leftOri = [leftOri; circ_mean(circ_dist(imOris(:),gabOris(:)))];
%             rightOri = [rightOri; circ_mean(circ_dist(distOris(:),gabOris(:)))];

        else
            
            leftOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(distOris(:),gabOris(:)).^2))];
            rightOri(trial + 600*(fileLoop-1),1) = [sqrt(sum(circ_dist(imOris(:),gabOris(:)).^2))];
            
            %             leftOri = [leftOri; circ_mean(circ_dist(distOris(:),gabOris(:)))];
%             rightOri = [rightOri; circ_mean(circ_dist(imOris(:),gabOris(:)))];

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

leftRightDiff2 = (leftMatch2 - rightMatch2);
logLeftMatch2 = (log(abs(leftMatch2)).*sign(leftMatch2));
logRightMatch2 = (log(abs(rightMatch2)).*sign(rightMatch2));
leftRightRatio2 = (logLeftMatch2-logRightMatch2);


leftOriTrans = log(leftOri);%.*(1/2);
rightOriTrans = log(rightOri);%.*(1/2);
leftRightOriDiff = (leftOri - rightOri);
leftRightOriRatio = leftOri./rightOri;

allData = table(participant,block,sameDiff,nGabors,targInt,resp,... % experiment factors
    targCat, distCat, ... % random effects of image category
    leftMatch,rightMatch,leftRightDiff,leftRightDiff2,leftMatch2,rightMatch2,leftRightRatio, leftRightRatio2,... % spatial image statistics
    leftOri, rightOri,leftOriTrans, rightOriTrans, leftRightOriDiff, leftRightOriRatio); % orientation image statistics

allData.absLeftOri = (abs(allData.leftOri));
allData.absRightOri = (abs(allData.rightOri));
allData.absLeftOri(allData.absLeftOri>pi/2) = -allData.absLeftOri(allData.absLeftOri>pi/2) + pi;
allData.absRightOri(allData.absRightOri>pi/2) = -allData.absRightOri(allData.absRightOri>pi/2) + pi;
allData.absLeftOri = (allData.absLeftOri).^(1/2);
allData.absRightOri = (allData.absRightOri).^(1/2);


qStats([allData.resp allData.resp==allData.targInt], table2array(allData(:,2:5)),2)

%%

% modSpec = ['leftRightDiff ~ 1  + targInt + '... % intercepts
%     'targInt:block + targInt:nGabors + targInt:sameDiff + '... % main effects
%     'targInt:block:nGabors + targInt:block:sameDiff + targInt:sameDiff:nGabors +'... % two way interactions
%     'targInt:block:nGabors:sameDiff + (1 | participant)']; % three-way and random effects
% imStatModel = fitglme(allData, modSpec);
% 
% mus = qStats(imStatModel.predict, table2array(allData(:,2:4)),2)
%%
clc
modSpec = ['targInt ~ leftRightOriRatio + leftMatch2*rightMatch2'];
dataModel = fitglme(allData, modSpec,'Distribution', 'Binomial','Link','probit')
allData.modelResp = dataModel.predict;
modPerformance = [qStats(round(allData.modelResp) == allData.targInt, table2array(allData(:,1:4)),2)];
humanPerformance = [qStats(allData.resp == allData.targInt, table2array(allData(:,1:4)),2)];

% regress model onto human
x1 = modPerformance(:,5);
y = humanPerformance(:,5);
g1 = humanPerformance(:,1);
dat = table(x1,y,g1);

pred = 'x1^2';
% f1 = fitglme(dat, ['y ~ ' pred ' + (' pred ' | g1)'])
f1 = fitglme(dat, ['y ~ ' pred ' + (1 | g1)'])
fittedModel = f1.predict; %x1*-.7652 + (x1.^2)*.782 + .7050;% f1.predict;

% summarise
results = qStats(fittedModel, humanPerformance(:,2:4))
[results(results(:,1) == 1 & results(:,2) == 0,4) results(results(:,1) == 1 & results(:,2) ==1,4) results(results(:,1) == 1 & results(:,2) == 0,4) - results(results(:,1) == 1 & results(:,2) ==1,4)]

%%
totalTime = toc
save('allFits.mat')