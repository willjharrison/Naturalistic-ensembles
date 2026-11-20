clear
close all
clc

%% sim info
nReps = 10;
oriNoiseLevels = [linspace(0,.1,10)]';
posNoiseLevels = [logspace(1,log10(20),10)]';
[oriNoiseLevels, posNoiseLevels] = BalanceFactors(1,0,oriNoiseLevels,posNoiseLevels);

nNoiseLevels = length(oriNoiseLevels)

%% find all images
encodingGamma = 2;
thingsVariablesDir = '~/OneDrive - University of the Sunshine Coast 2/Useful Functions/ThingsVariables/';
load(fullfile(thingsVariablesDir, 'image_concept_index.mat'))
load(fullfile(thingsVariablesDir, 'image_paths.mat'))
load(fullfile(thingsVariablesDir, 'unique_id.mat'))
% load('../../../../SteeredFilterPlacement/THINGS/Variables/image_paths.mat')
% load('../../../../SteeredFilterPlacement/THINGS/Variables/unique_id.mat')

%% find files

files = dir('../../Output/');
files = files(~ismember({files.name}, {'.','..','.DS_Store'}));
nFiles = length(files);

simScore = @(imA, imB) dot(imA(:), imB(:))/(norm(imA(:))*norm(imB(:)));


obsAv = [];

rng(1042); % for reproducibility

% Preallocate variables to be filled in the parfor loop
allData = cell(nNoiseLevels,nReps);
startTime = GetSecs;
parfor noiseLoop = 1:nNoiseLevels  % Change this to parfor
    
    
    thisOriNoise = oriNoiseLevels(noiseLoop);
    thisPosNoise = posNoiseLevels(noiseLoop);
    
    for repLoop = 1:nReps
        thisRep = repLoop;
        
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
        gabSD = zeros(600,1);
        
        targCat = [];
        distCat = [];
        
        for fileLoop = 1:nFiles
            
            ff = load([files(fileLoop).folder '/' files(fileLoop).name],'task','dv','subject_name','stimuli');
            thingsDir = '~/OneDrive - University of the Sunshine Coast/Research/THINGSimages/';
            
            % get the block data - is this aligned or randomised?
            whichBlock = contains(files(fileLoop).name,'align');
            
            % store experiment variables
            participant = [participant; repmat(str2double(ff.subject_name),600,1)];
            block = [block; repmat(whichBlock,600,1)];
            nGabors = [nGabors; ff.task.nGabors];
            sameDiff = [sameDiff; ff.task.sameDiff];
            targInt = [targInt; ff.task.targPos - 1];
            targCat = [targCat; ff.stimuli.allCats(:,1)];
            distCat = [distCat; ff.stimuli.allCats(:,2)];
            resp = [resp; ff.dv.chosenSquare-1];
            
            
            %% start trial loop
            for trial = 1:ff.task.nTrials
                
                %% load image files
                stimFileName = sprintf(['modIms/gabIm_trial_%03d_fileLoop_%02d.mat'], trial, fileLoop);
                lf = load(stimFileName, 'gabIm', 'targIm', 'distIm','gabEn','filtEn','distEn','imOris','distOris','theta','gabOris');
                %
                % first make measures phase-invariant
                imOris = lf.imOris;
                imOris(imOris<0) = imOris(imOris<0) + pi;
                imOris = wrap(imOris*2);
                
                distOris = lf.distOris;
                distOris(distOris<0) = distOris(distOris<0) + pi;
                distOris = wrap(distOris*2);
                
                gabOris = lf.gabOris - pi; % correct for image generation (see above)
                gabOris(gabOris<0) = gabOris(gabOris<0) + pi;
                gabOris = wrap(gabOris*2);
                
                
                %% add orientation noise
                oriNoiseLevel = pi*thisOriNoise;
                imOris = imOris + randn(size(imOris))*oriNoiseLevel;
                distOris = distOris + randn(size(distOris))*oriNoiseLevel;
                gabOris = gabOris + randn(size(gabOris))*oriNoiseLevel;
                
                %% add position noise
                posNoiseLevel = thisPosNoise;
                filtEn = lf.filtEn + imgaussfilt(randn(size(lf.filtEn)),4)*posNoiseLevel;
                distEn = lf.distEn + imgaussfilt(randn(size(lf.distEn)),4)*posNoiseLevel;
                gabEn = lf.gabEn + imgaussfilt(randn(size(lf.gabEn)),4)*posNoiseLevel;
                
                %% gab En similarity
                
                if ff.task.targPos(trial) == 1
                    leftMatch(trial + 600*(fileLoop-1),1) = dot(filtEn(:), gabEn(:));% max(xcorr2(filtEn, gabEn),[],[1 2]);
                    rightMatch(trial + 600*(fileLoop-1),1) = dot(distEn(:), gabEn(:)); %max(xcorr2(distEn, gabEn),[],[1 2]);
                    
                    
                else
                    
                    leftMatch(trial + 600*(fileLoop-1),1) = dot(distEn(:), gabEn(:)); %max(xcorr2(distEn, gabEn),[],[1 2]);
                    rightMatch(trial + 600*(fileLoop-1),1) = dot(filtEn(:), gabEn(:)); % max(xcorr2(filtEn, gabEn),[],[1 2]);
                end
                
                
                %% orientation similarity
                if ff.task.targPos(trial) == 1
                    
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
                
                if ff.task.targPos(trial) == 1
                    
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
                
                if ff.task.targPos(trial) == 1
                    
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
                
                %                 if mod(trial,300) == 0
                %
                %                     if exist('l')
                %                         fprintf(repmat('\b',1,l))
                %                     end
                %                     l = fprintf('File: %s \rTrials: %s\rElapsed: %.1f mins\rProgress: %2.1f%%\r', files(fileLoop).name(17:24), string(trial), toc/60,(fileLoop-1)/nFiles*100);
                %                     pause(.1);
                %                 end
                
                
                %%
            end
            
            clc
            fprintf('Noiseloop: %d/%d \nRepLoop: %d/%d\nRuntime: %.1f mins\n\n', noiseLoop, nNoiseLevels, repLoop, nReps, (GetSecs - startTime)/60)
            drawnow('update');
            
        end % end fileLoop
        
        %% summarise predictors
        
        leftMatchTrans = zscore(log(abs(leftMatch)).*sign(leftMatch));
        rightMatchTrans = zscore(log(abs(rightMatch)).*sign(rightMatch));
        
        leftOriTrans = zscore(log(leftOri));
        rightOriTrans = zscore(log(rightOri));
        leftRightOriDiff = leftOriTrans - rightOriTrans;
        
        leftSDDiff = gabSD - leftSD;
        rightSDDiff = gabSD - rightSD;
        
        allData{noiseLoop,repLoop} = table(participant,block,sameDiff,nGabors,targInt,resp,targCat,distCat,... % experiment factors
            leftMatch,rightMatch,leftMatchTrans, rightMatchTrans,... % spatial-only image statistics
            leftOri, rightOri, leftOriTrans, rightOriTrans,leftRightOriDiff, leftOriMean, rightOriMean,... % spatially-dependent orientation image statistics leftMu, rightMu, gabSD, leftSD, rightSD, leftSDDiff, rightSDDiff,
            gabSD, leftMu, leftSD, leftSDDiff, rightMu, rightSD, rightSDDiff, leftMuUnweighted, rightMuUnweighted, leftSDUnweighted, rightSDUnweighted); % spatially-invariant orientation image statistics
        
        %         save(sprintf('noiseModelOutput_fine/allData_noiseLoop_%02d_repLoop_%02d.mat', ...
        %             noiseLoop,repLoop), ...
        %             'allData','oriNoiseLevels', 'posNoiseLevels')
%         fprintf(repmat('\b',1,60))
%         fprintf('Noiseloop: %d/%d \nRepLoop: %d/%d\nRuntime: %.1f mins\n\n', noiseLoop, nNoiseLevels, repLoop, nReps, toc/60)
%         drawnow('update');
        
    end
    
end
%%
totalTime = (GetSecs - startTime)/60;
fprintf('Total processing time was: %.1f mins\n', totalTime)
save('allFits_cells.mat','-v7.3')