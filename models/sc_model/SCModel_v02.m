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

tic;
for fileLoop = 1:nFiles
    
    load([files(fileLoop).folder '/' files(fileLoop).name],'task','dv','subject_name','mask','stimuli')
    thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';

    % get the block data - is this aligned or randomised?
    whichBlock = contains(files(fileLoop).name,'align');   
    
    %% start trial loop
    
    for trial = 1:task.nTrials
        
        %% load image files
        stimFileName = sprintf(['modIms/gabIm_trial_%03d_fileLoop_%02d.mat'], trial, fileLoop);
        load(stimFileName, 'gabIm', 'targIm', 'distIm','gabEn','filtEn','distEn','imOris','distOris','theta','gabOris');
        theta = theta - pi; 
        
        %%
   
        normPoolSD = 16;
        noiseSD = 2;
        
%         gabStim = gabEn./imgaussfilt(gabEn,normPoolSD) + noiseSD*randn(size(gabEn));
%         targStim = filtEn./imgaussfilt(filtEn,normPoolSD) + noiseSD*randn(size(gabEn));
%         distStim = distEn./imgaussfilt(distEn,normPoolSD) + noiseSD*randn(size(gabEn));
%         
%         gabStim = gabEn./imgaussfilt(gabEn,normPoolSD);% + noiseSD*randn(size(gabEn));
%         targStim = filtEn./imgaussfilt(filtEn,normPoolSD);% + noiseSD*randn(size(gabEn));
%         distStim = distEn./imgaussfilt(distEn,normPoolSD);% + noiseSD*randn(size(gabEn));
% 
%         % correct nans due to divide by zero error
%         gabStim(isnan(gabStim(:))) = 0;
%         targStim(isnan(targStim(:))) = 0;
%         distStim(isnan(distStim(:))) = 0;
        
        gabStim = gabEn + noiseSD*std(gabEn(:))*randn(size(gabEn));
        targStim = filtEn + noiseSD*std(filtEn(:))*randn(size(gabEn));
        distStim = distEn + noiseSD*std(distEn(:))*randn(size(gabEn));
        
        
        targDist = corr(gabIm(:), targIm(:).*mask(:));
        distDist = corr(gabIm(:), distIm(:).*mask(:));
        
        
        targOri = circ_mean(circ_dist(imOris(:),gabOris(:)),gabStim(:).*targStim(:));
        distOri = circ_mean(circ_dist(distOris(:),gabOris(:)),gabStim(:).*distStim(:));
  
        simTest(trial,:,fileLoop) = [str2double(subject_name) whichBlock task.sameDiff(trial) task.nGabors(trial) dv.accuracy(trial) ...
            (dot(zscore(targStim(:)),zscore(gabStim(:)))) (dot(zscore(distStim(:)),zscore(gabStim(:)))) ...
            targDist distDist ...
            targOri distOri];
%         
%         [trial simTest(trial,[3 6:end],fileLoop)]
%         
%         imshow([gabIm targIm.*mask distIm.*mask gabIm + targIm.*mask gabIm + distIm.*mask],[])
%         pause(.1);
        
        if mod(trial,300) == 0
            
            if exist('l')
                fprintf(repmat('\b',1,l))
            end
            l = fprintf('File: %s \rTrials: %s\rElapsed: %.1f mins\rProgress: %2.1f%%\r', files(fileLoop).name(17:24), string(trial), toc/60,(fileLoop-1)/nFiles*100);
            pause(.1);
%             [m,g] = grpstats([simData.gabTargFFTSim > simData.gabDistFFTSim], [task.sameDiff(1:trial) task.nGabors(1:trial)],{'mean','gname'});
%             [str2double(g) m]
        end
        
        
        %%
    end
%     fprintf(repmat('\b',1,l))
%     fprintf('Progress: %2.1f%%\r', fileLoop/nFiles*100);
        
%     save(['SCoutput/' sprintf([files(fileLoop).name(1:end-4) '_SCmodel_gf_%s.mat'], string(gSize))], 'simData')
    clear simData
end % end fileLoop

save('simTest.mat','simTest')
%%
allData = reshape(permute(simTest,[3 1 2]),12000,size(simTest,2));

allData(:,12:14) = allData(:,[6:2:end])./allData(:,[7:2:end]);


% allData(:,end+1) = allData(:,end) > 1;
%%
imStats = qStats(abs(allData(:,[12:end])), allData(:,2:4),2)

%%
enStat = reshape(imStats(11:end,4),5,2)
stimStat = reshape(imStats(11:end,5),5,2)
oriStat = reshape(imStats(11:end,6),5,2)

%%

participant = allData(:,1);
block = allData(:,2);
sameDiff = allData(:,3);
nGabors = allData(:,4);
accuracy = allData(:,5);
targCorr = sqrt(abs(allData(:,6)));
distCorr = sqrt(abs(allData(:,7)));
absTargOri = (abs(allData(:,10)));
absDistOri = (abs(allData(:,11)));
absTargOri(absTargOri>pi/2) = -absTargOri(absTargOri>pi/2) + pi;
absDistOri(absDistOri>pi/2) = -absDistOri(absDistOri>pi/2) + pi;
absTargOri = sqrt(absTargOri);
absDistOri = sqrt(absDistOri);

corrRatio = abs(targCorr)./abs(distCorr);
targDist = (allData(:,8));
distDist = (allData(:,9));
distRatio = targDist - distDist;
% corrDiff = double(targCorr > distCorr);

data = table(participant,block,sameDiff,nGabors,accuracy,targCorr,distCorr, corrRatio, targDist, distDist, distRatio, absTargOri, absDistOri);
% data = data(data.block ==1,:);
% f1 = fitglme(data, 'accuracy ~ targCorr*distCorr + (1 | participant)', 'Distribution', 'Binomial')
f2 = fitglme(data, 'accuracy ~ targCorr*distCorr*targDist*distDist + (1 | participant)', 'Distribution', 'Binomial')
f3 = fitglme(data, 'accuracy ~ targCorr*distCorr*targDist*distDist*absTargOri^2*absDistOri^2 + (1 | participant)', 'Distribution', 'Binomial')

f1 = fitglme(data, 'accuracy ~ targCorr^2*distCorr^2 + (1 | participant)', 'Distribution', 'Binomial')

%%
mus = qStats((f1.predict), table2array(data(:,1:4)));
qStats(mus(:,5), mus(:,2:4))

if exist('f2','var')
mus = qStats((f2.predict), table2array(data(:,1:4)));
qStats(mus(:,5), mus(:,2:4))
end

if exist('f3','var')
mus = qStats((f3.predict), table2array(data(:,1:4)));
qStats(mus(:,5), mus(:,2:4))
end
