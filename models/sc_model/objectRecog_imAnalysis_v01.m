clear
close all
clc


%% find all images
encodingGamma = 2;
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_concept_index.mat')
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_paths.mat')
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/unique_id.mat')

% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\image_concept_index.mat')
% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\image_paths.mat')
% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\unique_id.mat')

% thingsDir = 'C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Main\';


%% find files

files = dir('../../Output/');
files = files(~ismember({files.name}, {'.','..'}));
nFiles = length(files);

simScore = @(imA, imB) dot(imA(:), imB(:))/(norm(imA(:))*norm(imB(:)));

tic;
for fileLoop = 1:nFiles
    
    load([files(fileLoop).folder '/' files(fileLoop).name])
    thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';

    % get the block data - is this aligned or randomised?
    whichBlock = contains(files(fileLoop).name,'align');
    
    
    %% start trial loop
    
    for trial = 1:task.nTrials
        
        %% load images
        
        for imLoop = 1:stimuli.nItems
            
            sourceIms(:,:,imLoop) = imresize((rgb2gray(double(imread([thingsDir targIms{trial, imLoop}]))/255).^(encodingGamma))*2 - 1,[stimuli.imSize stimuli.imSize]);
            
        end
        
        sourceIms(sourceIms>1) = 1;
        sourceIms(sourceIms<-1) = -1;
        
        % only im for filtering is always the first image - it's drawn randomly
        % so it doesn't matter that it's always first
        targIm = ((sourceIms(:,:,1)/2+.5).^(1/encodingGamma))*2 - 1;
        distIm = ((sourceIms(:,:,2)/2+.5).^(1/encodingGamma))*2 - 1;
        
        %% filter target
        imFFT = fft2(targIm);
        for sfLoop = 1:nSFs
            for oriLoop = 1:nOris
                thisIm = (ifft2(imFFT.*filt(:,:,oriLoop,sfLoop),'symmetric'));
                filtIm(:,:,oriLoop,sfLoop) = thisIm;
            end
        end
        
        % organise filtered image
        for sfLoop = 1:nSFs
            for oriLoop = 1:nOris
                filtIm(:,:,oriLoop,sfLoop) = fftshift(filtIm(:,:,oriLoop,sfLoop));
                filtMag(:,:,oriLoop,sfLoop) = (filtIm(:,:,oriLoop,sfLoop)); % needs to be un-rectified for correct interpolation below
            end
        end
        
        % interpolate orientations at each pixel
        for sfLoop = 1:nSFs
            for oriLoop = 1:nOris
                thisOri = allOris(oriLoop);
                oriSin(:,:,oriLoop,sfLoop) = filtMag(:,:,oriLoop,sfLoop).*sin(deg2rad(thisOri));
                oriCos(:,:,oriLoop,sfLoop) = filtMag(:,:,oriLoop,sfLoop).*cos(deg2rad(thisOri));
            end
            imOris(:,:,sfLoop) = atan2(sum(oriSin(:,:,:,sfLoop),3),sum(oriCos(:,:,:,sfLoop),3));
            
        end
        
        % normalise energy and threshold
        
        for sfLoop = 1:nSFs
            filtEn(:,:,sfLoop) = abs(complex(filtMag(:,:,1,sfLoop),filtMag(:,:,2,sfLoop))).*mask;
            filtEnNorm(:,:,sfLoop) = filtEn(:,:,sfLoop)/sum(sum(filtEn(:,:,sfLoop)));
        end
        
        %% filter distractor
%         distIm = imrotate(kron([0;1],ones(stimuli.imSize/2,stimuli.imSize)),45,'crop');% sanity check
        imFFT = fft2(distIm);
        for sfLoop = 1:nSFs
            for oriLoop = 1:nOris
                thisIm = (ifft2(imFFT.*filt(:,:,oriLoop,sfLoop),'symmetric'));
                filtIm(:,:,oriLoop,sfLoop) = thisIm;
            end
        end
        
        % organise filtered image
        for sfLoop = 1:nSFs
            for oriLoop = 1:nOris
                filtIm(:,:,oriLoop,sfLoop) = fftshift(filtIm(:,:,oriLoop,sfLoop));
                filtMag(:,:,oriLoop,sfLoop) = (filtIm(:,:,oriLoop,sfLoop)); % needs to be un-rectified for correct interpolation below
            end
        end
        
        % interpolate orientations at each pixel
        for sfLoop = 1:nSFs
            for oriLoop = 1:nOris
                thisOri = allOris(oriLoop);
                oriSin(:,:,oriLoop,sfLoop) = filtMag(:,:,oriLoop,sfLoop).*sin(deg2rad(thisOri));
                oriCos(:,:,oriLoop,sfLoop) = filtMag(:,:,oriLoop,sfLoop).*cos(deg2rad(thisOri));
            end
            distOris(:,:,sfLoop) = atan2(sum(oriSin(:,:,:,sfLoop),3),sum(oriCos(:,:,:,sfLoop),3));
            
        end
        
        % normalise energy and threshold
        
        for sfLoop = 1:nSFs
            distEn(:,:,sfLoop) = abs(complex(filtMag(:,:,1,sfLoop),filtMag(:,:,2,sfLoop))).*mask;
            filtEnNorm(:,:,sfLoop) = filtEn(:,:,sfLoop)/sum(sum(filtEn(:,:,sfLoop)));
        end
        
        %% do steered filter stuff
        outIm = zeros(stimuli.imSize*2);
        
        for sfLoop = 1:nSFs
            
            %% grab just this SF stuff
            nPlaces = task.nGabors(trial);
            
            if task.noiseLevel(trial) == 0 % no change
                
                thisOris = imOris(:,:,sfLoop) + pi;
                
            elseif task.noiseLevel(trial) > 0 && task.meanOrVar(trial) == 1 % mean change
                
                % randomly shift the mean change as positive or negative
                task.noiseLevel(trial) = task.noiseLevel(trial)*round(rand)*2 - 1;
                thisOris = imOris(:,:,sfLoop) + pi + deg2rad(task.noiseLevel(trial));
                
            elseif task.noiseLevel(trial) > 0 && task.meanOrVar(trial) == 2 % variance change
                
                thisOris = imOris(:,:,sfLoop) + pi + deg2rad(stimuli.checkerboard*task.noiseLevel(trial));
                
            end
            
            % add noise to orientations
            %         orientationNoise = randn(size(thisOris))*orientationNoiseMultiplier;
            %         thisOris = thisOris + orientationNoise;
            sd = stimuli.filtSize(sfLoop);
            
            %% set up energy map for where to select places from
            thisSFMag = filtEn(:,:,sfLoop);
            
            % find total energy for normalisation and limit locations
            %         totalEn = sum(thisSFMag(:));
            %         thisSFMag = thisSFMag.*grid; % can limit number of possible Gab locations
            
            x = masterX;
            y = masterY;
            
            theta = [];
            allPlaces = [];
            
            %% determine places and add gabors
            
            for locLoop = 1:nPlaces
                
                [gabCont(locLoop,1),allPlaces(locLoop,1)] = max(thisSFMag(:));
                
                thisPlace = allPlaces(locLoop,1);
                
                %             [xNoise,yNoise] = pol2cart(rand*2*pi,rand*task.positionNoise(trial)*sd);
                %             [xNoise,yNoise] = pol2cart(rand*2*pi,0*sd);
                [xNoise,yNoise] = pol2cart(rand*2*pi, rand*stimuli.imSize/2); % task.positionNoise(trial)*sd);
                
                
                gabCoords = [rows(thisPlace) cols(thisPlace)];% + stimuli.imSize/2;
                
                % update energy map
                thisSFMag = thisSFMag.*(1 - Gaussian2D(sd*stimuli.gaborSeparation, gabCoords, stimuli.imSize));
%                 imshow(thisSFMag,[])
%                 pause(.1)
%                 

                if ~whichBlock % then add distortion to gabor locations
                    gabCoords = [xNoise yNoise] + stimuli.imSize/2;
                end

%                 % in case random positions fall off canvas, clip Gabor pos
                gabCoords(gabCoords<1) = 1;
                gabCoords(gabCoords>stimuli.imSize*2) = stimuli.imSize*2;
%                 
%                 
                theta(locLoop,1) = thisOris(thisPlace);
%                 
%                 % update steered filter info
                x = (x-gabCoords(1))/sd;
                y = (y-gabCoords(2))/sd;
                
                thisFilt = cos(theta(locLoop,1))*(-2*x.*exp(-(x.^2+y.^2))) + ...
                    sin(theta(locLoop,1))*(-2*y.*exp(-(x.^2+y.^2)));
                
                % reset x/y maps
                x = x*sd+gabCoords(1);
                y = y*sd+gabCoords(2);
                
                % modulate contrast - still can't decide the best method, but
                % switching off and then controlling global RMS is fine.
                %             thisFilt = thisFilt;% *gabCont(locLoop,1);
                
                outIm = outIm + thisFilt;
                
            end
            %         imshow(maxContrast(outIm)/2+.5);
            
        end
        
        %% finalise im as per exp (without creating texture for PTB)
        gabIm = maxContrast(outIm)*.75;
        gabIm = gabIm(1:stimuli.imSize, 1:stimuli.imSize);
        
        simData.gabTargSim(trial,1) = simScore(gabIm, targIm.*mask);
        simData.gabDistSim(trial,1) = simScore(gabIm, distIm.*mask);

        simData.gabTargEnSim(trial,1) = simScore(gabIm.^2, filtEn.^2.*mask);
        simData.gabDistEnSim(trial,1) = simScore(gabIm.^2, distEn.^2.*mask);

        
        simData.gabTargFFTSim(trial,1) = simScore(mask.*fftshift(abs(fft2(gabIm))), mask.*fftshift(abs(fft2(targIm.*mask))));
        simData.gabDistFFTSim(trial,1) = simScore(mask.*fftshift(abs(fft2(gabIm))), mask.*fftshift(abs(fft2(distIm.*mask))));
%         simData
        
%         gabTargEnergySim = simScore(gabIm.^2, filtEn'.^2)
    
%         [task.sameDiff(1:trial) task.nGabors(1:trial) simData.gabTargSim > simData.gabDistSim ...
%             simData.gabTargEnSim > simData.gabDistEnSim ...
%             simData.gabTargFFTSim > simData.gabDistFFTSim]
%         a = 1
       
%         
        if mod(trial,300) == 0
            
            sprintf('File: %s %.1f mins', files(fileLoop).name(17:24), toc/60)
            [m,g] = grpstats([simData.gabTargFFTSim > simData.gabDistFFTSim], [task.sameDiff(1:trial) task.nGabors(1:trial)],{'mean','gname'});
            [str2double(g) m]
        end
        
        %%
    end
    
    save(['SCoutput/' files(fileLoop).name(1:end-4) '_SCmodel.mat'], 'simData')
    clear simData
end % end fileLoop



