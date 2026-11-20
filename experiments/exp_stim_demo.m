clear
close all
clc

%% setup filters
% filtSize = ([imSize/128]);
stimuli.imSize = 256;
stimuli.filtSize = 8; %([2 4 8 16]);
nSFs = length(stimuli.filtSize);
allOris = [0 90];
nOris = length(allOris);

for sfLoop = 1:nSFs
    filt1(:,:,sfLoop) = gaussDerZero(stimuli.imSize,stimuli.filtSize(sfLoop),[stimuli.imSize/2 stimuli.imSize/2]); % /2 because PTB wants things on the -.5 to .5 scale
    filt2(:,:,sfLoop) = gaussDerNinety(stimuli.imSize,stimuli.filtSize(sfLoop),[stimuli.imSize/2 stimuli.imSize/2]);
    
end

filt(:,:,1,:) = fft2(filt1);
filt(:,:,2,:) = fft2(filt2);

% for later mapping, we set this up once
[masterY, masterX] = meshgrid(1:stimuli.imSize);

% create aperture
[angDist, radDist] = polarDistFun(stimuli.imSize);
mask = ...
    createLowPassFilter(radDist,stimuli.imSize/2 - 6,6);

% store coordinates of locations for indexing later
idx = 1:stimuli.imSize^2; % a list of every point in the image
[rows,cols] = ind2sub(stimuli.imSize,idx); % a list of all coords

stimuli.gaborSeparation = 0.75; % determines the closest distance between neighbour elements (but see specific line of code)
%% find all images
encodingGamma = 1;
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_concept_index.mat')
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_paths.mat')
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/unique_id.mat')

% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\image_concept_index.mat')
% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\image_paths.mat')
% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\unique_id.mat')

% thingsDir = 'C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Main\';

thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';


whichCat = 135; % 20 and 2 is good
whichIm = 1;
allTargIms = image_paths(image_concept_index == whichCat);
sourceIms = imresize((rgb2gray(double(imread([thingsDir allTargIms{whichIm}]))/255).^(encodingGamma))*2 - 1,[stimuli.imSize stimuli.imSize]);

distCat = 25;
distIm = 12;
allTargIms = image_paths(image_concept_index == distCat);
sourceIms(:,:,2) = imresize((rgb2gray(double(imread([thingsDir allTargIms{distIm}]))/255).^(encodingGamma))*2 - 1,[stimuli.imSize stimuli.imSize]);


sourceIms(sourceIms>1) = 1;
sourceIms(sourceIms<-1) = -1;
%     sourceIms = sourceIms.*mask;

% only im for filtering is always the first image - it's drawn randomly
% so it doesn't matter that it's always first
%     sourceIms = ((sourceIms/2+.5).^(1/encodingGamma))*2 - 1;
targIm = sourceIms(:,:,1);

%% filter image
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

%% do steered filter stuff
outIm = zeros(stimuli.imSize);

for sfLoop = 1:nSFs
    
    %% grab just this SF stuff
    nPlaces = 16; %task.nGabors(trial)/(2^(sfLoop-1));
    
    %         if task.noiseLevel(trial) == 0 % no change
    
    thisOris = imOris(:,:,sfLoop) + pi;
    
    %         elseif task.noiseLevel(trial) > 0 && task.meanOrVar(trial) == 1 % mean change
    %
    %             % randomly shift the mean change as positive or negative
    %             task.noiseLevel(trial) = task.noiseLevel(trial)*round(rand)*2 - 1;
    %             thisOris = imOris(:,:,sfLoop) + pi + deg2rad(task.noiseLevel(trial));
    %
    %         elseif task.noiseLevel(trial) > 0 && task.meanOrVar(trial) == 2 % variance change
    %
    %             thisOris = imOris(:,:,sfLoop) + pi + deg2rad(stimuli.checkerboard*task.noiseLevel(trial));
    %
    %         end
    
    % add noise to orientations
    %         orientationNoise = randn(size(thisOris))*orientationNoiseMultiplier;
    %         thisOris = thisOris + orientationNoise;
    sd = stimuli.filtSize(sfLoop);
    
    %% set up energy map for where to select places from
    thisSFMag = filtEn(:,:,sfLoop).*mask;
    
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
        
        gabCoords = [rows(thisPlace) cols(thisPlace)];
        
        % update energy map
        thisSFMag = thisSFMag.*(1 - Gaussian2D(sd*stimuli.gaborSeparation, gabCoords, stimuli.imSize));
        
        % in case random positions fall off canvas, clip Gabor pos
        gabCoords(gabCoords<1) = 1;
        gabCoords(gabCoords>stimuli.imSize) = stimuli.imSize;
        
        
        theta(locLoop,1) = thisOris(thisPlace);
        
        % update steered filter info
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

%% prepare stim as textures

trial = 1;
task.noiseLevel = .4;
stimuli.baseRMS = .2; % sets the background image contrast
outImNorm = maxContrast(outIm)/2*task.noiseLevel(trial);

% blend gabors either with source im or the non-source im
targImNorm = sourceIms(:,:,1).*mask;
targImNorm2 = sourceIms(:,:,2).*mask;
distractorIm = sourceIms(:,:,2).*mask;

targImNorm = targImNorm/std(targImNorm(:))/2*stimuli.baseRMS;
targImNorm2 = targImNorm2/std(targImNorm2(:))/2*stimuli.baseRMS;

blendedIm = (outImNorm) + .5;
blendedIm(blendedIm>1) = 1;
blendedIm(blendedIm<0) = 0;

blendedIm2 = (outImNorm + targImNorm2) + .5;
blendedIm2(blendedIm2>1) = 1;
blendedIm2(blendedIm2<0) = 0;

distImNorm = distractorIm/std(distractorIm(:))/2*stimuli.baseRMS + .5;
distImNorm(distImNorm > 1) = 1;
distImNorm(distImNorm < 0) = 0;
colWheel = colouriseOris(deg2rad(angDist));%imrotate(createColourWheel_v04(stimuli.imSize,0),45,'crop');
colWheel(repmat(radDist,1,1,3)>stimuli.imSize/2 | ...
    repmat(radDist,1,1,3)<stimuli.imSize/4) = 0;
figure;
imshow([ ...
    repmat((targIm/2+.5).^(1/encodingGamma),1,1,3) ...
    colouriseOris_withAmp(imOris,filtEn) ...
    colWheel ...
    repmat(blendedIm.^(1/encodingGamma),1,1,3) ...
    repmat(blendedIm2.^(1/encodingGamma),1,1,3) ...
    repmat((outImNorm/2 + .5).^(1/encodingGamma),1,1,3) ...
])


% cd ..
% cd Figures/
% % imwrite((sourceIms(:,:,1)/2 + .5).^(1/encodingGamma),'targIm.png')
% % % imwrite((sourceIms(:,:,2)/2 + .5).^(1/encodingGamma),'distIm.png')
% % imwrite(colouriseOris(imOris), 'oriIm.png')
% % imwrite(maxContrast(filtEn), 'energyIm.png')
% % imwrite(colouriseOris_withAmp(imOris,filtEn), 'ori_en_Im.png')
% % imwrite(outImNorm/2 + .5, 'filterIm.png')
% % imwrite(blendedIm.^(1/encodingGamma),'blendedIm.png')
% % imwrite(blendedIm2.^(1/encodingGamma),'dist_blendedIm.png')
% % imwrite(distImNorm.^(1/encodingGamma),'distIm.png')
% imwrite(colWheel, 'colWheel.png')
% imwrite(filt1/2 + .5, 'zeroFilt.png')
% imwrite(filt2/2 + .5, 'ninetyFilt.png')
% % 
% % % imwrite(distImNorm.^(1/encodingGamma),'distIm.png')
% % 
% cd ..
% cd Experiment