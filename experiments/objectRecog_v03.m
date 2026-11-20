clear 
close all
clc

% function objectRecog_v01
% this is built off steeredFilterPlace_v07_diffPlacement.m

Screen('Preference', 'SkipSyncTests', 1);
commandwindow;

% reset and shuffle seed
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);
randomSeed = rng;

%% get participant info
prompt = {'Enter participant code', 'Age', 'Gender', 'Corrected Vision?'};
dlg_title = 'Input';
num_lines= 1;
def     = {'XX','00','x','NA'};
participantInfo = inputdlg(prompt,dlg_title,num_lines,def);
subject_name = participantInfo{1};

% find out when this test is being run to add it to saved file names
newFileName=datestr(clock,30); % unique time stamp
[~,expName,~] = fileparts(mfilename); % name of this exp version

HideCursor;
ListenChar(2);

%% setup PTB

screens=Screen('Screens');
screenNumber=max(screens);

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma'); % gamma correction

[DrawingWindow, rect] = PsychImaging('OpenWindow', screenNumber, 127,[0 0 256*2 256*2],[],[],[],2); %[0 0 imSize*2 imSize*2],[],[],[],2);
% [DrawingWindow, rect] = PsychImaging('OpenWindow', screenNumber, 127,[],[],[],[],2); %[0 0 imSize*2 imSize*2],[],[],[],2);
screen.res = rect(3:4)';

Screen('ColorRange', DrawingWindow, 1, 1, 1);

[X_middle, Y_middle] = RectCenter(rect);
xy = [X_middle, Y_middle];
screen.centre =[X_middle; Y_middle];

% gamma correction
screen.gamma = 2;
PsychColorCorrection('SetEncodingGamma', DrawingWindow, 1./screen.gamma); % implement gamma correction
Screen('FillRect',DrawingWindow,.5);
Screen('Flip',DrawingWindow);

%% set up some basic experiment features
stimuli.imSize = 256;
stimuli.nItems = 2;

%% workout screen dimensions for tiling prototypical images
stimuli.protoRect = [1 1 stimuli.imSize stimuli.imSize];
stimuli.protoSize = stimuli.imSize;
stimuli.protoX = [1 stimuli.imSize] - stimuli.imSize/2 + xy(1); % centre x coords
stimuli.protoY = [xy(2) xy(2)];
stimuli.protoRect = CenterRectOnPoint(stimuli.protoRect, stimuli.protoX(:), stimuli.protoY(:));
% %%
% [stimuli.protoX,stimuli.protoY] = meshgrid(stimuli.protoCoords, stimuli.protoCoords);
% stimuli.protoX = stimuli.protoX + xy(1) - stimuli.imSize;
% stimuli.protoY = stimuli.protoY + xy(2) - stimuli.imSize/2;
% 
% stimuli.protoRect = CenterRectOnPoint(stimuli.protoRect, stimuli.protoX(:), stimuli.protoY(:));
stimuli.targRect = CenterRectOnPoint([1 1 stimuli.imSize stimuli.imSize],xy(1),xy(2));

% set up checkerboard for controlling variance
stimuli.checkerboard = (checkerboard(1,stimuli.imSize/2,stimuli.imSize/2)>.5)*2 - 1;

stimuli.strokeWidth = 14; % line width for cues/probe circle
stimuli.cueColour = [255 0 0]/255;

stimuli.targDur = .2;
stimuli.blankDur = .5;

stimuli.gaborSeparation = 0.75; % determines the closest distance between neighbour elements (but see specific line of code)
%% set up task design
noiseLevels = 0;% linspace(0,45,5);
noiseLevels = [noiseLevels noiseLevels(2:end)]; % need to repeat all numbers except 0 to balance n-trials per cond
gaborLevels = [1 4 16 64 256]; %2.^(2:6);
[task.meanOrVar, task.noiseLevel, task.nGabors, task.sameDiff] = BalanceFactors(25,1, 1:2, noiseLevels, gaborLevels,[0 1]); %[0 1 2 4 8],512, 1:25); %2.^(4:8),1:25);
task.nTrials = length(task.meanOrVar);
task.targPos = zeros(size(task.meanOrVar));

%% setup accuracy
dv.accuracy = 0;
dv.chosenSquare = zeros(size(task.meanOrVar));

%% setup filters
% filtSize = ([imSize/128]);
stimuli.filtSize = ([8]);
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
[masterY, masterX] = meshgrid(1:stimuli.imSize*2);

% create aperture
[~, radDist] = polarDistFun(stimuli.imSize);
mask = ...
    createLowPassFilter(radDist,stimuli.imSize/2 - 6,6);

% store coordinates of locations for indexing later
idx = 1:stimuli.imSize^2; % a list of every point in the image
[rows,cols] = ind2sub(stimuli.imSize,idx); % a list of all coords


%% find all images
encodingGamma = 2;
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_concept_index.mat')
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/image_paths.mat')
load('/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Variables/unique_id.mat')

% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\image_concept_index.mat')
% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\image_paths.mat')
% load('C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Variables\unique_id.mat')

% thingsDir = 'C:\Users\cogneuro_user\Desktop\Will Harrison\SteeredFilterPlacement\THINGS\Main\';

thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';

%% select all images
% pre-select concepts for every image
stimuli.allCats = randperm(length(unique_id),task.nTrials)'; % targets
distractorSet = 1:length(unique_id);
distractorSet(stimuli.allCats) = [];
chosenDistractors = distractorSet(randperm(length(distractorSet), task.nTrials))';
stimuli.allCats(:,2) = chosenDistractors;

for imAssignLoop = 1:task.nTrials
    
    % first concept is our target
    allTargIms = image_paths(image_concept_index == stimuli.allCats(imAssignLoop,1));
    
    % first index is the prototypical image we will show on response screen
    stimuli.imRef(imAssignLoop,1) = randperm(length(allTargIms),1);
    targIms{imAssignLoop,1} = allTargIms{stimuli.imRef(imAssignLoop,1)};
    
    % pick a distractor based on the same/diff condition for this trial
    if task.sameDiff(imAssignLoop) == 0 % same trial
        
        % select n-distractor concepts
        for i = 2:stimuli.nItems
            % just ignore the chosen distractor category, and just pick
            % another im from the target category
            allDistIms = image_paths(image_concept_index == stimuli.allCats(imAssignLoop,1));
            allDistIms(stimuli.imRef(imAssignLoop,1)) = [];
            
            stimuli.imRef(imAssignLoop,i) = randperm(length(allDistIms),1);
            targIms{imAssignLoop, i} = allDistIms{stimuli.imRef(imAssignLoop,i)};
        end
        
    elseif task.sameDiff(imAssignLoop) == 1 % diff trial
        
        % select n-distractor concepts
        for i = 2:stimuli.nItems
            allDistIms = image_paths(image_concept_index == stimuli.allCats(imAssignLoop,i));
            
            stimuli.imRef(imAssignLoop,i) = randperm(length(allDistIms),1);
            targIms{imAssignLoop, i} = allDistIms{stimuli.imRef(imAssignLoop,i)};
        end
    end
end

%% start trial loop
startTimer = GetSecs;

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
    outIm = zeros(stimuli.imSize*2);
    
    for sfLoop = 1:nSFs
        
        %% grab just this SF stuff
        nPlaces = task.nGabors(trial)%/(2^(sfLoop-1));
        
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
            imshow(thisSFMag,[])
            pause(.1)
            
            % then add distortion to gabor locations
%             gabCoords = [xNoise yNoise] + stimuli.imSize;
            
            % in case random positions fall off canvas, clip Gabor pos
            gabCoords(gabCoords<1) = 1;
            gabCoords(gabCoords>stimuli.imSize*2) = stimuli.imSize*2;
            
            
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
    
    %% dakin method to filter whole image
%     steeredIm = imgaussfilt(steeredImDakin(targIm,32,8).*mask/2 + .5,1);
%     steeredImNorm = normRMS(steeredIm,.05);
  
    
    
%     imshow(outIm,[])
    %% prepare stim as textures
    
    % max contrast and tidy bad values
%     outImNorm = normRMS(maxContrast(outIm)/2 + .5,.1);
    outImNorm = maxContrast(outIm)/2*.75 + .5;
    outImNorm(outImNorm>1) = 1;
    outImNorm(outImNorm<0) = 0;
    
    % create textures of all images
    targTex = Screen('MakeTexture',DrawingWindow,outImNorm,[],[],2);
    
    
    %% show targ pre-cue
    Screen('FrameOval',DrawingWindow,1, stimuli.targRect',stimuli.strokeWidth);
    Screen('Flip', DrawingWindow);
    
    WaitSecs(.2);
    
    %% draw targ, then n-AFC
    Screen('DrawTexture',DrawingWindow,targTex);
    Screen('Flip', DrawingWindow);
    
    % show targ followed by a delay
    WaitSecs(stimuli.targDur);
    
    Screen('Flip', DrawingWindow);
    WaitSecs(stimuli.blankDur);
    
    %% generate and show response screen
    stimuli.protoOrder(:,trial) = randperm(stimuli.nItems)';
    responseTex = Screen('MakeTexture',DrawingWindow,.5*ones([rect(4) rect(3)]),[],[],2);
    for i = 1:stimuli.nItems
        
        % keep track of where we put the target
        if stimuli.protoOrder(i,trial) == 1
            task.targPos(trial,1) = i;
        end
        
        thisIm = maxContrast(mask.*sourceIms(:,:,stimuli.protoOrder(i,trial)))/2 + .5;
        thisTex = Screen('MakeTexture', DrawingWindow, thisIm, [], [], 2);
        Screen('DrawTexture', responseTex, thisTex, [], stimuli.protoRect(i,:));
    end
    Screen('DrawTexture', DrawingWindow, responseTex);
    Screen('Flip', DrawingWindow);
    
    startRT = GetSecs;
    
    %% collect response
    takeResponse = -1;
    respSelected = 0;
    
    while takeResponse == -1
        
        [~,~,keyCode] = KbCheck(-1);
        if strcmp(KbName(keyCode),'LeftArrow')
            
            respSelected = 1;
            dv.chosenSquare(trial,1) = 1;
            
        elseif strcmp(KbName(keyCode),'RightArrow')
            
            respSelected = 1;
            dv.chosenSquare(trial,1) = 2;
            
        elseif strcmp(KbName(keyCode),'q')
            dv.accuracy(trial,1) = 0;
            break;
            
        end
        
        if respSelected == 1
            
            % highlight that grid location
%             Screen('DrawTexture', DrawingWindow, responseTex);
%             Screen('FrameOval',DrawingWindow,[100 255 100]/255, stimuli.protoRect(dv.chosenSquare(trial),:),stimuli.strokeWidth);
%             Screen('Flip',DrawingWindow);
            
%             [~,~,keyCode] = KbCheck(-1);
%             if strcmp(KbName(keyCode),'space')
                break
%             end
        end
        
    end
    
    dv.responseTime(trial,1) = GetSecs - startRT;
    dv.accuracy(trial,1) = task.targPos(trial) == dv.chosenSquare(trial);
    
%     SetMouse(xy(1),xy(2),DrawingWindow); % set the mouse to the random first orientation
%     takeResponse = -1; % keep the response  on a loop until a response is made
%     mouseMove = 0; % only show the mouse pointer and orientation if they move the mouse
%     while takeResponse == -1
%         
%         [mouse.X,mouse.Y,mouse.buttons] = GetMouse(DrawingWindow);
%         
%         % only show the mouse pointer and orientation if they move the mouse
%         if round(xy(1)) ~= mouse.X || round(xy(2)) ~= mouse.Y
%             mouseMove = 1;
%         end
%         
%         
%         if mouseMove == 1% show mouse-contingent probe
%             
%             % calculate mouse's grid location
%             mouse.gridDist = sqrt((stimuli.protoX - mouse.X).^2 + (stimuli.protoY - mouse.Y).^2);
%             [~,mouse.closestSquareIdx] = min(mouse.gridDist(:));
%             
%             ps = .5;
%             colDist1 = ((1 - (mouse.gridDist(:))/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^(1/ps))*127 + 127;
%             colDist2 = ((mouse.gridDist(:)/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^ps)*127;
%             colDist3 = ((mouse.gridDist(:)/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^ps)*127;
%             circStrokeWidth = ((1 - (mouse.gridDist(:))/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^(1/.25))*(stimuli.strokeWidth-2) + 1;
%             
%             % highlight that grid location
%             Screen('DrawTexture', DrawingWindow, responseTex);
%             Screen('FrameOval',DrawingWindow,[colDist2 colDist3 colDist1]'/255, stimuli.protoRect',circStrokeWidth);
%             Screen('DrawDots',DrawingWindow,[mouse.X mouse.Y],8,stimuli.cueColour);
%             Screen('Flip',DrawingWindow);
%         else
%             Screen('DrawTexture', DrawingWindow, responseTex);
%             Screen('Flip',DrawingWindow);
%         end
%         
%         if any(mouse.buttons) && mouseMove ~= 0
%             
%             takeResponse = 0;
%             firstClick = mouse.closestSquareIdx;
%             
%             while any(mouse.buttons)
%                 [mouse.X,mouse.Y,mouse.buttons] = GetMouse(DrawingWindow);
%             end
%             
%             while takeResponse == 0
%                 
%                 % calculate mouse's grid location
%                 [mouse.X,mouse.Y,mouse.buttons] = GetMouse(DrawingWindow);
%                 mouse.gridDist = sqrt((stimuli.protoX - mouse.X).^2 + (stimuli.protoY - mouse.Y).^2);
%                 [~,mouse.closestSquareIdx] = min(mouse.gridDist(:));
%                 colDist1 = ((1 - (mouse.gridDist(:))/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^(1/ps))*127 + 127;
%                 colDist2 = ((mouse.gridDist(:)/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^ps)*127;
%                 colDist3 = ((mouse.gridDist(:)/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^ps)*127;
%                 circStrokeWidth = ((1 - (mouse.gridDist(:))/(sqrt(min(rect(3))^2 + min(rect(4))^2))).^(1/.25))*(stimuli.strokeWidth-2) + 1;
%                 
%                 % highlight the selected circle in Green for the participant to
%                 % confirm
%                 Screen('DrawTexture', DrawingWindow, responseTex);
%                 Screen('FrameOval',DrawingWindow,[colDist2 colDist3 colDist1]'/255, stimuli.protoRect',circStrokeWidth);
%                 Screen('FrameOval',DrawingWindow,[100 255 100]/255, stimuli.protoRect(firstClick,:),stimuli.strokeWidth);
%                 Screen('DrawDots',DrawingWindow,[mouse.X mouse.Y],8,[0 255 0]/255);
%                 Screen('Flip',DrawingWindow);
%                 
%                 % quit the response routine when a response is made
%                 if any(mouse.buttons)
%                     secondClick = mouse.closestSquareIdx;
%                     
%                     if firstClick == secondClick
%                         
%                         dv.responseTime(trial,1) = GetSecs - startRT;
%                         
%                         HideCursor;
%                         takeResponse = 1;
%                         dv.chosenSquare(trial,1) = mouse.closestSquareIdx;
%                         
% %                         dv.chosenCeleb{trial} = dv.usableCelebs{stimuli.protoOrder(mouse.closestSquareIdx,trial)};
%                         dv.accuracy(trial,1) = mouse.closestSquareIdx == task.targPos(trial,1);
%                         
%                         [trial task.meanOrVar(trial) task.targPos(trial,1) mouse.closestSquareIdx dv.accuracy(trial)]
%                         
%                         %                     Eyelink('Message',['ResponseMade [' num2str(GetSecs*1000) ']']);
%                     else
%                         takeResponse = -1;
%                         
%                     end
%                     
%                 end
%                 
%                 % quit code
%                 [~,~,keyCode] = KbCheck(-1);
%                 if strcmp(KbName(keyCode),'q')
%                     dv.accuracy(trial,1) = 0;
%                     break;
%                 end
%                 
%             end
%         end
%         
%         % quit code
%         [~,~,keyCode] = KbCheck(-1);
%         if strcmp(KbName(keyCode),'q')
%             dv.accuracy(trial,1) = 0;
%             break;
%         end
%         
%     end
    
    %% give feedback    
    
    % show correct answer
    responseTex = Screen('MakeTexture',DrawingWindow,.5*ones([rect(4) rect(3)]),[],[],2);
    thisIm = maxContrast(mask.*sourceIms(:,:,1))/2 + .5;
    thisTex = Screen('MakeTexture', DrawingWindow, thisIm, [], [], 2);
    Screen('DrawTexture', responseTex, thisTex, [],stimuli.protoRect(task.targPos(trial,1),:));
    Screen('DrawTexture', DrawingWindow, responseTex);

    if dv.accuracy(trial,1) == 1
        
%         Screen('DrawTexture', DrawingWindow, responseTex);
        Screen('FrameOval',DrawingWindow,[100 255 100]/255, stimuli.protoRect(task.targPos(trial),:),stimuli.strokeWidth);
        Screen('Flip',DrawingWindow);
        WaitSecs(.5);
        
    % if an error, show correct im
    elseif dv.accuracy(trial,1) == 0
        
        
        Screen('FrameOval',DrawingWindow,stimuli.cueColour, stimuli.protoRect(task.targPos(trial,1),:),stimuli.strokeWidth);
        Screen('Flip',DrawingWindow);
        WaitSecs(1);
    end
    
    % clear the screen
    Screen('Flip', DrawingWindow);
    
    clear targTex responseTex thisTex
    
    % quit code
    [~,~,keyCode] = KbCheck(-1);
    if strcmp(KbName(keyCode),'q')
        break;
    end
    
    KbReleaseWait(-1);
    %%
end

%% tidy up
ListenChar(0)
sca
ShowCursor;


%% Output/Save file
newFileName=datestr(clock,30); % unique time stamp
experimentDuration = (GetSecs-startTimer)/60
cd Output/
save(sprintf('%s_%s_%s', expName,subject_name, newFileName))
cd ..




