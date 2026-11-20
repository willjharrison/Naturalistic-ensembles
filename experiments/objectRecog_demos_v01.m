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

thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';

saveDir = 'DemoFigures/';
%% setup im info
imSize = 256;

% create aperture
[~, radDist] = polarDistFun(imSize);
mask = ...
    createLowPassFilter(radDist,imSize/2 - 6,6);
%% load im

whichCat = 20; % 20 and 2 is good
whichIm = 4;

allTargIms = image_paths(image_concept_index == whichCat);
im = imresize((rgb2gray(double(imread([thingsDir allTargIms{whichIm}]))/255).^(encodingGamma))*2 - 1,[imSize imSize]);
padSize = 0;
% im = addPadding(mask.*im,padSize,0,0);
im(im<-1) = -1;
im(im>1) = 1;

figure;
imshow(im,[])

%%
nGabs = [1 4 16 64 256];
imSize = length(im); % update for padding

outIm = zeros(imSize,imSize,length(nGabs));
for i = 1:length(nGabs)
    
    outIm(:,:,i) = decompose_v01_rand(im, nGabs(i), 8, 1.5);
    
end

outIm(outIm(:)<0) = 0;
outIm(outIm(:)>1) = 1;

figure;

finalIm = sqrt([(mask.*im)/2 + .5 unpackIm(outIm,1)]);

imshow(finalIm)

%% save
imwrite(finalIm, sprintf([saveDir 'demoIm_%s_%s_rand.png'], num2str(whichCat), num2str(whichIm)))