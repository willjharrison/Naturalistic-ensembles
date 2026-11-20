clear
close all
clc

thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/images/';
tic;
imds = imageDatastore(thingsDir, 'IncludeSubfolders', true, 'LabelSource', 'foldernames', 'FileExtensions', '.jpg');
toc
imSize = 32;
nIms = 500; % numel(imds.Files)
images = cell(1,nIms);

imOrder = randperm(numel(imds.Files));
whichIms = imOrder(1:nIms);
parfor i = 1:nIms%numel(imds.Files)
    images{i} = double(rgb2gray(imresize(readimage(imds,whichIms(i)),[imSize imSize])))/255;
end
toc

% imSize = 64;
% [images,rn] = read(imds);
% images = imresize(images, [imSize imSize]);
% data = reshape(images, [nIms/(imSize^2) imSize imSize]);
hiddenSize = 32;
autoenc = trainAutoencoder(images, hiddenSize, ...
            'L2WeightRegularization', 0.004, ...
            'SparsityRegularization', 4, ...
            'SparsityProportion', 0.15);


whichIms = imOrder(nIms + (1:nIms));
imagesChunk = cell(1,nIms);
parfor i = 1:nIms%numel(imds.Files)

    imagesChunk{i} = double(rgb2gray(imresize(readimage(imds,whichIms(i)),[imSize imSize])))/255;
    
end
autoenc2 = train(network(autoenc),imagesChunk,imagesChunk);
        
%         trainAutoencoder(images, 50, 'EncoderTransferFunction', 'logsig', 'DecoderTransferFunction', 'logsig');
toc

%%
xReconstructed = predict(autoenc, images);
figure;
for i = 1:20
    subplot(4,5,i);
    imshow(images{(i)});
end
figure;
for i = 1:20
    subplot(4,5,i);
    imshow(xReconstructed{i});
end

%%

figure;
for i = 1:hiddenSize
    subplot(4,round(hiddenSize/4),i);
    imshow(reshape(autoenc.EncoderWeights(i,:),imSize,imSize),[]);
end

% bigImShow(unpackIm(reshape(autoenc.EncoderWeights',imSize,hiddenSize*imSize),0,1),10,1);