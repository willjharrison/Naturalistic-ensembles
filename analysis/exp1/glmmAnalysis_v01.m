clear
close all
% clc


%% find files

files = dir('../Output/');
files = files(~ismember({files.name}, {'.','..','.DS_Store'}));
nFiles = length(files);

participant = [];
block = [];
nGabors = [];
sameDiff = [];
targInt = [];
resp = [];
for fileLoop = 1:nFiles
    
    load([files(fileLoop).folder '/' files(fileLoop).name],'task','dv','subject_name','mask','stimuli')
    thingsDir = '/Users/uqwharr1/Dropbox/Documents/Experiments/QBI/SteeredFilterPlacement/THINGS/Main/';
    
    % get the block data - is this aligned or randomised?
    whichBlock = contains(files(fileLoop).name,'align');
    
    participant = [participant; repmat(str2double(subject_name),600,1)];
    block = [block; repmat(whichBlock,600,1)];
    nGabors = [nGabors; task.nGabors];
    sameDiff = [sameDiff; task.sameDiff];
    targInt = [targInt; task.targPos - 1.5];
    resp = [resp; dv.chosenSquare-1];
    
end % end fileLoop

% gabCompression = (1/8);
% nGabors = (nGabors.^gabCompression);
nGabors = categorical(nGabors);
data = table(participant, block, sameDiff, nGabors, targInt, resp);

%%
modSpec = ['resp ~ 1  + '... % intercepts
    'targInt:block + targInt:nGabors + targInt:sameDiff + '... % main effects
    'targInt:block:nGabors + targInt:block:sameDiff + targInt:sameDiff:nGabors +'... % two way interactions
    'targInt:block:nGabors:sameDiff + (1 | participant)']; % three-way and random effects
dataModel = fitglme(data, modSpec,'Distribution', 'Binomial')

%%
idx = data.targInt < 0;
if isa(data.nGabors,'categorical')
    data.nGabors = double(string(data.nGabors)); %(data.nGabors);
else
    data.nGabors = data.nGabors.^(1/gabCompression);
end

dataArray = table2array(data);
modPredictions = dataModel.predict;

mus = qStats((modPredictions(idx)), dataArray(idx,1:4));

mus2 = qStats((modPredictions(~idx)), dataArray(~idx,1:4));


qStats(mean([(1 - mus(:,5)) mus2(:,5)],2), mus(:,2:4))
