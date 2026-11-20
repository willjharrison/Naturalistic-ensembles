clear
close all
clc

allData = prepareData;

%% summarise mean/SD info

gtDiff = circ_dist(allData.targMean, allData.gabMean);
gdDiff = circ_dist(allData.distMean, allData.gabMean);

gtSDDiff = allData.targSD - allData.gabSD;
gdSDDiff = allData.distSD - allData.gabSD;

%% combine mean and sd info to guess response on each trial
simp = qStats((abs(gtDiff).^2 + abs(gtSDDiff).^2) <  (abs(gdDiff).^2 + abs(gdSDDiff).^2),[allData.participant, allData.block allData.sameDiff allData.nGabors],2);
out = qStats(simp(:,5), simp(:,2:4),3);

ensembleModel = splitvars(table(out),'out','NewVariableNames',{'block','sameDiff','nGabors','propCorrect','sem'})

% writetable(ensembleModel,'ensembleModel.csv')
