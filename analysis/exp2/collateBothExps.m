clear
close all
clc

exp2Data = (prepareData);
exp2Data.participant = exp2Data.participant+1;
exp2Data = rmfield(exp2Data, 'allCats');
exp2Data = splitvars(struct2table(exp2Data));
td = cd;

cd ../../Analysis
exp1Data = prepareData; %splitvars(struct2table(prepareData));
exp1Data = rmfield(exp1Data, 'allCats');
exp1Data = splitvars(struct2table(exp1Data));
cd(td)
%%
% re-label thingOrMat so we can concatenate with sameDiff from Exp 1
exp2Data.thingOrMat(exp2Data.thingOrMat==1) = exp2Data.thingOrMat(exp2Data.thingOrMat==1) + 2;
exp2Data = renamevars(exp2Data, 'thingOrMat','sameDiff');

allData = exp1Data;%(:,[1:3 8]);
allData = [allData; exp2Data];%(:,[1:3 8])];

%%
% average within each participant
[m,c] = grpstats(allData.accuracy, ...
    [allData.participant, ...
    allData.block allData.sameDiff, allData.nGabors],{'mean','gname'});
c = str2double(c);
individualData = [c m];
individualData = array2table(individualData,'VariableNames',{'participant','block','sameDiff','nGabors','meanAccuracy'});

% writetable(individualData,'individualData.csv')


%%
% average over participants
[mu,c,se] = grpstats(m, c(:,2:end),{'mean','gname','sem'});
meanData = [str2double(c) mu se];
meanData = array2table(meanData,'VariableNames',{'block','sameDiff','nGabors','meanAccuracy','standardError'})
% meanData.sameDiff=meanData.sameDiff
writetable(meanData,'collatedData_meanData.csv')
