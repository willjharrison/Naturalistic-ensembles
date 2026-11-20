function [individualData,meanData] = meanAccuracy(allData)

% average within each participant
[m,c] = grpstats(allData.accuracy, ...
    [allData.participant, ...
    allData.block allData.sameDiff, allData.nGabors],{'mean','gname'});
c = str2double(c);
individualData = [c m];
individualData = array2table(individualData,'VariableNames',{'participant','block','sameDiff','nGabors','meanAccuracy'});

writetable(individualData,'individualData.csv')


%%
% average over participants
[mu,c,se] = grpstats(m, c(:,2:end),{'mean','gname','sem'});
meanData = [str2double(c) mu se];
meanData = array2table(meanData,'VariableNames',{'block','sameDiff','nGabors','meanAccuracy','standardError'});

writetable(meanData,'meanData.csv')


%% get data sorted for JASP

% for jasp, we need this in the wide format, which will look something like
% this:
JaspfullAnova = splitvars(table(reshape(table2array(individualData(:,end)),20,max(allData.participant))'));
% forJasp.Properties.VariableNames = {sprintf('block %d \r', ( meanData.block))}
writetable(JaspfullAnova,'JaspfullAnova.csv');

%% break data down into 2 lots of 2 x 5 ANOVAs

randoData = individualData(individualData.block == 0,[1 3 4 5]);

randANOVAjasp = array2table(reshape(randoData.meanAccuracy,10,10)');
randANOVAjasp = renamevars(randANOVAjasp, 1:5, append('similar_', string(unique(allData.nGabors))));
randANOVAjasp = renamevars(randANOVAjasp, 6:10, append('dissimilar_', string(unique(allData.nGabors))));
writetable(randANOVAjasp,'randANOVAjasp.csv');

alignedData = individualData(individualData.block == 1,[1 3 4 5]);
alignedANOVAjasp = array2table(reshape(alignedData.meanAccuracy,10,10)');
alignedANOVAjasp = renamevars(alignedANOVAjasp, 1:5, append('similar_', string(unique(allData.nGabors))));
alignedANOVAjasp = renamevars(alignedANOVAjasp, 6:10, append('dissimilar_', string(unique(allData.nGabors))));
writetable(alignedANOVAjasp,'alignedANOVAjasp.csv');


%% set up data for checking if participants are better than chance in either aligned or randomised condition
[m,c] = grpstats(allData.accuracy, ...
    [allData.participant, ...
    allData.block],{'mean','gname'});
c = str2double(c);
individualData = [c m];

JaspMainCondsVsChance = splitvars(table(reshape(individualData(:,3), 2,10)'));
writetable(JaspMainCondsVsChance,'JaspMainCondsVsChance.csv');


% %% some modelling examples
% x = allData.targPos - 1.5;
% y = allData.chosenSquare - 1;
% 
% x1 = x.*(allData.block*2-1);
% x2 = x.*(allData.sameDiff*2-1);
% x3 = x.*allData.nGabors;
% g1 = allData.participant;
% 
% xTable = table(x,x1,x2,x3,y,g1);
% fit = fitglme(xTable, 'y ~ x1 + x2 + x3 - 1 + (x1 + x2 + x3 - 1 | g1)','Distribution', 'Binomial', 'Link','probit')

