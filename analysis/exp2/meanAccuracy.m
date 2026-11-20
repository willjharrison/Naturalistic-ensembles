function [individualData,meanData] = meanAccuracy(allData)

% average within each participant
[m,c] = grpstats(allData.accuracy, ...
    [allData.participant, ...
    allData.block allData.thingOrMat, allData.nGabors],{'mean','gname'});
c = str2double(c);
individualData = [c m];
individualData = array2table(individualData,'VariableNames',{'participant','block','thingOrMat','nGabors','meanAccuracy'});

writetable(individualData,'individualData.csv')


%%
% average over participants
[mu,c,se] = grpstats(m, c(:,2:end),{'mean','gname','sem'});
meanData = [str2double(c) mu se];
meanData = array2table(meanData,'VariableNames',{'block','thingOrMat','nGabors','meanAccuracy','standardError'});

writetable(meanData,'meanData.csv')


%% get data sorted for JASP
% grab only relevant data
jData = individualData(:,[1 2 3 4 5]);
% chuck in a table without proper labels
JaspfullAnova = array2table(reshape(jData.meanAccuracy,20,10)');
% use PTB to create appropriately factorised labels
[c1,c2,c3] = BalanceFactors(1,0,{'rand_','aligned_'},{'object_','material_'},string(unique(allData.nGabors)));
condNames = strcat(c1,c2,c3);
% rename data in table and save
JaspfullAnova = renamevars(JaspfullAnova, 1:20, condNames); 
writetable(JaspfullAnova,'JaspfullAnova.csv');

%% break data down into 2 lots of 2 x 5 ANOVAs

randoData = individualData(individualData.block == 0,[1 3 4 5]);

randANOVAjasp = array2table(reshape(randoData.meanAccuracy,10,10)');
randANOVAjasp = renamevars(randANOVAjasp, [1:5], append('object_', string(unique(allData.nGabors))));
randANOVAjasp = renamevars(randANOVAjasp, 6:10, append('material_', string(unique(allData.nGabors))));
writetable(randANOVAjasp,'randANOVAjasp.csv');

alignedData = individualData(individualData.block == 1,[1 3 4 5]);
alignedANOVAjasp = array2table(reshape(alignedData.meanAccuracy,10,10)');
alignedANOVAjasp = renamevars(alignedANOVAjasp, 1:5, append('object_', string(unique(allData.nGabors))));
alignedANOVAjasp = renamevars(alignedANOVAjasp, 6:10, append('material_', string(unique(allData.nGabors))));
writetable(alignedANOVAjasp,'alignedANOVAjasp.csv');


%% set up data for checking if participants are better than chance in either aligned or randomised condition
[m,c] = grpstats(allData.accuracy, ...
    [allData.participant, ...
    allData.block],{'mean','gname'});
c = str2double(c);
individualData = [c m];

JaspMainCondsVsChance = splitvars(table(reshape(individualData(:,3), 2,length(unique(individualData(:,1))))'));
writetable(JaspMainCondsVsChance,'JaspMainCondsVsChance.csv');


%% some modelling examples
% x = allData.targPos - 1.5;
% y = allData.chosenSquare - 1;
% 
% x1 = x.*allData.block;
% x2 = x.*allData.sameDiff;
% x3 = x.*allData.nGabors;
% g1 = allData.participant;
% 
% xTable = table(x,x1,x2,x3,y,g1);
% fit = fitglme(xTable, 'y ~ x + x1 + (1 | g1)','Distribution', 'Binomial', 'Link','probit')

