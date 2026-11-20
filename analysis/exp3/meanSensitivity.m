% function dPrime = meanSensitivity(allData)
allData = prepareData;

% recode targ interval as [-.5 .5]
x = allData.targPos - 2 + .5;

% recode response as 0 1
y = allData.chosenSquare - 1;

g1 = allData.participant;

block = allData.block;

nGabs = log(allData.nGabors);

% data = table(x,y,block,g1,nGabs);
% 
% f1 = fitglme(data,'y ~ x + x:block + (x + x:block | g1)','Distribution','binomial', 'Link', 'probit')
% f2 = fitglme(data,'y ~ x + x:block + x:nGabs + (x + x:block + x:nGabs | g1)','Distribution','binomial', 'Link', 'probit')
% f3 = fitglme(data,'y ~ x:block + x:nGabs:block + (x:block + x:nGabs:block | g1)','Distribution','binomial', 'Link', 'probit')
% 
% %%
% f1fixef = f1.fixedEffects;
% f1ranef = f1.randomEffects;
% 
% b1Sens = f1fixef(2) + f1ranef(2:3:end);
% b2Sens = f1fixef(2) + f1fixef(3) + f1ranef(3:3:end);
% dMean = [b1Sens b2Sens]/sqrt(2);
% 
% [qStats([dMean(:,1)],zeros(10,1),3); qStats([dMean(:,2)],ones(10,1),3)]
% 
% %% CIs
% ci = f1.coefCI;
% ci = [ci(2,:); ci(2,:) + ci(3,:)]/sqrt(2)
% 
% %%
% 
% x1 = x(allData.block == 1);
% y1 = y(allData.block == 1);
% 
% data = table(x1,y1);
% fitglm(data,'y1 ~ x1 - 1','Distribution','binomial', 'Link', 'probit')
% xPred = 0
% predData = table(
% predict(f1, [0 0 1])

%% just use aligned data

idx = allData.block == 1;

x1 = x(idx);
nGabs1 = allData.nGabors(idx).^1;
y1 = y(idx);
g11 = g1(idx);
data2 = table(x1,y1,nGabs1,g11);

ff1 = fitglme(data2,'y1 ~ x1 + x1:nGabs1 + (x1 + x1:nGabs1 | g11)','Distribution','binomial', 'Link', 'probit')
%%
ff1fixef = ff1.fixedEffects;
d = (ff1fixef(2) + ff1fixef(3).*unique(nGabs1))/sqrt(2);

mu = qStats(allData.accuracy(idx), [allData.nGabors(idx)]);
dOld = (norminv(mu(:,2)) - norminv(1 - mu(:,2)))/sqrt(2);
pred = sqrt(unique(nGabs1)*d(1)^2); 
[d dOld pred]
