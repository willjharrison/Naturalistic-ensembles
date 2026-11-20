clear
close all
clc

allData = prepareData;

%% show the distributions of difference between gabors-target, vs gabors-distractor

gtDiff = circ_dist(allData.targMean, allData.gabMean);
gdDiff = circ_dist(allData.distMean, allData.gabMean);

figure;
subplot(1,2,1)
hist(gtDiff,-pi:pi/12:pi)
hold;

axis([-pi pi 0 max(hist(gtDiff,-pi:pi/12:pi))])

subplot(1,2,2)
hist(gdDiff,-pi:pi/12:pi)
axis([-pi pi 0 max(hist(gtDiff,-pi:pi/12:pi))])

%% simple mean model
% what proportion correct would people get if they just used the mean?

simp = qStats(abs(gtDiff) < abs(gdDiff), [allData.participant, allData.block allData.sameDiff allData.nGabors],2);

simpleMeans = qStats(simp(:,5), simp(:,2:4),3)

%% simple SD model

gtSDDiff = allData.targSD - allData.gabSD;
gdSDDiff = allData.distSD - allData.gabSD;

simp = qStats(abs(gtSDDiff) < abs(gdSDDiff), [allData.participant, allData.block allData.sameDiff allData.nGabors],2);
simpleSD = qStats(simp(:,5), simp(:,2:4),3)

%%
simp = qStats((abs(gtDiff).^2 + abs(gtSDDiff).^2) <  (abs(gdDiff).^2 + abs(gdSDDiff).^2),[allData.participant, allData.block allData.sameDiff allData.nGabors],2);
out = qStats(simp(:,5), simp(:,2:4),3);

ensembleModel = splitvars(table(out),'out','NewVariableNames',{'block','sameDiff','nGabors','propCorrect','sem'});

% writetable(ensembleModel,'ensembleModel.csv')

%% weight both

x = allData.targPos - 1.5;
y = allData.chosenSquare - 1;

x1 = log(abs(gdDiff)./abs(gtDiff)).*x; 
x2 = log(abs(gdSDDiff)./abs(gtSDDiff)); 

x3 = (allData.block);
x4 = allData.sameDiff;
x5 = log(allData.nGabors);

data = table(x,x1,x2,x3,x4,x5,y);
f1 = fitglm(data,'y ~ x1', 'Distribution','binomial','Link','probit')
f1.Rsquared


%%
yHat = f1.predict;
faIdx = x < 0;
hIdx = x > 0;

% full summary - useful for error bars
faSummary = qStats(yHat(faIdx), [allData.block(faIdx) allData.nGabors(faIdx)],2)
hSummary = qStats(yHat(hIdx), [allData.block(hIdx) allData.nGabors(hIdx)],1)
% 
% sum = faSummary;
% sum(:,end) = norminv(hSummary) - norminv(faSummary(:,end))
% 

% a = qStats(round(f1.predict) == (x+.5), [allData.participant, allData.block allData.sameDiff allData.nGabors],2);
% a = qStats(a(:,5), a(:,2:4),3)

%% SDT

sensAll = qStats(y, [allData.targPos - 1.5 allData.participant allData.block allData.nGabors],2);

n = length(sensAll);

sensAll = [sensAll(1:n/2,2:end-1) norminv(sensAll(n/2+1:end,end)) - norminv(sensAll(1:n/2,end))];

qStats(sensAll(:,end), sensAll(:,2:end-1))

