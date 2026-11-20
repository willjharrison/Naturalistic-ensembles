function [negLL, g, fittedData, dataModel] = simData(allData, k, noiseLoops, x)

if any(x < 0)
    negLL = 10e6;
else
    negLL = 0;
    oriNoiseLevel = x(1);
    matchNoiseLevel = x(2);
    
    for nLoop = 1:noiseLoops
        
        %% apply noise
            allData.leftOriNoise = zscore(allData.leftOri + randn(size(allData.leftOri))*oriNoiseLevel);
            allData.rightOriNoise = zscore(allData.rightOri + randn(size(allData.leftOri))*oriNoiseLevel);
%         allData.leftRightOriDiffNoise = zscore(allData.leftRightOriDiff + randn(size(allData.leftRightOriDiff))*oriNoiseLevel);
        allData.leftMatchNoise = zscore(allData.leftMatch + randn(size(allData.leftMatch))*matchNoiseLevel);
        allData.rightMatchNoise = zscore(allData.rightMatch + randn(size(allData.rightMatch))*matchNoiseLevel);
%         allData.leftRightMatchDiff = zscore((allData.leftMatch - allData.rightMatch) + randn(size((allData.leftMatch - allData.rightMatch)))*matchNoiseLevel);
        %% use k-fold cross validation
        
        % grab memory
        fittedData = zeros(12000,8);
        
        % shuffle data before fitting so we remove sequential effects from
        % different folds
        shuffledData = allData(randperm(size(allData,1)),:);
        
        %     modSpec = ['targInt ~ leftOriNoise*rightOriNoise+leftMatchNoise+rightMatchNoise'];
        modSpec = ['targInt ~ leftOriNoise + rightOriNoise+leftMatchNoise+rightMatchNoise'];
        dataModel = cell(k,1);
        for foldLoop = 1:k
            
            testData = shuffledData((1:size(allData,1))<=size(allData,1)/k*foldLoop & (1:size(allData,1)) > size(allData,1)/k*(foldLoop -1),:);
            trainData = shuffledData(~((1:size(allData,1))<=size(allData,1)/k*foldLoop) | ~((1:size(allData,1)) > size(allData,1)/k*(foldLoop -1)),:);
            
            dataModel{foldLoop} = fitglme(trainData, modSpec,'Distribution', 'Binomial');
            %         fittedData = [fittedData; [repmat(foldLoop, size(testData,1),1) repmat(noiseLoop, size(testData,1),1) table2array(testData(:,1:5)) predict(dataModel, testData)]];
            fittedData((1:size(allData,1)/k)+size(allData,1)/k*(foldLoop-1),:) = [repmat(foldLoop, size(testData,1),1) table2array(testData(:,1:6)) predict(dataModel{foldLoop}, testData)];
            
        end
        
        %% simple fitting for aggregate data
%         % aggregate across folds
%         foldAv = qStats(round(fittedData(:,7))==fittedData(:,6), fittedData(:,2:5),2);
%         
%         % average over observers
%         obsAv = [qStats(foldAv(:,5), foldAv(:,[2:4]))];
%         
%         % compare with real data
%         g = fitglm(obsAv(:,4), meanAccuracy.meanAccuracy, 'Distribution','binomial');
        %% complex fitting for raw data
        x = norminv(fittedData(:,8));
        y = fittedData(:,7);
        g1 = fittedData(:,2);
        dat = table(x,y,g1);
%         g = fitglme(dat,'y ~ x^2 + (x^2 | g1)','Distribution','binomial');
        g = fitglm(dat,'y ~ x','Distribution','binomial');
        negLL = negLL + g.LogLikelihood;
    end
%     negLL
end

