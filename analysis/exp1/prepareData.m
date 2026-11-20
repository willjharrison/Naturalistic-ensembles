function allData = prepareData

files = dir('../Output/');
files = files(~ismember({files.name}, {'.','..'}));
files = files(~contains({files.name}, {'imStats', 'wrong','DS','model'}));
nFiles = length(files);

for i = 1:nFiles
    
    initials{i} = files(i).name(23:24);
    partCode = str2double(initials{i});
    
    load([files(i).folder '/' files(i).name],'task', 'stimuli','dv','trial')    
    load(['imAnalysis/stats/' files(i).name(1:end-4) '_stats.mat' ],'stats')
    load(['SCmodel/SCoutput/' files(i).name(1:end-4) '_SCmodel_gf_7.mat'],'simData')
        
    if i == 1
        
        allData.participant = repmat(partCode,trial,1);
        allData.block = ones(trial,1)*strcmp(files(i).name(17:21),'align');
        allData.nGabors = task.nGabors;
        allData.sameDiff = task.sameDiff;
        allData.targPos = task.targPos;
        allData.allCats = stimuli.allCats;
        allData.imRef = stimuli.imRef;
        allData.accuracy = dv.accuracy;
        allData.chosenSquare = dv.chosenSquare;
        allData.responseTime = dv.responseTime;
        allData.targMean = stats.targMean;
        allData.targSD = stats.targSD;
        allData.distMean = stats.distMean;
        allData.distSD = stats.distSD;
        allData.gabMean = stats.gabMean;
        allData.gabSD = stats.gabSD;
        allData.gabTargSim = simData.gabTargFFTSim;
        allData.gabDistSim = simData.gabDistFFTSim;

    else
           
        allData.participant = [allData.participant;repmat(partCode ,trial,1)];
        allData.block = [allData.block; ones(trial,1)*strcmp(files(i).name(17:21),'align')];
        allData.nGabors = [allData.nGabors;task.nGabors];
        allData.sameDiff = [allData.sameDiff;task.sameDiff];
        allData.targPos = [allData.targPos;task.targPos];
        allData.allCats = [allData.allCats;stimuli.allCats];
        allData.imRef = [allData.imRef;stimuli.imRef];
        allData.accuracy = [allData.accuracy;dv.accuracy];
        allData.chosenSquare = [allData.chosenSquare;dv.chosenSquare];
        allData.responseTime = [allData.responseTime;dv.responseTime];
        allData.targMean = [allData.targMean; stats.targMean];
        allData.targSD = [allData.targSD;stats.targSD];
        allData.distMean = [allData.distMean;stats.distMean];
        allData.distSD = [allData.distSD;stats.distSD];
        allData.gabMean = [allData.gabMean;stats.gabMean];
        allData.gabSD = [allData.gabSD;stats.gabSD];
        allData.gabTargSim = [allData.gabTargSim;simData.gabTargFFTSim];
        allData.gabDistSim = [allData.gabDistSim;simData.gabDistFFTSim];

    end
    
end

exp1Data = struct2table(allData);
writetable(exp1Data, 'exp1Data.csv');