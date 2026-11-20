function allData = prepareData

md = cd;

cd ../Output/

files = dir;
files = files(~ismember({files.name}, {'.','..','.DS_Store'}));
files = files(~contains({files.name}, {'exclusion'}));
nFiles = length(files);

for i = 1:nFiles
    
    initials{i} = files(i).name(28:29);
    partCode = str2double(initials{i});
    
    load(files(i).name, 'task', 'stimuli','dv','trial')
    
    if i == 1
        
        allData.participant = repmat(partCode,trial,1);
        allData.block = ones(trial,1)*strcmp(files(i).name(22:26),'align');
        allData.nGabors = task.nGabors;
        allData.thingOrMat = task.thingOrMat;
        allData.targPos = task.targPos;
        allData.allCats = stimuli.allCats;
        allData.imRef = stimuli.imRef;
        allData.accuracy = dv.accuracy;
        allData.chosenSquare = dv.chosenSquare;
        allData.responseTime = dv.responseTime;
        
    else
           
        allData.participant = [allData.participant;repmat(partCode ,trial,1)];
        allData.block = [allData.block; ones(trial,1)*strcmp(files(i).name(22:26),'align')];
        allData.nGabors = [allData.nGabors;task.nGabors];
        allData.thingOrMat = [allData.thingOrMat;task.thingOrMat];
        allData.targPos = [allData.targPos;task.targPos];
        allData.allCats = [allData.allCats;stimuli.allCats];
        allData.imRef = [allData.imRef;stimuli.imRef];
        allData.accuracy = [allData.accuracy;dv.accuracy];
        allData.chosenSquare = [allData.chosenSquare;dv.chosenSquare];
        allData.responseTime = [allData.responseTime;dv.responseTime];

    end
    
end

cd(md)