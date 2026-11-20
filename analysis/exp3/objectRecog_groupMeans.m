
md = cd;

cd ../Output/

files = dir;
files = files(~ismember({files.name}, {'.','..'}));
files = files(~contains({files.name}, {'imStats', 'wrong','DS','model'}));
nFiles = length(files);

files.name


partCode = 1;
for i = 1:nFiles
    
    initials{i} = files(i).name(23:24);
    
    if i > 1
        if ~strcmp(initials{i},initials{i -1})
            partCode = partCode + 1;
        end
    end
    
    load(files(i).name, 'task', 'stimuli','dv')
    
    if i == 1
        
        allData.participant = repmat(partCode,trial,1);
        allData.noiseLevel = task.noiseLevel;
        allData.nGabors = task.nGabors;
        allData.sameDiff = task.sameDiff;
        allData.targPos = task.targPos;
        allData.allCats = stimuli.allCats;
        allData.imRef = stimuli.imRef;
        allData.similarity = similarity;
        allData.accuracy = dv.accuracy;
        allData.chosenSquare = dv.chosenSquare;
        allData.responseTime = dv.responseTime;
        allData.targMaxDiff = minMax.targMax - minMax.targMin;
        allData.distMaxDiff = minMax.distMax - minMax.distMin;
        allData.targRMS = targRMS;
        
        allData.cMap = [repmat(partCode,length(classIm.muContrast(:)),1) ...
            (1:length(classIm.muContrast(:)))'...
            classIm.muContrast(:)];
        
        allData.total = [repmat(partCode,length(classIm.total(:)),1) ...
            (1:length(classIm.total(:)))'...
            classIm.total(:)];

        allData.filtEn = [repmat(partCode,length(classIm.filtEn(:)),1) ...
            (1:length(classIm.filtEn(:)))'...
            classIm.filtEn(:)];
        
        allData.model.accuracy = dv.model.accuracy;
    else
           
        allData.participant = [allData.participant;repmat(partCode ,trial,1)];
        allData.noiseLevel = [allData.noiseLevel;task.noiseLevel];
        allData.nGabors = [allData.nGabors;task.nGabors];
        allData.sameDiff = [allData.sameDiff;task.sameDiff];
        allData.targPos = [allData.targPos;task.targPos];
        allData.allCats = [allData.allCats;stimuli.allCats];
        allData.imRef = [allData.imRef;stimuli.imRef];
        allData.similarity = [allData.similarity; similarity];
        allData.accuracy = [allData.accuracy;dv.accuracy];
        allData.chosenSquare = [allData.chosenSquare;dv.chosenSquare];
        allData.responseTime = [allData.responseTime;dv.responseTime];
        allData.targMaxDiff = [allData.targMaxDiff; minMax.targMax - minMax.targMin];
        allData.distMaxDiff = [allData.distMaxDiff; minMax.distMax - minMax.distMin];
        allData.targRMS = [allData.targRMS; targRMS];
        
        allData.cMap = [allData.cMap; repmat(partCode,length(classIm.muContrast(:)),1) ...
            (1:length(classIm.muContrast(:)))'...
            classIm.muContrast(:)];
        
        allData.total = [allData.total; repmat(partCode,length(classIm.total(:)),1) ...
            (1:length(classIm.total(:)))'...
            classIm.total(:)];

        allData.filtEn = [allData.filtEn; repmat(partCode,length(classIm.filtEn(:)),1) ...
            (1:length(classIm.filtEn(:)))'...
            classIm.filtEn(:)];
        
        allData.model.accuracy = [allData.model.accuracy;dv.model.accuracy];
    end
    
end

cd(md)