% ArrowDemographics
% Get ppID from result files, match to ppID in xlsx, get age, gender from
% that

clc;clear;

% demoFilePath = '../../../../OneDrive - Nexus365/Comp_Neurology_Lab/John/Experiments/Forms/WM/Pp_Info_WM_Motiv.xlsx';
demoFilePath = '../../Experiments/Forms/WM/Pp_Info_WM_Motiv.xlsx';
%%

% get ppIDs for each experiment
matFiles = {'ArrowMotivAnalysis.mat';
            'ArrowMotivCueBlockedAnalysis.mat';
            'ArrowMotivRetrocueAnalysis.mat';
            'ArrowMotivReallocateAnalysis.mat';
            'ArrowMotivMotorAnalysis_A.mat';
            'ArrowMotivMotorAnalysis_B.mat';
            };
        
nExpt = length(matFiles);
dataFolder = './Data/';

[exptFiles, exptNames, ppID] = deal(cell(nExpt,1));
for i = 1:nExpt
    
    d = load(fullfile(dataFolder, matFiles{i}), 'dataFiles');
    
    exptFiles{i} = d.dataFiles;
    
    exptNames{i} = cellfun(@(x) x(regexp(x,'Cue\\')+4:regexp(x,'_')-1), exptFiles{i},'UniformOutput',0);
    ppID{i} = cellfun(@(x) x(regexp(x,'_')+1:end-4), exptFiles{i},'UniformOutput',0);
          
end

  
%% match

t = readtable(demoFilePath);

v = t.Properties.VariableNames;

% get info
for i = 1:nExpt
    
    [~,xlsInds{i}] = ismember(ppID{i},t.PtID); % get IDs 
    
    ages{i} = t.Age(xlsInds{i});
    sex{i} = t.Sex(xlsInds{i});
    
end

ages = array2table(sq(nancat(ages)));
sex = array2table(sq(nancat(sex)));

%% summary table

v = ages.Properties.VariableNames;
experiment = {'1'; '2'; '3'; '4'; '5a'; '5b'};

sumTab = table(experiment, cell(nExpt,1), cell(nExpt,1), cell(nExpt,1));
sumTab.Properties.VariableNames = {'Experiment', 'N','Ages','Gender'};
for i = 1:nExpt
    
    sumTab.N{i} = sum(~isnan(ages.(v{i})));
    sumTab.Ages{i} = sprintf('%.2f (%.2f)', nanmean(ages.(v{i})), nanstd(ages.(v{i})));
    sumTab.Gender{i} = sprintf('%d:%d', sum(strcmp(sex.(v{i}),'M')), sum(strcmp(sex.(v{i}),'F')));
    
end

disp(sumTab)