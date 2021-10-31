% ArrowModelComparison

%% fit 3 models to all expts, compare bic

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files = what(dataFolder);%find files there


isOutput = regexp(files.mat,'ArrowMotiv');%keep only files with this in
isOutput = ~cellfun(@isempty,isOutput);
files.mat = files.mat(isOutput);

dataInds = 1:length(files.mat);%change this to only analyse certain files
dataFiles = fullfile(dataFolder, files.mat(dataInds));%make path to files

dataStruct = ArrowDataLoad(dataFiles);%analyse each file
n = length(dataStruct);%num pps

%% get to model

modelStruct = ArrowModelLoad(dataStruct); % get data for model fits

for i = 1:n
    modelStruct(i).data.prevResp = ones(size(modelStruct(i).data.errors))*180;%mod(modelStruct(i).data.targets + 360,360)-180;    
end

%% fit

models = {StandardMixtureModel(), SwapModel(), SwapModelPrev(), SwapModelNoGuess() };
nModels = length(models);
modelNames = {'standard', 'misbinding','swap180'};

splitBy = 'allTrials'; % fit each condition separately
useMemFit = 0; % use MLE non MemFit

fits = cell(nModels,1);
for i = 1:nModels
    fits{i} = ArrowModelCall(modelStruct, models{i}, splitBy, useMemFit);
end

fits = nancat(fits);
%% BIC

bic = nancat(2, fits.bic);

meanBIC = nanmean(bic);
disp(meanBIC);

[m,i] = min(meanBIC);

fprintf('\nbest model is %s\n', modelNames{i});

disp(meanBIC - m);


%% remove motor?

isMotor = cellRegexpi(files.mat, 'Motor') ~= 0;

bic2 = bic(~isMotor,:);

meanBIC2 = nanmean(bic2);
disp(meanBIC2);

%% by expt

exptNames = {'MotivCue_', 'CueBlocked', 'Retrocue', 'Reallocate', 'MotorControl'};
nExpts = length(exptNames);

for i = 1:nExpts
    exptNum(:,i) = cellRegexpi(files.mat, exptNames{i}) ~= 0;
end
exptNum = exptNum * [1:nExpts]';


bicByExpt = groupMeans(bic, 1, exptNum,'dim');

meanBICByExpt = nanmean(bicByExpt, 3)';
disp(meanBICByExpt - min(meanBICByExpt,[],2));