% ArrowRaceModelFittingAll
% Run race model fitting on each experiment (at end of script)

clear
close all
clc

%% MotivCue

% disp('ArrowMotivCue')
load('../Data/ArrowMotivAnalysis.mat','allStartRT','allEndRT','allRT');
load('../Data/ArrowMotivModelAnalysis.mat','modelStruct','allSwap');

% fit model - collapsing across participants
allRT = nancat(4, allStartRT, allEndRT, allRT);

expt{1}.rt = allRT;
expt{1}.modelStruct = modelStruct;
expt{1}.pars = allSwap.pars;
expt{1}.splitBy = 'trialTypes';
expt{1}.saveName = 'ArrowMotivCueRace';

% output = ArrowRaceModelFitting(allRT, dataStruct, allSwap.pars, 'trialTypes');

% save('../Data/ArrowMotivCueRace.mat','output')
%% MotivCueBlocked
% clear
% disp('ArrowMotivCueBlocked')
load('../Data/ArrowMotivCueBlockedAnalysis.mat','allStartRT','allEndRT','allTotalRT');
load('../Data/ArrowMotivCueBlockedModelling.mat','modelStruct','allSwap');

% fit model - collapsing across participants
allRT = nancat(4, allStartRT, allEndRT, allTotalRT);

expt{2}.rt = allRT;
expt{2}.modelStruct = modelStruct;
expt{2}.pars = allSwap.pars;
expt{2}.splitBy = 'trialTypes2';
expt{2}.saveName = 'ArrowMotivCueBlockedRace';
% output = ArrowRaceModelFitting(allRT, modelStruct, allSwap.pars, 'trialTypes2');

% save('../Data/ArrowMotivCueBlockedRace.mat','output')

%% MotivReallocate
% clear
% disp('ArrowMotivReallocate')
load('../Data/ArrowMotivReallocateAnalysis.mat','allStartRT','allEndRT','allRT');
load('../Data/ArrowMotivReallocateModelling.mat','modelStruct','allSwap');

% fit model - collapsing across participants
allRT = nancat(4, allStartRT, allEndRT, allRT);

expt{3}.rt = allRT;
expt{3}.modelStruct = modelStruct;
expt{3}.pars = allSwap.pars;
expt{3}.splitBy = 'trialTypes';
expt{3}.saveName = 'ArrowMotivReallocateRace';

% output = ArrowRaceModelFitting(allRT, modelStruct, allSwap.pars);

% save('../Data/ArrowMotivReallocateRace.mat','output')


%% MotivRetrocue
% 
% clear
% disp('ArrowMotivRetrocue')
load('../Data/ArrowMotivRetrocueAnalysis.mat','allStartRT','allEndRT','allRT');
load('../Data/ArrowMotivRetrocueModelling.mat','modelStruct','allSwap');

% fit model - collapsing across participants
allRT = nancat(4, allStartRT, allEndRT, allRT);

expt{4}.rt = allRT;
expt{4}.modelStruct = modelStruct;
expt{4}.pars = allSwap.pars;
expt{4}.splitBy = 'trialTypes';
expt{4}.saveName = 'ArrowMotivRetrocueRace';

% output = ArrowRaceModelFitting(allRT, modelStruct, allSwap.pars);

% save('../Data/ArrowMotivRetrocueRace.mat','output')
% don't do this, as too few incorrect responses
% %% MotivMotor 1
% 
% % clear
% % disp('ArrowMotivRetrocue')
% load('../Data/ArrowMotivMotorAnalysis_A.mat','allRewMeasures');
% load('../Data/ArrowMotivMotorModelling_A.mat','modelStruct','allSwap');
% 
% % fit model - collapsing across participants
% allRT = allRewMeasures(:,:,:,2:end);
% 
% expt{5}.rt = allRT;
% expt{5}.modelStruct = modelStruct;
% expt{5}.pars = allSwap.pars;
% expt{5}.splitBy = 'trialTypes';
% expt{5}.saveName = 'ArrowMotivMotorARace';
% 
% %% MotivMotor 2
% 
% % clear
% % disp('ArrowMotivRetrocue')
% load('../Data/ArrowMotivMotorAnalysis_B.mat','allRewMeasures');
% load('../Data/ArrowMotivMotorModelling_B.mat','modelStruct','allSwap');
% 
% % fit model - collapsing across participants
% allRT = allRewMeasures(:,:,:,2:end);
% 
% expt{6}.rt = allRT;
% expt{6}.modelStruct = modelStruct;
% expt{6}.pars = allSwap.pars;
% expt{6}.splitBy = 'trialTypes';
% expt{6}.saveName = 'ArrowMotivMotorBRace';
%%
parfor i = 1:4
    outputAll(i) = ArrowRaceModelFitting(expt{i});
end

%%

for i = 1:4
    output = outputAll(i);
    save(['../Data/' expt{i}.saveName], 'output');
end

%%

save('../Data/ArrowRaceModelFittingAll.mat','outputAll','expt')