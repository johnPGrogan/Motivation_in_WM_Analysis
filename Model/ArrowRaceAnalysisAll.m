% ArrowRaceAnalysisAll

clear
close all

%% defaults
drawKs = 0;
drawPars = 1;
rtCols = [3];
drawRTs = 0;
%% 

a = load('../Data/ArrowMotivCueRace.mat','output');

motivFit = a.output;

motivFit.condCols = {[1 2 3 4; 5 6 7 8]; [1 2 5 6;3 4 7 8]; [1 3 5 7; 2 4 6 8]};
motivFit.rtLabels = {'startRT','endRT','totalRT'}; %name of rt(s) to draw
motivFit.rtCols = rtCols; % which ones to draw/analyse
% conds = {'lowPreShort','lowPreLong','lowPostShort','lowPostLong','hiPreShort','hiPreLong','hiPostShort','hiPostLong'}
% factors for anova
% rew = permute(repmat([0 0 0 0 1 1 1 1],30,1,1),[1,3,2]);
% cue = permute(repmat([0 0 1 1 0 0 1 1],30,1,1),[1,3,2]);
% del = permute(repmat([0 1 0 1 0 1 0 1],30,1,1),[1,3,2]);
% pp  = permute(repmat(1:30, 1,1,8), [2,1,3]);
% motivFit.modelMat = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0;
%     1 1 0 0; 1 0 1 0; 0 1 1 0; 1 1 1 0];
% motivFit.factorNames = {'rew','cue','del','pp'};
% motivFit.factors = [col(rew),col(cue),col(del),col(pp)];

% for rmanova
motivFit.factorNames = {'pp','delay','cue','rew'};
motivFit.reshapeDims = [motivFit.nPP 2 2 2]; % how to reshape each parameter 

motivFit.condInds = {[1 5 3 7]; [2 6 4 8]};
motivFit.xTickLabels = {'Low','High','Low','High'}; % for figure
motivFit.xLabel = 'Pre                Post';
motivFit.legendLabels = {'Short','Long'};
motivFit.factorOrder = [3 2 1]; % order of factors to use for * on figures

motivFit.drawKs = drawKs;
motivFit.drawPars = drawPars;
motivFit.drawRTs = drawRTs;

motivFit = ArrowRaceAnalysis(motivFit); % analyse + draw figs

%%

a = load('../Data/ArrowMotivCueBlockedRace.mat','output');

blockedFit = a.output;

blockedFit.condCols = {[1 3; 2 4]; [1 2; 3 4]}; % split cond columns by these 
blockedFit.rtLabels = {'startRT','endRT','totalRT'}; %name of rt(s) to draw
blockedFit.rtCols = rtCols; % which ones to draw/analyse
% conds = {'lowPre','highPre','lowPost','highPost'};

% factors for anova
% rew = permute(repmat([0 0 1 1],blockedFit.nPP,1,1),[1,3,2]);
% cue = permute(repmat([0 1 0 1],blockedFit.nPP,1,1),[1,3,2]);
% pp  = permute(repmat(1:blockedFit.nPP, 1,1,blockedFit.nConds), [2,1,3]);
% blockedFit.factors = {col(rew), col(cue), col(pp)};
% blockedFit.factorNames = {'reward','cueTime','pp'};
% blockedFit.modelMat = [0 0 1; 1 0 0; 0 1 0; 1 1 0;];

% for rmanova
blockedFit.factorNames = {'pp','rew','cue'};
blockedFit.reshapeDims = [blockedFit.nPP 2 2]; % how to reshape each parameter 


blockedFit.condInds = {[1 3]; [2 4]}; % how to split conditions for figure
blockedFit.xTickLabels = {'Low','High'}; % for figure
blockedFit.xLabel = 'Reward';
blockedFit.legendLabels = {'Pre','Post'};
blockedFit.factorOrder = [1 2]; % order of factors to use for * on figures

blockedFit.drawKs = drawKs;
blockedFit.drawPars = drawPars;
blockedFit.drawRTs = drawRTs;

blockedFit = ArrowRaceAnalysis(blockedFit); % analyse + draw figs


%%

a = load('../Data/ArrowMotivReallocateRace.mat','output');

reallocateFit = a.output;

% rearrange the pars here, so that 1st dim in rmanova is reward
fn = {'RT', 'Corr', 'pars', 'data'};
for i = 1:length(fn)
    reallocateFit.(fn{i}) = reallocateFit.(fn{i})(:,[1 4 5 8 2 3 6 7],:);
end
% change into:
% conds = {'motivPreLow', 'motivPreHigh', 'motivPostLow', 'motivPostHigh', 'realPreLow', 'realPreHigh', 'realPostLow', 'realPostHigh'};


% reallocateFit.condCols = {[1 2 3 4; 5 6 7 8]; [1 2 5 6;3 4 7 8];  [1 4 5 8; 2 3 6 7]}; % split cond columns by these 
reallocateFit.condCols = {[1 3 5 7; 2 4 6 8]; [1 2 5 6; 3 4 7 8]; [1 2 3 4; 5 6 7 8]};
reallocateFit.rtLabels = {'startRT','endRT','totalRT'}; %name of rt(s) to draw
reallocateFit.rtCols = rtCols; % which ones to draw/analyse

% cue = permute(repmat([0 0 0 0 1 1 1 1],reallocateFit.nPP,1,1),[1,3,2]);
% rew = permute(repmat([0 0 1 1 0 0 1 1],reallocateFit.nPP,1,1),[1,3,2]);
% motiv = permute(repmat([1 0 0 1 1 0 0 1],reallocateFit.nPP,1,1),[1,3,2]);
% pp  = permute(repmat(1:reallocateFit.nPP, 1,1,reallocateFit.nConds), [2,1,3]);
% reallocateFit.modelMat = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0;
%     1 1 0 0; 1 0 1 0; 0 1 1 0; 1 1 1 0];
% reallocateFit.factorNames = {'rew','cue','motiv','pp'};
% reallocateFit.factors = {col(rew),col(cue),col(motiv),col(pp)};

% for rmanova
reallocateFit.factorNames = {'pp','rew','cue','motiv'};
reallocateFit.reshapeDims = [reallocateFit.nPP 2 2 2]; % how to reshape each parameter 

reallocateFit.condInds = {[1 2 5 6]; [3 4 7 8]};
reallocateFit.xTickLabels = {'Low','High','Low','High'};
reallocateFit.xLabel = 'Motivate          Reallocate';
reallocateFit.legendLabels = {'Pre','Post'};
reallocateFit.factorOrder = [1 2 3];

reallocateFit.drawKs = drawKs;
reallocateFit.drawPars = drawPars;
reallocateFit.drawRTs = drawRTs;

reallocateFit.splitAnova = 1; % also split into motiv and realloc for sep anova

reallocateFit = ArrowRaceAnalysis(reallocateFit);


%%

a = load('../Data/ArrowMotivRetrocueRace.mat','output');

retrocueFit = a.output;

retrocueFit.condCols =  {[1 2; 3 4]; [1 3; 2 4]};
retrocueFit.rtLabels = {'startRT','endRT','totalRT'}; %name of rt(s) to draw
retrocueFit.rtCols = rtCols; % which ones to draw/analyse
% conds = {'lowIncon','lowCon','highIncon','highCon'};

% rew = permute(repmat([0 0 1 1],retrocueFit.nPP,1,1),[1,3,2]);
% cue = permute(repmat([0 1 0 1],retrocueFit.nPP,1,1),[1,3,2]);
% pp  = permute(repmat(1:retrocueFit.nPP, 1,1,retrocueFit.nConds), [2,1,3]);
% 
% retrocueFit.modelMat = [0 0 1; 1 0 0; 0 1 0;
%     1 1 0;];
% retrocueFit.factorNames = {'rew','cue','pp'};
% retrocueFit.factors = {col(rew),col(cue),col(pp)};

% for rmanova
retrocueFit.factorNames = {'pp','congr','rew'};
retrocueFit.reshapeDims = [retrocueFit.nPP 2 2]; % how to reshape each parameter 

retrocueFit.condInds = {[1 3]; [2 4]};
retrocueFit.xTickLabels = {'Low','High'};
retrocueFit.xLabel = 'Reward';
retrocueFit.legendLabels = {'Incongruent','Congruent'};
retrocueFit.factorOrder = [2 1];

retrocueFit.drawKs = drawKs;
retrocueFit.drawPars = drawPars;
retrocueFit.drawRTs = drawRTs;

retrocueFit = ArrowRaceAnalysis(retrocueFit);

%% motor 1
% 
% a = load('../Data/ArrowMotivMotorARace.mat','output');
% 
% motorAFit = a.output;
% 
% motorAFit.condCols =  {[1; 2]};
% motorAFit.rtLabels = {'startRT','endRT','totalRT'}; %name of rt(s) to draw
% motorAFit.rtCols = rtCols; % which ones to draw/analyse
% % conds = {'lo','hi'};
% % rew = permute(repmat([0 1],motor1Fit.nPP,1,1),[1,3,2]);
% % pp  = permute(repmat(1:motor1Fit.nPP, 1,1,motor1Fit.nConds), [2,1,3]);
% % 
% % motor1Fit.modelMat = [0 1; 1 0;];
% % motor1Fit.factorNames = {'rew','pp'};
% % motor1Fit.factors = {col(rew), col(pp)};
% 
% % for rmanova
% motorAFit.factorNames = {'pp','rew'};
% motorAFit.reshapeDims = [motorAFit.nPP 2]; % how to reshape each parameter 
% 
% motorAFit.condInds = {[1]; [2]};
% motorAFit.xTickLabels = {'Low','High'};
% motorAFit.xLabel = 'Reward';
% % motor1Fit.legendLabels = {'Low', 'High'};
% motorAFit.factorOrder = [1];
% 
% motorAFit.drawKs = drawKs;
% motorAFit.drawPars = drawPars;
% motorAFit.drawRTs = drawRTs;
% 
% motorAFit = ArrowRaceAnalysis(motorAFit);
% 
% %% motor 2
% 
% a = load('../Data/ArrowMotivMotorBRace.mat','output');
% 
% motorBFit = a.output;
% 
% motorBFit.condCols =  {[1; 2]};
% motorBFit.rtLabels = {'startRT','endRT','totalRT'}; %name of rt(s) to draw
% motorBFit.rtCols = rtCols; % which ones to draw/analyse
% % conds = {'lo','hi'};
% 
% % rew = permute(repmat([0 1],motor2Fit.nPP,1,1),[1,3,2]);
% % pp  = permute(repmat(1:motor2Fit.nPP, 1,1,motor2Fit.nConds), [2,1,3]);
% % 
% % motor2Fit.modelMat = [0 1; 1 0;];
% % motor2Fit.factorNames = {'rew','pp'};
% % motor2Fit.factors = {col(rew), col(pp)};
% 
% % for rmanova
% motorBFit.factorNames = {'pp','rew'};
% motorBFit.reshapeDims = [motorBFit.nPP 2]; % how to reshape each parameter 
% 
% motorBFit.condInds = {[1]; [2]};
% motorBFit.xTickLabels = {'Low','High'};
% motorBFit.xLabel = 'Reward';
% % motor2Fit.legendLabels = {'Low', 'High'};
% motorBFit.factorOrder = [1];
% 
% motorBFit.drawKs = drawKs;
% motorBFit.drawPars = drawPars;
% motorBFit.drawRTs = drawRTs;
% 
% motorBFit = ArrowRaceAnalysis(motorBFit);


%%

save('../Data/ArrowRaceAnalysisAll.mat','motivFit','blockedFit','reallocateFit','retrocueFit')%, 'motorAFit','motorBFit')