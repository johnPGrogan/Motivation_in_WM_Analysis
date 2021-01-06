function [allSwap, interpModelData] =  ArrowModelFitsTrajNormalise(dataFiles, dataStruct, nInterps, isBlocked)
% clear
% close all

traj = ArrowTrajAnalysis(dataFiles, dataStruct, 0, nInterps); % load up the trajectories

[nPP, nTr, nPoints] = size(traj.interpTraj); % this is the interpolated one


%% get model formatted data
modelStruct = ArrowModelLoad(dataStruct);%get standard data in model format

if exist('isBlocked','var') && isBlocked % if isblocked expt, collapse across catch trials
    %% for expt 2, combine across catch trials
    for i = 1:nPP % combine catchTrials and nonCatchTrials together as it shouldn't affect anything

        % make trialtypes 1:4 be prelow, postlow, prehi, posthi
        trialType = modelStruct(i).rewardLevel + modelStruct(i).cueTime; % 2=preLow, 3=postLow, 51=preHi, 52=postHi
        v = [2 3 51 52];
        for j = 1:4 % replace with 1:4
            trialType(trialType==v(j)) = j;
        end

        modelStruct(i).trialTypes = trialType;
    end
end

%% prepare to fit the model at each time point

for iT = 1:nPoints
    for iPP = 1:nPP
        resps = sq(traj.interpAngles(iPP,:,iT))';
        resps = mod(resps,2*pi);
        
        errors = mod(resps - traj.targAngles(iPP,:) + pi, 2*pi) - pi;
        
        distractors =  mod( permute(traj.nonTargAngles(iPP,:,:),[3,2,1]) - traj.targAngles(iPP,:) + pi, 2*pi) - pi;
        
        interpModelData(iT,iPP).data.errors = rad2deg(errors);
        interpModelData(iT,iPP).data.distractors = rad2deg(distractors);
        interpModelData(iT,iPP).data.trialTypes = modelStruct(iPP).trialTypes;
        interpModelData(iT,iPP).trialType = modelStruct(iPP).trialTypes;
        interpModelData(iT,iPP).allTrials = modelStruct(iPP).allTrials;
    end
end

%% run model fitting

splitBy = 'trialType';%trialType or allTrials
useMemFit = 0; % 0 = MLE, 1 = MemFit

parfor iT = 1:nPoints
    fprintf('time = %d/%d\n',iT,nPoints)
    allSwap(iT,1) = ArrowModelCall(interpModelData(iT,:), SwapModel(), splitBy,useMemFit);
%     allMix(iT,1) = ArrowModelCall(trajInterp(iT,:), StandardMixtureModel(), splitBy,useMemFit);
end

%% save

% save('./Data/ArrowModelFitsTrajNormalise.mat','allSwap','allMix','nPoints','nPP',...
%     'timePoints','interpCoords','interpAngles','initTime','trajInterp')