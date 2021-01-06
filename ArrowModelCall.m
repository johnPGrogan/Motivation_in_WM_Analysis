function output = ArrowModelCall(dataStruct, model, splitBy, useMemFit)
% output = ArrowModelCall(dataStruct, model, splitBy, useMemFit)
% take in formatted modell data
% fit model to each pp, splitting by splitBy

%%
if ~exist('model','var')
    model = SwapModel(); % default model to use
end
if ~exist('splitBy','var')
    splitBy = 'trialTypes';%run each trial type separately - default
end
if ~exist('useMemFit','var')
    useMemFit = 0;%use MLE by default
end


nPP = length(dataStruct); % number of pps/sessions
nPars = length(model.upperbound); % number of parameters

trialTypesVec = dataStruct(1).(splitBy);%trial types
trialTypes = unique(trialTypesVec);%num uniques
nTrialTypes = length(trialTypes); % number of trial types to split data by

nTrialsEachType = NaN(1,1,nTrialTypes);
for i = 1:nTrialTypes%num of each trial type
    nTrialsEachType(1,1,i) = sum(trialTypesVec == trialTypes(i));
end


%% fitting
fits = cell(nPP,1,nTrialTypes); % initialise
pars = NaN(nPP,nPars,nTrialTypes);
negLogLike = NaN(nPP,1,nTrialTypes);

fprintf('Fitting model: %s',model.name)
for iPP = 1:nPP
    thisSession = dataStruct(iPP);
    fprintf('\npp: %d/%d:     ',iPP,nPP)
    
    % get splitbby field
    thisSession.data.(splitBy) = thisSession.(splitBy);
    
    %  split by that
    [datasets, conditionOrder] = SplitDataByField(thisSession.data,splitBy);
    nTrialTypes = length(conditionOrder);
    
    for iTrialType = 1:nTrialTypes%for each trial type
        fprintf('\b\b\b%d/%d',iTrialType, nTrialTypes)
        
        if useMemFit
            %fit using memtoolbox - gives logLikelihood
            fits{iPP,1,iTrialType} = MemFit(datasets(iTrialType),model,'Verbosity',0);
            pars(iPP,:,iTrialType) = fits{iPP,iTrialType}.maxPosterior;%max posteriors
            negLogLike(iPP,1,iTrialType) = -max(fits{iPP,iTrialType}.posteriorSamples.like);%neg max like
            
        else
            
            [pars(iPP,:,iTrialType), logLike] = MLE(datasets{iTrialType},model);
            negLogLike(iPP,1,iTrialType) = -logLike;%make negative
        end
    end
    
end
fprintf('\n')

%% 

%store in output
output.fits = fits;
output.pars = pars;
output.negLogLike = negLogLike;
%calculate AIC and BIC
aic = 2 * (negLogLike + nPars);
bic = 2 * (negLogLike) + nPars .* log(repmat(nTrialsEachType,nPP,1,1));
output.aic = aic;
output.bic = bic;
output.nPars = nPars;
output.splitBy = splitBy;
output.nTrialTypes = nTrialTypes;
output.nTrialsEachType = nTrialsEachType;

end%ArrowModelFits
