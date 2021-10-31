function output = ArrowDataLoad(dataFiles)

% load up result file
n = length(dataFiles);
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
for i = 1:n
    d = load(dataFiles{i});
    output1 = ExtractSummary(d.result);
    if i>1
        [output, output1] = ensureStructsAssignable(output, output1,false);
    end
    output(i) = output1;
end
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle')

end

function output = ExtractSummary(result)
% extract trial data and info, store in vectors, get conditions etc
% also split precision and RT by trial types

    %% expt info
    
    nTrials = length(result.data); % number of trials
    output.nTrials = nTrials;
    
    % store screen dims as complex co-ordinates
    output.screenDims = result.params.screenSize(1) + result.params.screenSize(2)*1j;
  
    output.blocks = [result.data.block];
    output.trialNums = [result.data.trialIndex];
    
    
    %% timings
    
    timeLabels = result.params.dataFields.times;
    
    % keep only used labels
    isLabel = cellfun(@(x) isfield(result.data(1), x), timeLabels);
    timeLabels = timeLabels(isLabel);
    
    times = NaN(length(timeLabels),nTrials);
    for iTime = 1:length(timeLabels)
        for iTrial = 1:nTrials
             times(iTime,iTrial) = result.data(iTrial).(timeLabels{iTime});
        end
    end  
    
    output.timeStamps =  times - repmat(times(1,:),size(times,1),1);%time since trial start
    output.timeDiffs = diff(times,1,2);%time since last timestamp
    
    % get duration of entire experiment
    if isfield(result, 'practiceResult') % if practice
        exptStart = result.practiceResult(1).trialStart; % start of first practice trial
        if isnan(exptStart) % if this was skipped, use first main trial
            exptStart = result.data(1).trialStart;
        end
    else
        exptStart = result.data(1).trialStart;
    end
    
    exptEnd = max(times,[],'all');
    
    output.exptDuration = exptEnd - exptStart;
    
    %% check which variable names are used
    
    if isfield(result.data(1), 'probeStart')
        probeStart = 'probeStart';
    elseif isfield(result.data(1), 'startProbe')
        probeStart = 'startProbe';
    else
        probeStart = '';
    end
    if isfield(result.data(1), 'responseStart')
        respStart = 'responseStart';
    elseif isfield(result.data(1), 'startResponse')
        respStart = 'startResponse';
    else
        respStart = '';
    end
    if isfield(result.data(1), 'responseEnd')
        respEnd = 'responseEnd';
    elseif isfield(result.data(1), 'endChoice')
        respEnd = 'endChoice';
    else
        respEnd = '';
    end
    
    % get RTs
    if ~isempty(respStart) && ~isempty(probeStart)
        output.startRTs = ([result.data.(respStart)] - [result.data.(probeStart)]) .*1000;%ms
    end
    if ~isempty(respEnd) && ~isempty(respStart)
        output.endRTs = ([result.data.(respEnd)] - [result.data.(respStart)]) .*1000;%ms
    end
    
    % retrocue
    if isfield(result.data(1), 'ResponseRetrocue')
        cueResp = 'ResponseRetrocue';
    else
        cueResp = '';
    end
    if isfield(result.data(1), 'ShowRetrocue')
        cueStart = 'ShowRetrocue';
    else
        cueStart = '';
    end
    if ~isempty(cueStart) && ~isempty(cueResp)
        output.cueRTs = ([result.data.(cueResp)] - [result.data.(cueStart)]) .*1000;%ms
    end
    
    % catch trial RT
    if isfield(result.data(1), 'CatchResponse') && isfield(result.data(1),'CatchShown')
        output.catchRTs = ([result.data.CatchResponse] - [result.data.CatchShown]) .* 1000;
    end
        
    
    
    %% get single variables
    
    singleVars = {'precision','reward','won','response','button',...
        'targetIndex','targetColour','catchCorrect','RetrocueResponseCorrect'}; % possible measures to get
    
    for iVar = 1:length(singleVars)
        if isfield(result.data(1), singleVars{iVar}) % extract if exists
            output.(singleVars{iVar}) = [result.data.(singleVars{iVar})];
        end
    end
    
    
    %% multivars - mutliple values per trial, e.g. stimuli info
    
    multiVars = {'angles','colourOrder'};
    
    for iVar = 1:length(multiVars)
        if isfield(result.data(1), multiVars{iVar})
            output.(multiVars{iVar}) = reshape([result.data.(multiVars{iVar})],[],nTrials);
        end
    end
        
    output.loci = [result.data.locations];%stim locations
    output.stimLoci = output.loci(:,1:2:end) + output.loci(:,2:2:end)*1j;%store in complex
    
        
    %% trial info - conditions
       
    if isfield(result.params, 'blockVariables')
        trialVars = [fieldnames(result.params.blockVariables);
                     fieldnames(result.params.trialVariables)];
    else
        trialVars = fieldnames(result.params.trialVariables);
    end
    
    output.trialVals = NaN(length(trialVars), nTrials);
    for iVar = 1:length(trialVars)
        if isfield(result.data(1), trialVars{iVar})
            output.(trialVars{iVar}) = [result.data.(trialVars{iVar})];
            output.trialVals(iVar,:) = output.(trialVars{iVar});
        end
    end
    
    output.trialVarNames = trialVars;
    
    %% sort by trial type
    
    output.uniqueTrialVals = unique(output.trialVals','rows'); % get unique trials
    output.nTrialTypes = size(output.uniqueTrialVals,1);
    
    % turn it into vector
    for i = 1:output.nTrialTypes 
        trialTypes(:,i) = all(output.trialVals' == output.uniqueTrialVals(i,:),2);
    end
    output.trialTypes = (trialTypes * [1:output.nTrialTypes]')';
    
    % get the trial vars for the ones that vary
    levelInds = structfun(@length, result.params.trialVariables) > 1;
    output.colTrialTypes = output.uniqueTrialVals(:,levelInds)';
    
    %% split some vars by trial type
    
    varsToSplit = {'precision','prec';
                   'startRTs','startRT';
                   'endRTs','endRT';
                   'cueRTs','cueRT';
                   'trialTypes','tt';
                   'catchCorrect','catchCorr';
                   'catchRTs','catchRT';
                   'RetrocueResponseCorrect','retrocueCorr'};
     
    for i = 1:size(varsToSplit,1)
        if isfield(output, varsToSplit{i,1})
            output.(varsToSplit{i,2}) = groupMeans(output.(varsToSplit{i,1}), 2, output.trialTypes,'dim');
        end
    end

    %% 
    
    output.allTrials = ones(1,nTrials);%1 for all 
    output.result = result;
end
