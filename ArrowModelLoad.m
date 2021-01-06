function [dataStruct] = ArrowModelLoad(dataStruct, removeRetrocueIncorr)
% [dataStruct] = ArrowModelLoad(dataStruct)
% Load up data from each structure in dataStruct, format for memtoolbox fitting
% and convert to degrees.
% The new data is stored in dataStruct.data

%% load up beh data for each participant, convert to degrees and format for
% use in memtoolbox


nPP = length(dataStruct);%number of result files
% fprintf('Formatting data for model fits...\n')
for iPP = 1:nPP%for each result file
%     fprintf('%d/%d\n',iPP,nPP)
    
    thisSession = dataStruct(iPP); %get just this result
    
    nTrials = length(thisSession.trialTypes);%number of trials
    nArrows = thisSession.result.data(1).nTargets; % assumes same number of tar
    
    %get beh data
    targAngles = [thisSession.result.data.targetAngle];%target angles
    nonTargAngles = nancat(1, thisSession.result.data.angles)';%distractor angles
    nonTargAngles = nonTargAngles(2:end,:);%remove targets
    
    
    %get response angle
    responseAngles = [thisSession.result.data.response];
    
    %put into structure for memtoolbox
    %errors is response - target
    data = struct();
    data.errors = mod(responseAngles - targAngles + pi, pi*2) - pi;%-pi:pi
    for i = 1:size(nonTargAngles,1)%same for distractors
        %this is distance from target to distractors
        data.distractors(i,:) = mod(nonTargAngles(i,:) - targAngles + pi, pi*2) - pi;
    end
    
    %%%% add in previous response as a distractor
    prevRespAngles = [NaN,responseAngles(1:end-1)];   
    data.prevResp = rad2deg(mod(prevRespAngles - targAngles + pi, pi*2) - pi);
    
    
    data.resps = rad2deg(responseAngles);%store response angles
    %make pdf of resps
    data.respPdf = fitdist(data.resps','Kernel','Kernel','normal','Width',10);%fit dist of resps - 10deg bins
    
    
    data.errors = rad2deg(data.errors);%convert to degrees
    data.distractors = rad2deg(data.distractors);
    data.trialTypes = thisSession.trialTypes;
    
    data.targets = rad2deg(targAngles);

    if exist('removeRetrocueIncorr','var') && removeRetrocueIncorr
    	data.errors(dataStruct(iPP).RetrocueResponseCorrect==0) = NaN; 
    end

    dataStruct(iPP).data = data;%store in output
    

end


end%ArrowModelLoad