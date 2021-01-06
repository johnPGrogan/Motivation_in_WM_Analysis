function o = ArrowTrajAnalysis(dataFiles, dataStruct, doPlots, nInterps)
% Plot trajectories on some trials, then polar hist init/final/all
% angles/errors
% 

% clear
% close all
% 
% %load datafiles paths
% load('./Data/ArrowMotivAnalysis.mat','dataFiles','dataStruct')

%%
%load traj data
traj = ArrowTrajLoad(dataFiles);

%%
nStim = length(traj(1).trials(1).stimAngles); % number of stimuli
nTrialsPPs = [dataStruct.nTrials]; % number of trials per pp
nTrials = max(nTrialsPPs);
nPP = length(dataStruct);



%% check precision is the same
[finalPrec, finalAngles, targAngles, initAngles, initPrec,...
    nonTargAngles, prec, resps, targs, nonTargs] = deal(NaN(nPP, nTrials));

for iPP=1:nPP
    
    %from traj data
    %final angle precision
    finalPrec(iPP,1:nTrialsPPs(iPP)) = [traj(iPP).trials.finalPrec]; % get prec from traj
    finalAngles(iPP,1:nTrialsPPs(iPP)) = [traj(iPP).trials.finalAngle]; % final angle   
    
    trajStims = reshape([traj(iPP).trials.stimAngles],nStim,[]); % and targ angle from traj
    targAngles(iPP,1:nTrialsPPs(iPP)) = trajStims(1,:);
   
    initAngles(iPP,1:nTrialsPPs(iPP)) = [traj(iPP).trials.initAngle]; % init angle
    initPrec(iPP,1:nTrialsPPs(iPP)) = [traj(iPP).trials.initPrec]; % init angle precision
    
    
    %data
    prec(iPP,1:nTrialsPPs(iPP)) = [dataStruct(iPP).precision]; % and from dataanalysis
    resps(iPP,1:nTrialsPPs(iPP)) = [dataStruct(iPP).response]; % response angles
    stims = reshape([dataStruct(iPP).angles],nStim,[]); % stim angles
    targs(iPP,1:nTrialsPPs(iPP)) = stims(1,:); % target angles
    
    
    %get nontarg angles
    for j = 1:nStim-1
        nonTargAngles(iPP,1:nTrialsPPs(iPP),j) = trajStims(j+1,:);
        nonTargs(iPP,1:nTrialsPPs(iPP),j) = stims(j+1,:);
    end
    
end

%% calc distances

newPrec = mod( resps - targs + pi, pi*2) - pi;%from data


% nonTargDist

finalNonTargDists = mod(finalAngles - nonTargAngles + pi, pi*2) - pi;
initNonTargDists = mod(initAngles - nonTargAngles + pi, pi*2) - pi;

% prev resp
prevRespAngles = [NaN(nPP,1),finalAngles(:,1:end-1)];

finalPrevRespDists = mod(finalAngles - prevRespAngles + pi, pi*2) - pi;
initPrevRespDists = mod(initAngles - prevRespAngles + pi, pi*2) - pi;
%% check traj data against other data

all(round(finalPrec,5) == round(prec,5) | isnan(finalPrec),'all')

all(newPrec == prec| isnan(finalPrec),'all')

all(finalPrec == finalPrec| isnan(finalPrec),'all')

all(abs(finalPrec) == abs(finalPrec)| isnan(finalPrec),'all')


%% look at all angles within traj

[allAngles, allTrajTimes, allTraj] = deal([]);
% fprintf('\n     ')
for iPP = 1:nPP
%     fprintf('\b\b\b\b\b%02d/%02d',iPP,nPP)
    allAngsPp = [];
    %for each person
    %for each trial
    for  iTrial = 1:nTrialsPPs(iPP)
        %get angles
        allAngsTr = reshape([traj(iPP).trials(iTrial).angles],1,1,[]);
        allAngsPp = nancat(2, allAngsPp, allAngsTr);
        
    end
    allAngles = nancat(1, allAngles, allAngsPp);    
    
    allTrajTimes = nancat(3, allTrajTimes, nancat(2, traj(iPP).trials.trajTime));
    
    allTraj = nancat(3, allTraj, nancat(2, traj(iPP).trials.traj));
end

%mod
allAngles = mod(allAngles,2*pi);
allTrajTimes = permute(allTrajTimes, [3,2,1]);
allTraj = permute(allTraj, [3,2,1]);
%%

allTargDists = mod(allAngles - repmat(targAngles,1,1,size(allAngles,3)) + pi,2*pi) - pi;
allPrevRespDists = mod(allAngles - repmat(prevRespAngles,1,1,size(allAngles,3)) + pi, 2*pi) - pi;

for i = 1:nStim-1
    allNonTargDists(:,:,:,i) = mod(allAngles - repmat(nonTargAngles(:,:,i),1,1,size(allAngles,3)) + pi, 2*pi) - pi;
end



%% look at responses in relation to stimulus locations

stimLoci = NaN(nStim, nTrials, nPP);
for iPP = 1:nPP
    
    stimLoci(:,1:nTrialsPPs(iPP),iPP) = reshape([traj(iPP).trials.stimLocCompl],nStim,[]); % and targ angle from traj
    
end

stimLoci = permute(stimLoci,[3,2,1]);
stimLoci = stimLoci ./ [960 + 540*1j];

stimLociAngles = mod(angle(stimLoci),2*pi);

%%
finalStimDists = mod(finalAngles - stimLociAngles + pi,2*pi) - pi;
initStimDists = mod(initAngles - stimLociAngles + pi,2*pi) - pi;

for i = 1:nStim
    allStimDists(:,:,:,i) = mod(allAngles - repmat(stimLociAngles(:,:,i),1,1,size(allAngles,3)) + pi, 2*pi) - pi;
end


if exist('doPlots','var') && doPlots
    drawPlots();
end

%% interpolate trajectories within movements (ignore InitRT)

if ~exist('nInterps','var') || isempty(nInterps)
    nInterps = 100;
end

allTrajTimesRel = allTrajTimes - allTrajTimes(:,:,1);
endTimes = max(allTrajTimesRel,[],3);


[interpTimePoints, interpTargDists, interpTraj, interpAngles] = deal(NaN(nPP,nTrials,nInterps));
% initTime = zeros(nPP,nTrials);
% trajTime = cell(nPP,nTrials);


for iPP = 1:nPP
    for iTrial = 1:nTrials
        
        interpTimePoints(iPP,iTrial,:) = linspace(0, endTimes(iPP,iTrial), nInterps); % points to interpolate along
        
        nPts = min([find(isnan(allTargDists(iPP,iTrial,:)),1,'first')-1, size(allTargDists,3)]); % number of non-nan points
        if nPts > 0
            interpTargDists(iPP,iTrial,:) = interp1(sq(allTrajTimesRel(iPP,iTrial,1:nPts)), sq(allTargDists(iPP,iTrial,1:nPts)), interpTimePoints(iPP,iTrial,:));
            interpTraj(iPP,iTrial,:) = interp1(sq(allTrajTimesRel(iPP,iTrial,1:nPts)), sq(allTraj(iPP,iTrial,1:nPts)), interpTimePoints(iPP,iTrial,:));
            interpAngles(iPP,iTrial,:) = interp1(sq(allTrajTimesRel(iPP,iTrial,1:nPts)), sq(allAngles(iPP,iTrial,1:nPts)), interpTimePoints(iPP,iTrial,:));

        end
%         trajTime{iPP,iTrial} = traj(iPP).trials(iTrial).trajTime; % get time of each mouse sample relative to first
%         initTime(iPP,iTrial) = trajTime{iPP,iTrial}(1);
        
%         timePoints(iPP,iTrial,:) = linspace(0,max(trajTime{iPP,iTrial}),nPoints); % get time points to sample
        
%         interpCoords(iPP,iTrial,:) = interp1(trajTime{iPP,iTrial}, traj(iPP).trials(iTrial).traj, timePoints,'linear'); % angle
        
%         interpAngles(iPP,iTrial,:) = interp1(trajTime{iPP,iTrial}, traj(iPP).trials(iTrial).angles,timePoints,'linear'); % angle from centre
    end
end


%%

o = workspace2struct();



function drawPlots()


%% plot some sample trials

figure()
for iPP = 1:9
    subplot(3,3,iPP)
    ArrowTrajPlotTrial(traj(iPP), 10,0)
end


%% plot
nBins = 100;
figure()

subplot(3,4,1)
polarhistogram(initAngles(:),nBins)
title('init angle')

subplot(3,4,2)
polarhistogram(initPrec(:),nBins)
title('init precision')

subplot(3,4,3)
polarhistogram(initNonTargDists(:),nBins)
title('init nTargDists')

subplot(3,4,4)
polarhistogram(initPrevRespDists(:),nBins)
title('init dist prev resp')

subplot(3,4,5)
polarhistogram(finalAngles(:),nBins)
title('final angle')

subplot(3,4,6)
polarhistogram(finalPrec(:),nBins)
title('final precision')

subplot(3,4,7)
polarhistogram(finalNonTargDists(:),nBins)
title('final nTargDists')

subplot(3,4,8)
polarhistogram(finalPrevRespDists(:),nBins)
title('final dist prev resp')


% %% hist

% figure()

subplot(3,4,9)
polarhistogram(allAngles(:),nBins)
title('all traj angles')

subplot(3,4,10)
polarhistogram(allTargDists(:),nBins)
title('all traj TargDists')

subplot(3,4,11)
polarhistogram(allNonTargDists(:),nBins)
title('all traj nonTargDists')

subplot(3,4,12)
polarhistogram(allPrevRespDists(:),nBins)
title('all traj prevRespDists')
%%

figure()

subplot(2,2,1)
polarhistogram(finalStimDists(:),nBins)
title('final angle dist to stim loci')

subplot(2,2,2)
polarhistogram(initStimDists(:),nBins)
title('initial angle dist to stim loci')

subplot(2,2,3)
polarhistogram(allStimDists(:),nBins)
title('all angles dist to stim loci')
end

end