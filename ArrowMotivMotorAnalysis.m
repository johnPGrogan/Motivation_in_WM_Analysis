function ArrowMotivMotorAnalysis(exptVersion)
% load ArrowMotivMotorControl data 
% plot average precision and RT over trials, and between pps
% look at speed acc tradeoff
% look at difference between reward trials

if ~(strcmp(exptVersion,'A') || strcmp(exptVersion,'B'))
    error('exptVersion must be A or B')
end

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files = what(dataFolder);%find files there

isOutput = regexp(files.mat,'ArrowMotivMotorControl_');%keep only files with this in
isOutput = ~cellfun(@isempty,isOutput);
files.mat = files.mat(isOutput);


dataInds = 1:length(files.mat);%change this to only analyse certain files
dataFiles = fullfile(dataFolder, files.mat(dataInds));%make path to files

dataStruct = ArrowDataLoad(dataFiles);%analyse each file

%% analyse 1st or 2nd version


nTr = [dataStruct.nTrials]';

if exptVersion=='A'
    dataStruct(nTr~=24) = [];
    dataFiles (nTr~=24) = [];
    
elseif exptVersion=='B'
    dataStruct(nTr~=60) = [];
    dataFiles (nTr~=60) = [];
end


n = length(dataStruct);%num pps
nTrials = max([dataStruct.nTrials]);
nTT = dataStruct(1).nTrialTypes;




%% which vars to look at

inds = [1 4];%which measures to plot 1:4 or [1 4]
nInds = length(inds);

%% format vars across pps

prec = abs(rad2deg(nancat(1, dataStruct.precision)));
startRT = nancat(1, dataStruct.startRTs);
endRT = nancat(1, dataStruct.endRTs);
allRT = startRT + endRT;
allMeasures = nancat(3, prec, startRT, endRT, allRT);

if exptVersion=='a' % remove first 2 trials as they have huge errors and RT - practice effects
    allMeasures(:,1:2, :) = NaN;

else % remove 1 outlier trial with huge error
    outlier = prec==max(prec,[],'all'); % 126 degrees
%     allMeasures(repmat(outlier,1,1,4)) = NaN; % remove
end

dvLabels = {'precision','startRT','endRT','totalRT'};
unitLabels = {'abs err (deg)', 'start RT (ms)', 'end RT (ms)', 'total RT (ms)'};
%% plot acc over trials


figure()
for i = 1:nInds
    if nInds>1
        subplot(ceil(nInds/2), 2, i)
    end
    errorBarPlot(allMeasures(:,:,inds(i)), 'type','line','area',1,...
        'plotargs',{'LineWidth',2},'plotIndividuals',1);
    title(dvLabels{inds(i)})
    xlabel('trial')
    ylabel(unitLabels{inds(i)})
    box off
end

%% look at mean acc

figure()
for i = 1:nInds
    if nInds>1
        subplot(ceil(nInds/2), 2, i)
    end
    bar(nanmean(allMeasures(:,:,inds(i)),2))
    title(dvLabels{inds(i)})
    xlabel('person')
    ylabel(unitLabels{inds(i)})
    box off
end



%% difference in reward levels

rewardLevel = nancat(1, dataStruct.rewardLevel);

allMeasures2 = allMeasures;
% allMeasures2(:,1:2,:) = NaN; % remove first 5 trials
allRewMeasures = permute(groupMeans(allMeasures2, 2, repmat(rewardLevel, [1 1 4]),'dim'), [1 2 4 3]);

allRewMeasuresLog = allRewMeasures;
allRewMeasuresLog(:,:,:,2:4) = log(allRewMeasures(:,:,:,2:4));
%%

x = ones(n,2,nTrials/nTT);%zeros
rew = x;
rew(:,1,:) = 0;
rewTrial = col(rew);

ppInd = repmat([1:n]', 1, 2, nTrials./2);
ppIndTrial = col(ppInd);
% ppInd = repmat([1:n]',nTT,1);%participant indices 1:n repeated for each condition
modelMat = [0 1; 1 0 ;]; % random pp + reward main

factorsTrial = {rewTrial,  ppIndTrial};%combine factors into cells - IVs and participant
factors = {col(nanmean(rew,3)), col(nanmean(ppInd,3))}; % average across trials
factorLabels = {'rewardLevel','pp'};%labels for factors



for i = 1:nInds
    
    DVVec =  col(allRewMeasures(:,:,:,inds(i)));%make a matrix into a vector
    %anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms
    [statsTrial(i).p,statsTrial(i).tab,statsTrial(i).s] = anovan(DVVec,factorsTrial,'varnames',factorLabels,'display','off','model','full','random',2,'model',modelMat);
    
    if inds(i)==1
        allRewSummary(:,:,i) = sq(nanmean(allRewMeasures(:,:,:,inds(i)),3));
    else
        allRewSummary(:,:,i) = sq(nanmedian(allRewMeasures(:,:,:,inds(i)),3));
    end
    rmStats{i} = rmanova(allRewSummary(:,:,i), {'pp','rew'});
    pVals(:,i) = rmStats{i}.pValue;

    % average across trials and repeat anova
    DVVec = col(nanmean(allRewMeasures(:,:,:,inds(i)), 3));
    [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',2,'model',modelMat);
end

%%

% pVals = [stats.p];
pValsTrial = [statsTrial.p];
alphas = [.05, .01, .001, .0001];

[nStars,nStarsTrial] = deal(zeros(1,nInds));
for i = 1:nInds
    nStars(:,i) = sum(pVals(2:end,i) < alphas,2);
    nStarsTrial(:,i) = sum(pValsTrial(2:end,i) < alphas,2);
    for j = 1
        stars{j,i} = repmat('*',1,nStars(j,i));
        starsTrial{j,i} = repmat('*',1,nStarsTrial(j,i));
    end
end

allStats = table();
allStats1 = table();


for i = 1:nInds
    allStats1.measure = repmat(dvLabels(inds(i)),1,1);
    allStats1.effect = rmStats{i}.Term(2:end);
    allStats1.df = rmStats{i}.DF1(2:end);
    allStats1.dfErr = rmStats{i}.DF2(2:end);
    allStats1.F = rmStats{i}.FStat(2:end);
    allStats1.p = rmStats{i}.pValue(2:end);
    allStats1.stars = stars(:,i);

    allStats1.dfTrial = [statsTrial(i).tab{3,3}]';
    allStats1.dfErrTrial = [statsTrial(i).tab{3,11}]';
    allStats1.FTrial = [statsTrial(i).tab{3,6}]';
    allStats1.pTrial = [statsTrial(i).tab{3,7}]';
    allStats1.starsTrial = starsTrial(:,i);

    allStats = vertcat(allStats, allStats1);
end  
    

disp(allStats(:,[1,2,6,7]))%,11,12]))

%%
markers = {'o','^'};
figure()
for i = 1:nInds
    if nInds>1
        subplot(ceil(nInds/2), 2, i)
    end
    for j = 1
        if i==1
            h(j) = errorBarPlot(nanmean(allRewMeasures(:,:,:,inds(i)),3),...
                'type','line','plotargs',{'LineWidth',2.5,'Marker',markers{j}});
            hold on
        else
            h(j) = errorBarPlot(nanmedian(allRewMeasures(:,:,:,inds(i)),3),...
                'type','line','plotargs',{'LineWidth',2.5,'Marker',markers{j}});
            hold on
        end
%         h(j).XData(3-j) = NaN;
    end

    title(dvLabels{inds(i)})
    set(gca,'XTick',1:2,'XTickLabel',{'Low','High'})
    xlabel('Reward')
    ylabel(unitLabels{inds(i)})
    box off
    xlim([0.5 2.5])
end


%% effect size per pp
effect = diff(nanmean(allRewMeasures,3), 1, 2);
effect(:,:,:,2:4) = diff(nanmedian(allRewMeasures(:,:,:,2:4),3), 1, 2);
effect = sq(effect);

figure();
for i = 1:nInds
    subplot(nInds/2,2,i)
    % diff in rew levels per pp, averaged over trials
    bar(effect(:,i))
    xlabel('pp')
    ylabel('hi - low rew')
    box off
    ylim(ylim .* 1.2)
    title(dvLabels(inds(i)))
end

%% sensitivity analysis on removing beginning trials

if exptVersion=='a'
    p = [];
    for i = 1:10
        precByRew = groupMeans(prec(:,i:end), 2, rewardLevel(:,i:end));
        m = rmanova(precByRew,{'pp','rew'});
        p(1,i) = m.pValue(end);
        
        startRTByRew = groupMeans(startRT(:,i:end), 2, rewardLevel(:,i:end));
        m = rmanova(startRTByRew,{'pp','rew'});
        p(2,i) = m.pValue(end);
        
        endRTByRew = groupMeans(endRT(:,i:end), 2, rewardLevel(:,i:end));
        m = rmanova(endRTByRew,{'pp','rew'});
        p(3,i) = m.pValue(end);
    end
end
        
        
%%

save(sprintf('./Data/ArrowMotivMotorAnalysis_%c.mat', exptVersion))