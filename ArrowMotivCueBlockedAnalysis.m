% ArrowMotivCueBlockedAnalysis

% check dimensions so it will run the same on just 1 person
clear
close all

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files = what(dataFolder);%find files there

isOutput = regexp(files.mat,'ArrowMotivCueBlocked_');%keep only files with this in
isOutput = ~cellfun(@isempty,isOutput);
files.mat = files.mat(isOutput);

dataInds = 1:length(files.mat);%change this to only analyse certain files
dataFiles = fullfile(dataFolder, files.mat(dataInds));%make path to files

dataStruct = ArrowDataLoad(dataFiles);%analyse each file
n = length(dataStruct);%num pps
nTrials = dataStruct(1).nTrials;
nTT = dataStruct(1).nTrialTypes;
%% which measures to plot

inds = [1 4];%which measures to plot 1:4 or [1 4]
nInds = length(inds);

%% format

allPrecTrials = abs(rad2deg(nancat(1,dataStruct.precision)));
allStartRTTrials = nancat(1, dataStruct.startRTs);
allEndRTTrials = nancat(1, dataStruct.endRTs);
allTotalRTTrials = allStartRTTrials + allEndRTTrials;

cueTime = nancat(1, dataStruct.cueTime);
rewLevel = nancat(1, dataStruct.rewardLevel);
isCatchTrial = nancat(1, dataStruct.isCatchTrial); % only shown after response, should not affect

trialType = rewLevel + cueTime; % 2=preLow, 3=postLow, 51=preHi, 52=postHi
v = [2 3 51 52];
for i = 1:4 % replace with 1:4
    trialType(trialType==v(i)) = i;
end

% split each variable by trial type
allPrec = groupMeans(allPrecTrials, 2, trialType, 'dim');
allStartRT = groupMeans(allStartRTTrials, 2, trialType, 'dim');
allEndRT = groupMeans(allEndRTTrials, 2, trialType, 'dim');
allTotalRT = groupMeans(allTotalRTTrials, 2, trialType, 'dim');

allMeasuresTrialSep = nancat(4, allPrec, allStartRT, allEndRT, allTotalRT);
allMeasures = sq(nanmedian(allMeasuresTrialSep,3));
allMeasures(:,:,1) = sq(nanmean(allMeasuresTrialSep(:,:,:,1),3));

allMeasuresLog = sq(nanmean(log(allMeasuresTrialSep),3));
allMeasuresLog(:,:,1) = sq(nanmean(allMeasuresTrialSep(:,:,:,1),3));

allMeasuresTrial = nancat(3, abs(rad2deg([dataStruct.precision]))', [dataStruct.startRTs]',...
    [dataStruct.endRTs]', [dataStruct.startRTs]' + [dataStruct.endRTs]');


allCatchCorrTrials = nancat(1, dataStruct.catchCorrect);
allCatchCorr = groupMeans(allCatchCorrTrials, 2, trialType,'dim');
allCatchRTs = nancat(1, dataStruct.catchRTs);
allCatchRT = groupMeans(allCatchRTs, 2, trialType, 'dim');

nTrialsPerPP = sum(sum(~isnan(allPrec),3),2);
ppNums = [];
for i = 1:n
    ppNums = [ppNums; repmat(i,nTrialsPerPP(i),1)];
end

dvLabels = {'error','startRT','endRT','totalRT'};

%% anova

% allMeasures = nancat(3, nanmean(allPrec,3), nanmean(allStartRT,3), nanmean(allEndRT,3), nanmean(allRT,3));%conc
% allMeasuresTrial = nancat(3, rad2deg(abs([dataStruct.precision]')), ...
%     [dataStruct.startRTs]', [dataStruct.endRTs]', [dataStruct.startRTs]' + [dataStruct.endRTs]');

x = zeros(size(allMeasures(:,:,1)));%zeros
rew = x;
rew(:,3:4) = 1;
postcue = x;
postcue(:,[2 4]) = 1;

ppInd = repmat([1:n]',1,4);%participant indices 1:n repeated for each condition
modelMat = [0 0 1; 1 0 0;0 1 0;...
    1 1 0;];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors = {col(rew), col(postcue), col(ppInd)};%combine factors into cells - IVs and participant
factorLabels = {'rewardLevel','cueTime','pp'};%labels for factors


for i = 1:nInds
%     DVVec =  col(allMeasuresLog(:,:,i));%make a matrix into a vector
%     [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',3,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms

    rmStats{i} = rmanova(reshape(allMeasuresLog(:,:,inds(i)), [],2,2), {'pp','cue','reward'});
    pVals(:,i) = rmStats{i}.pValue;
end

%% trials

trialRew = [dataStruct.rewardLevel]';
trialCue = [dataStruct.cueTime]';
trialPP = ppNums;%reshape(repmat([1:n],[nTrials,1]),[],1);
modelMat = [0 0 1; 1 0 0; 0 1 0;...
    1 1 0; ];
trialFactors = {trialRew, trialCue, trialPP};
factorLabels = {'rewardLevel','cueTime','pp'};


for i = 1:nInds
    DVVec = [allMeasuresTrial(:,:,inds(i))];
    if regexp(dvLabels{inds(i)},'RT')
        DVVec = log(DVVec);
    end
    [statsTrial(i).p,statsTrial(i).tab,statsTrial(i).s] = anovan(DVVec,trialFactors,'varnames',factorLabels,'display','off','model','full','random',3,'model',modelMat);
end



%% combine stats into one table and display

% pVals = [stats.p];
pValsTrial = [statsTrial.p];
alphas = [.05, .01, .001, .0001];

nStars = zeros(3,nInds);
nStarsTrial = nStars;
for i = 1:nInds
    nStars(:,i) = sum(pVals(2:end,i) < alphas,2);
    nStarsTrial(:,i) = sum(pValsTrial(2:end,i) < alphas,2);
    for j = 1:3
        stars{j,i} = repmat('*',1,nStars(j,i));
        starsTrial{j,i} = repmat('*',1,nStarsTrial(j,i));
    end
end

allStats = table();
allStats1 = table();


for i = 1:nInds
    allStats1.measure = repmat(dvLabels(inds(i)),3,1);
    allStats1.effect = rmStats{i}.Term(2:end);
    allStats1.df = rmStats{i}.DF1(2:end);
    allStats1.dfErr = rmStats{i}.DF2(2:end);
    allStats1.F = rmStats{i}.FStat(2:end);
    allStats1.p = rmStats{i}.pValue(2:end);
    allStats1.stars = stars(:,i);
      
    allStats1.dfTrial = [statsTrial(i).tab{3:5,3}]';
    allStats1.dfErrTrial = [statsTrial(i).tab{3:5,11}]';
    allStats1.FTrial = [statsTrial(i).tab{3:5,6}]';
    allStats1.pTrial = [statsTrial(i).tab{3:5,7}]';
    allStats1.starsTrial = starsTrial(:,i);
    
    allStats = vertcat(allStats, allStats1);
end  
    

disp(allStats(:,[1,2,6,7]))%,11,12]))



%% main figs
yLabs = {'error (deg)', 'start RT (ms)', 'end RT (ms)', 'total RT (ms)'};

figure()
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 3; 2 4];
for i = 1:nInds
    if nInds > 1
        subplot(nInds/2,2,i)
    end
    h = errorBarPlot(nancat(3, allMeasures(:,condInds(1,:),inds(i)), ...
        allMeasures(:,condInds(2,:),inds(i))),'type','line');
    hold on
    for j = 1:2
        h(j).LineWidth = 2;
        h(j).Marker = markers{j};
        h(j).Color = lineColours(j,:);
    end
    set(gca,'XTick',1:4,'XTickLabel',{'Low','High'})
    xlabel('Reward')
    xlim([0.5 2.5])
    ylabel(yLabs{inds(i)})
    if i==1
        legend(fliplr(h), fliplr({'Pre','Post'}),'Location','Best')
    end
    box off

end

%%
allMeanConds = NaN(n,2,2,4);
cols = {[1 3]; [2 4];...
        [1 2]; [3 4];...
        };
for i = 1:4%for each DV
    for j = 1:2%for each factor
        %get means
        allMeanConds(:,:,j,i) = [nanmean(allMeasures(:,cols{j*2-1},i),2), nanmean(allMeasures(:,cols{j*2},i),2)];
    end
    
    %get ylims for each measure
    yLims (i,:) = [min(nanmean(allMeanConds(:,:,:,i),1),[],'all'),max(nanmean(allMeanConds(:,:,:,i),1),[],'all')] .* [.9 1.1];
end


%% plot main effects
xFacLabels = {'pre','post';...
            'low','high';...
            };
figure()

for i = 1:nInds
    for j = 1:2
        subplot(nInds,2,(i-1)*2+j)
        errorBarPlot(allMeanConds(:,:,j,inds(i)),'type','bar');
        set(gca,'XTick',1:2,'XTickLabel',xFacLabels(j,:))
        if i==1
            title(factorLabels{3-j})
        end
        if j==1 
            ylabel(dvLabels{inds(i)})
        end
        ylim(yLims(inds(i),:))
        text(0.5, 0.8, stars(j,i),'Units','normalized','FontSize',14)

    end
end



%% effect size for pps

figure()

effects = allMeanConds(:,2,:,:) - allMeanConds(:,1,:,:);%high error - low error
%

for i = 1:nInds
    
   yLims(inds(i),:) = repmat( max(abs([min(effects(:,:,:,inds(i)),[],'all'), max(effects(:,:,:,inds(i)),[],'all')])),1,2) .* [-1.1 1.1];%get ylimits 

    for j = 1:2
        
%         yLims = repmat( max(abs([min(effects,[],'all'), max(effects,[],'all')])),1,2) .* [-1.1 1.1];%get ylimits
        subplot(nInds,2,(i-1)*2+j)
        bar(sort(effects(:,:,j,inds(i)),'ascend'));
        if j==1
            ylabel(['\Delta ' dvLabels{inds(i)}])
        end
        if i==1
            title(factorLabels{3-j})
        end
        ylim(yLims(inds(i),:))
    end
end




%% catch trial accuracy

% precCatchCorr = groupMeans(allPrec, 3, allCatchCorr);
% meanCatchCorr = nanmean(allCatchCorr, 3) .* 100;
% [preLo, postLo, preHi, postHi]
catchMeasures = nancat(3, nanmean(allCatchCorr,3), nanmedian(allCatchRT,3));
for i = 1:2
    catchStats{i} = rmanova(reshape(catchMeasures(:,:,i),[],2,2),{'pp','cue','reward'});
end

[~,p] = ttest(catchMeasures(:,:,1)-.5); %is catch acc sig diff from chance?

% draw fig

yLabs = {'% coin catch accuracy', 'median catch RT (ms)'};
figure();
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 2; 3 4];
for i = 1:2
    subplot(1,2,i)
    h = errorBarPlot(nancat(3, catchMeasures(:,condInds(1,:),i), ...
        catchMeasures(:,condInds(2,:),i)),'type','line');
    hold on
    for j = 1:2
        h(j).LineWidth = 2;
        h(j).Marker = markers{j};
        h(j).Color = lineColours(j,:);
    end
    set(gca,'XTick',1:2,'XTickLabel',{'pre','post'})
    xlabel('cue time')
    xlim([0.5 2.5])
    ylabel(yLabs{i})
    if i==1
        legend({'low','high'},'Location','Best')
    end
    box off
end

%% does acc differ by catch accuracy?
tt = nancat(1, dataStruct.trialTypes);
precCond = groupMeans(allPrecTrials, 2, tt, 'dim');
catchCond = groupMeans(allCatchCorrTrials, 2, tt, 'dim');

precByCatch = groupMeans(precCond, 3, catchCond, 'dim');
precByCatch(:,1:2:end,:,:) = []; % remove noCatch conds
% conds = [preLo, preHi, postLo, postHi]

meanPrecByCatch = nanmean(precByCatch, 4);

precByCatchStats = rmanova(reshape(meanPrecByCatch,[],2,2,2), {'pp','cue','rew','catchAcc'});

figure();
titles = {'catch incorrect', 'catch correct'};
for i = 1:2
    subplot(1,2,i)
    h = errorBarPlot(nancat(3, meanPrecByCatch(:,condInds(1,:),i), ...
        meanPrecByCatch(:,condInds(2,:),i)),'type','line');
    hold on
    for j = 1:2
        h(j).LineWidth = 2;
        h(j).Marker = markers{j};
        h(j).Color = lineColours(j,:);
    end
    set(gca,'XTick',1:2,'XTickLabel',{'low','high'})
    xlabel('reward')
    xlim([0.5 2.5])
    ylabel('error (deg)')
    if i==1
        legend({'pre','post'},'Location','North')
    end
    box off
    title(titles{i});
end

makeSubplotScalesEqual(1,2);
%%

save('./Data/ArrowMotivCueBlockedAnalysis.mat')

% end