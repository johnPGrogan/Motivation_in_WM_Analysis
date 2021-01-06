% ArrowMotivRetrocueAnalysis

clear
close all

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files = what(dataFolder);%find files there

isOutput = regexp(files.mat,'ArrowMotivRetrocue_');%keep only files with this in
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
% conds = inconLow, conLow, inconHigh, conHigh
condLabels = {'incongruent','congruent'; 'low','high'};
dvLabels = {'error','startRT','endRT','totalRT','cueCorr','cueRT'};

allPrec = permute( abs(rad2deg(nancat(3,dataStruct.prec))), [3,2,1]) ;
allStartRT = permute( nancat(3,dataStruct.startRT), [3,2,1]);
allEndRT = permute( nancat(3,dataStruct.endRT), [3,2,1]);
allRT = allStartRT + allEndRT;

allRetrocueCorr = permute(nancat(3, dataStruct.retrocueCorr),[3,2,1]);
allCueRT = permute( nancat(3,dataStruct.cueRT), [3,2,1]);


%% decide on whether to remove pps or trials by retrocue correct

allRetrocueCorr = permute(nancat(3, dataStruct.retrocueCorr),[3,2,1]);

figure();
bar(sort(nanmean(allRetrocueCorr,[3,2])));
xlabel('person')
ylabel('retrocue accuracy')
% ylim([.9 1]);
box off

disp(nansum(~allRetrocueCorr,[3,2]));

%% does motiv affect acc

retrocueByRew = nancat(3, sq(nanmean( nancat(4, allRetrocueCorr(:,1:2,:),...
    allRetrocueCorr(:,3:4,:)),[3,2])), sq(nanmean( nancat(4, allCueRT(:,1:2,:),...
    allCueRT(:,3:4,:)),[3,2])));



figure();
for i = 1:2
    
    retStats{i} = rmanova(retrocueByRew(:,:,i), {'pp','rew'});
    
    subplot(1,2,i)
    h = errorBarPlot(retrocueByRew(:,:,i),'type','line','plotargs',{'LineWidth',2});
    set(gca,'XTick',1:2,'XTickLabel',condLabels(2,:))
    xlabel('reward')
    xlim([0.5 2.5])
    ylabel(dvLabels{i+4})
    box off
    
    if retStats{i}.pValue(2) < .05
        text(0.5, 0.8, '*','Units','normalized','FontSize',14)
    end

end


%% remove incorrect trials

toRemove = ~allRetrocueCorr;
allPrec(toRemove) = NaN;
allStartRT(toRemove) = NaN;
allEndRT(toRemove) = NaN;
allTotalRT(toRemove) = NaN;
allCueRT(toRemove) = NaN;
% allRetrocueCorr(toRemove) = NaN;

% %% remove people with low acc 
% 
% lowAcc = nanmean(allRetrocueCorr,[3,2]) <.95;
% 
% allPrec(lowAcc,:,:) = NaN;
% allStartRT(lowAcc,:,:) = NaN;
% allEndRT(lowAcc,:,:) = NaN;
% allTotalRT(lowAcc,:,:) = NaN;
% allCueRT(lowAcc,:,:) = NaN;

%%

allMeasuresSep = nancat(4, allPrec, allStartRT, allEndRT, allRT, allRetrocueCorr, allCueRT);
allMeasures(:,:,[1 5]) = sq(nanmean(allMeasuresSep(:,:,:,[1 5]),3));
allMeasures(:,:,[2 3 4 6]) = sq(nanmedian(allMeasuresSep(:,:,:,[2 3 4 6]),3));

allMeasuresLog = allMeasures;
allMeasuresLog(:,:,[2 3 4 6]) = sq(nanmean(log(allMeasuresSep(:,:,:,[2 3 4 6])),3));

allMeasuresTrial = nancat(3, rad2deg(abs([dataStruct.precision]')), ...
    [dataStruct.startRTs]', [dataStruct.endRTs]',...
    [dataStruct.startRTs]' + [dataStruct.endRTs]',...
    [dataStruct.RetrocueResponseCorrect]',[dataStruct.cueRTs]');

trialTypes = nancat(2, [dataStruct.trialTypes])';
retrocueCongruent = nancat(2, [dataStruct.RetrocueValid])';
reward = nancat(2, [dataStruct.rewardLevel])';
trialFactors = [retrocueCongruent, reward, trialTypes];


nTrialsPerPP = sum(sum(~isnan(allPrec),3),2);
ppNums = kron([1:n],ones(1,nTrials))'; % each ppnum repeated


allMeasuresTrial(repmat(~allMeasuresTrial(:,5),[1,1,nInds])) = NaN; % remove retrocue incorrect
%% anova

% allMeasures = nancat(3, nanmean(allPrec,3), nanmean(allStartRT,3), nanmean(allEndRT,3), nanmean(allRT,3));%conc
% allMeasuresTrial = nancat(3, rad2deg(abs([dataStruct.precision]')), ...
%     [dataStruct.startRTs]', [dataStruct.endRTs]', [dataStruct.startRTs]' + [dataStruct.endRTs]');

x = zeros(size(allMeasures(:,:,1)));%zeros
rew = x;
rew(:,3:4) = 1;
cueCongr = x;
cueCongr(:,[2 4]) = 1;

ppInd = repmat([1:n]',1,4);%participant indices 1:n repeated for each condition
modelMat = [0 0 1; 1 0 0;0 1 0;...
    1 1 0;];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors = {col(rew), col(cueCongr), col(ppInd)};%combine factors into cells - IVs and participant
factorLabels = {'rewardLevel','cueCongr','pp'};%labels for factors


for i = 1:nInds
%     DVVec =  col(allMeasuresLog(:,:,inds(i)));%make a matrix into a vector
%     [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',3,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms
    
    rmStats{i} = rmanova(reshape(allMeasures(:,:,inds(i)), n,2,2), {'pp','cong','rew'});
    pVals(:,i) = rmStats{i}.pValue;
end

%% trials

trialRew = [dataStruct.rewardLevel]';
trialCue = [dataStruct.RetrocueValid]';
trialPP = ppNums;%reshape(repmat([1:n],[nTrials,1]),[],1);
modelMat = [0 0 1; 1 0 0; 0 1 0;...
    1 1 0; ];
trialFactors = {trialRew, trialCue, trialPP};
factorLabels = {'rewardLevel','cueCongr','pp'};


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

nStars = zeros(3,4);
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

figure()
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 3; 2 4];
for i = 1:nInds
    if nInds>1
        subplot(ceil(nInds/2),2,i)
    end

    h = errorBarPlot(nancat(3, allMeasures(:,condInds(1,:),inds(i)), ...
        allMeasures(:,condInds(2,:),inds(i))),'type','line');
    hold on
    for j = 1:2
        h(j).LineWidth = 2;
        h(j).Marker = markers{j};
        h(j).Color = lineColours(j,:);
    end
    set(gca,'XTick',1:4,'XTickLabel',condLabels(2,:))
    xlabel('reward')
    xlim([0.5 2.5])
    ylabel(dvLabels{inds(i)})
    if i==nInds
        legend(condLabels(1,:),'Location','Best')
    end
    box off

end

%%
allMeanConds = NaN(n,2,2,nInds);
cols = {[1 2]; [3 4];...
    [1 3]; [2 4];...
    };
for i = 1:nInds%for each DV
    for j = 1:2%for each factor
        %get means
        allMeanConds(:,:,j,i) = [nanmean(allMeasures(:,cols{j*2-1},inds(i)),2), nanmean(allMeasures(:,cols{j*2},inds(i)),2)];
    end
    
    %get ylims for each measure
    yLims (i,:) = [min(nanmean(allMeanConds(:,:,:,i),1),[],'all'),max(nanmean(allMeanConds(:,:,:,i),1),[],'all')] .* [.9 1.1];
end


%% plot main effects
figure()

for i = 1:nInds
    for j = 1:2
        subplot(nInds,2,(i-1)*2+j)
        errorBarPlot(allMeanConds(:,:,j,i),'type','bar');
        set(gca,'XTick',1:2,'XTickLabel',condLabels(3-j,:))
        if i==1
            title(factorLabels{j})
        end
        if j==1
            ylabel(dvLabels{inds(i)})
        end
        ylim(yLims(i,:))
        text(0.5, 0.8, stars(3-j,i),'Units','normalized','FontSize',14)

    end
end



%% effect size for pps

figure()

effects = allMeanConds(:,2,:,:) - allMeanConds(:,1,:,:);%high error - low error
%

for i = 1:nInds
    
   yLims(inds(i),:) = repmat( max(abs([min(effects(:,:,:,i),[],'all'), max(effects(:,:,:,i),[],'all')])),1,2) .* [-1.1 1.1];%get ylimits 

    for j = 1:2
        
%         yLims = repmat( max(abs([min(effects,[],'all'), max(effects,[],'all')])),1,2) .* [-1.1 1.1];%get ylimits
        subplot(nInds,2,(i-1)*2+j)
        bar(sort(effects(:,:,j,i),'ascend'));
        if j==1
            ylabel(['\Delta ' dvLabels{inds(i)}])
        end
        if i==1
            title(factorLabels{j})
        end
        ylim(yLims(inds(i),:))
    end
end

%%

save('./Data/ArrowMotivRetrocueAnalysis.mat')

% end