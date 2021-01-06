% ArrowMotivReallocateAnalysis

% check dimensions so it will run the same on just 1 person
clear
close all

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files = what(dataFolder);%find files there

isOutput = regexp(files.mat,'ArrowMotivReallocate');%keep only files with this in
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
allMeasuresTrial = nancat(3, [dataStruct.precision]', [dataStruct.startRTs]',...
    [dataStruct.endRTs]', [dataStruct.startRTs]' + [dataStruct.endRTs]');

% split each variable by trial type
allPrec = permute( abs(rad2deg(nancat(3,dataStruct.prec))), [3,2,1]) ;
allStartRT = permute( nancat(3,dataStruct.startRT), [3,2,1]);
allEndRT = permute( nancat(3,dataStruct.endRT), [3,2,1]);
allRT = allStartRT + allEndRT;

allMeasures = nancat(3, nanmean(allPrec,3), nanmedian(allStartRT,3), nanmedian(allEndRT,3), nanmedian(allRT,3));%conc

allMeasuresLog = nancat(3, nanmean(allPrec,3), nanmean(log(allStartRT),3), nanmean(log(allEndRT),3), nanmean(log(allRT),3));%conc

cueTime = nancat(2, dataStruct.cueTime)';
targRew = nancat(2, dataStruct.targRewLevel)'; % rew level of target
distRew = nancat(2, dataStruct.distRewLevel)'; 
isMotiv = targRew == distRew;

trialTypes = nancat(2, [dataStruct.trialTypes])'; % [preLL, preLH, preHL, preHH, post...] (TD)
trialFactors = [cueTime, targRew, distRew, trialTypes];

nTrialsPerPP = sum(sum(~isnan(allPrec),3),2);
ppNums = [];
for i = 1:n
    ppNums = [ppNums; repmat(i,nTrialsPerPP(i),1)];
end

dvLabels = {'error','startRT','endRT','totalRT'};
%% anova

x = zeros(n,1);%zeros
motiv = [x+1;x;x;x+1;x+1;x;x;x+1];%1 = motiv, 0=reallo
rew = [x;x;x+1;x+1;x;x;x+1;x+1];%0=low, 1=hi
cue = [x;x;x;x;x+1;x+1;x+1;x+1];%0=pre, 1=post
ppInd = repmat([1:n]',nTT,1);%participant indices 1:n repeated for each condition
modelMat = [0 0 0 1; 1 0 0 0;0 1 0 0; 0 0 1 0;...
    1 1 0 0; 1 0 1 0; 0 1 1 0;...
    1 1 1 0];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors = {rew, cue, motiv, ppInd};%combine factors into cells - IVs and participant
factorLabels = {'rewardLevel','cueTime','isEqual','pp'};%labels for factors


allMeasuresRM = allMeasuresLog(:,[2 3 1 4 6 7 5 8],:);
for i = 1:nInds
%     DVVec =  reshape(allMeasures(:,:,i),1,nTT*n)';%make a matrix into a vector
%     [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',4,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms

%     [motivPreLo, realloPreLo, realloPreHi, motivPreHi, motivPostLo,
%     reallocPostLo, realloPostHi, motivPostHi]
%  need to get it [rew, cue, motiv]
    % reorder columns
    rmStats{i} = rmanova(reshape(allMeasuresRM(:,:,i), [],2,2,2), {'pp','reward','equal','cue'});
    pVals(:,i) = rmStats{i}.pValue;

end

%% trials


trialPP = ppNums;
modelMat = [0 0 0 1; 1 0 0 0;0 1 0 0; 0 0 1 0;...
    1 1 0 0; 1 0 1 0; 0 1 1 0;...
    1 1 1 0];
trialFactors = {targRew, cueTime, isMotiv, trialPP};


for i = 1:nInds
    DVVec = [allMeasuresTrial(:,:,i)];
%     if regexp(dvLabels{inds(i)},'RT')
%         DVVec = log(DVVec);
%     end
    [statsTrial(i).p,statsTrial(i).tab,statsTrial(i).s] = anovan(DVVec,trialFactors,'varnames',factorLabels,'display','off','model','full','random',4,'model',modelMat);
end



%% combine stats into one table and display

% pVals = [stats.p];
pValsTrial = [statsTrial.p];
alphas = [.05, .01, .001, .0001];

nStars = zeros(7,nInds);
nStarsTrial = nStars;
for i = 1:nInds
    nStars(:,i) = sum(pVals(2:end,i) < alphas,2);
    nStarsTrial(:,i) = sum(pValsTrial(2:end,i) < alphas,2);
    for j = 1:7
        stars{j,i} = repmat('*',1,nStars(j,i));
        starsTrial{j,i} = repmat('*',1,nStarsTrial(j,i));
    end
end

allStats = table();
allStats1 = table();


for i = 1:nInds
    allStats1.measure = repmat(dvLabels(inds(i)),7,1);
    allStats1.effect = rmStats{i}.Term(2:end);
    allStats1.df = rmStats{i}.DF1(2:end);
    allStats1.dfErr = rmStats{i}.DF2(2:end);
    allStats1.F = rmStats{i}.FStat(2:end);
    allStats1.p = rmStats{i}.pValue(2:end);
    allStats1.stars = stars(:,i);
    
    allStats1.dfTrial = [statsTrial(i).tab{3:9,3}]';
    allStats1.dfErrTrial = [statsTrial(i).tab{3:9,11}]';
    allStats1.FTrial = [statsTrial(i).tab{3:9,6}]';
    allStats1.pTrial = [statsTrial(i).tab{3:9,7}]';
    allStats1.starsTrial = starsTrial(:,i);
    
    allStats = vertcat(allStats, allStats1);
end  
    

disp(allStats(:,[1,2,6,7]))%,11,12]))


%% also do stats on equal/unequal separately
% anova
allMeasuresConds = nancat(4, allMeasures(:,[1 4 5 8],:), allMeasures(:,[2 3 6 7],:));

x = zeros(n,1);%zeros
rew = [x;x+1;x;x+1];%1=low, 2=hi
cue = [x;x;x+1;x+1];%1=pre, 2=post
ppInd = repmat([1:n]',nTT./2,1);%participant indices 1:n repeated for each condition
modelMat2 = [0 0 1; 1 0 0;0 1 0; ...
    1 1 0 ; ];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors2 = {rew, cue, ppInd};%combine factors into cells - IVs and participant
factorLabels2 = {'rewardLevel','cueTime','pp'};%labels for factors

for j = 1:2
    for i = 1:nInds
%         DVVec =  reshape(allMeasuresConds(:,:,i,j),1,nTT*n/2)';%make a matrix into a vector
%         [stats2(i,j).p,stats2(i,j).tab,stats2(i,j).s] = anovan(DVVec,factors2,'varnames',factorLabels2,'display','off','model','full','random',3,'model',modelMat2);

        rmStats2{i,j} = rmanova(reshape(sq(allMeasuresConds(:,:,inds(i),j)),[],2,2), {'pp','reward','cue'});
        pVals2(:,i,j) = rmStats2{i,j}.pValue;
    end
end

%% combine stats into one table and display

% pVals = nancat(3, [stats2(:,1).p],[stats2(:,2).p]);
alphas = [.05, .01, .001, .0001];

nStars2 = zeros(3,nInds);
stars2 = cell(3,nInds);
for k = 1:2
    for i = 1:nInds
        nStars2(:,i,k) = sum(pVals2(2:end,i,k) < alphas,2);
        for j = 1:3
            stars2{j,i,k} = repmat('*',1,nStars2(j,i,k));
        end
    end
end

%%
allStats2 = table();
allStats1 = table();

condLabels = {'equal','unequal'; 'low','high';'pre','post'};
for j = 1:2
    for i = 1:nInds
        allStats1.cond = repmat(condLabels(1,j),3,1);
        allStats1.measure = repmat(dvLabels(inds(i)),3,1);
        allStats1.effect = rmStats2{i,j}.Term(2:end);
        allStats1.df = rmStats2{i,j}.DF1(2:end);
        allStats1.dfErr = rmStats2{i,j}.DF2(2:end);
        allStats1.F = rmStats2{i,j}.FStat(2:end);
        allStats1.p = rmStats2{i,j}.pValue(2:end);
        allStats1.stars = stars2(:,i,j);
        
        allStats2 = vertcat(allStats2, allStats1);
    end  
end

disp(allStats2(:,[1,2,3,7,8]))%,11,12]))
%% main figs
yLabs = {'error (deg)', 'start RT (ms)', 'end RT (ms)', 'total RT (ms)'};

figure()
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 4 2 3; 5 8 6 7];
for i = 1:nInds
    if nInds>1
        subplot(ceil(nInds/2),2,i)
    end

    h = errorBarPlot(nancat(3, allMeasures(:,condInds(1,:),inds(i)), ...
                    allMeasures(:,condInds(2,:),inds(i))),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    hold on
    JoinUpErrorPoints(h, [1 2; 3 4]);
    
    set(gca,'XTick',1:4,'XTickLabel',{'Low','High','Low','High'})
    xlabel('Equal                Unequal')
    xlim([0.5 4.5])
    ylabel(yLabs{inds(i)})
    if i==1
        legend(fliplr(h), fliplr({'Pre','Post'}),'Location','Best')
    end
    box off

end


%%
allMeanConds = NaN(n,2,3,nInds);
cols = {[1 2 5 6]; [3 4 7 8];... % low vs high
    [1 2 3 4]; [5 6 7 8];... % pre vs post
    [1 4 5 8];[2 3 6 7];}; % motiv vs reallocate
for i = 1:nInds%for each DV
    for j = 1:3%for each factor
        %get means
        allMeanConds(:,:,j,i) = [nanmean(allMeasures(:,cols{j*2-1},inds(i)),2), nanmean(allMeasures(:,cols{j*2},inds(i)),2)];
    end
    
    %get ylims for each measure
    yLims (i,:) = [min(nanmean(allMeanConds(:,:,:,i),1),[],'all'),max(nanmean(allMeanConds(:,:,:,i),1),[],'all')] .* [.9 1.1];
end


%% plot main effects
xFacLabels = {'low','high';...
    'pre','post';...
    'equal','unequal'};
figure()

for i = 1:nInds
    for j = 1:3
        subplot(nInds,3,(i-1)*3+j)
        errorBarPlot(allMeanConds(:,:,j,i),'type','bar');
        set(gca,'XTick',1:2,'XTickLabel',xFacLabels(j,:))
        if i==1
            title(factorLabels{j})
        end
        if j==1
            ylabel(dvLabels{inds(i)})
        end
        ylim(yLims(i,:))
        text(0.5, 0.8, stars(j,i),'Units','normalized','FontSize',14)

    end
end



%% plot interaction of cue*delay
intCols =  {   [1 2],[5 6],[3 4],[7 8];...%rew*cue
        [1 5], [2 6], [4 8], [3 7];...%rew*motiv
        [1 4], [2 3], [5 8], [6 7];};%cue*motiv

allInts = [];
for i = 1:nInds
    for j = 1:3
        allInts(:,:,:,j,i) = nancat(3, [nanmean(allMeasures(:,intCols{j,1},inds(i)),2),nanmean(allMeasures(:,intCols{j,2},inds(i)),2)],[nanmean(allMeasures(:,intCols{j,3},inds(i)),2),nanmean(allMeasures(:,intCols{j,4},inds(i)),2)]);
    end
%     allInts = nancat(4, allInts, precInt); 
    yLims (i,:) = [min(nanmean(allInts(:,:,:,:,i),1),[],'all'),max(nanmean(allInts(:,:,:,:,i),1),[],'all')] .* [.9 1.1];
end
%%
intLabels = {'low','high';...
    'low','high';...
    'pre','post'};
intLegends = {'pre','post';...
    'equal','unequal';...
    'equal','unequal'};
intTitles = {'rew*cue','rew*equal','cue*equal'};

intCols = {[.1 .8 .2], [.85, .32, .98];...
    [1 .8 .2], [0 .5 .8];...
    [1 .3 0], [.6, .2, .6];};
figure()


for i = 1:nInds
    
    for j = 1:3
        subplot(nInds,3,(i-1)*3+j)
        h = errorBarPlot(permute(allInts(:,:,:,j,i),[1,3,2]), 'type', 'line','plotargs',{'LineWidth',2});
%         h(1).Color = intCols{j,1};
%         h(2).Color = intCols{j,2};
        if i==nInds
            set(gca,'XTick',1:2,'XTickLabel',intLegends(j,:))
        else
            set(gca,'XTickLabels','');
        end
        xlim([.5 2.5])
        ylim([yLims(i,:)])
        text(0.5, 0.8, stars(j+3,i),'Units','normalized','FontSize',14)
        if i==1
            title(intTitles{j})
        end
        if i==nInds && j==3
            legend(intLabels(1,:),'Location','Best')
        end
        if j==1
            ylabel(dvLabels{inds(i)})
        end
        
    end
end



%% effect size for pps

figure()

effects = allMeanConds(:,2,:,:) - allMeanConds(:,1,:,:);%high error - low error
%

for i = 1:nInds
    
   yLims(i,:) = repmat( max(abs([min(effects(:,:,:,i),[],'all'), max(effects(:,:,:,i),[],'all')])),1,2) .* [-1.1 1.1];%get ylimits 

    for j = 1:3
        
%         yLims = repmat( max(abs([min(effects,[],'all'), max(effects,[],'all')])),1,2) .* [-1.1 1.1];%get ylimits
        subplot(nInds,3,(i-1)*3+j)
        bar(sort(effects(:,:,j,i),'ascend'));
        if j==1
            ylabel(['\Delta ' dvLabels{inds(i)}])
        end
        if i==1
            title(factorLabels{j})
        end
        ylim(yLims(i,:))
    end
end



%%

save('./Data/ArrowMotivReallocateAnalysis.mat')

% end