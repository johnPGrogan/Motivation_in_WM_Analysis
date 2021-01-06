% ArrowMotivAnalysis

% check dimensions so it will run the same on just 1 person
clear
close all

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files = what(dataFolder);%find files there

isOutput = regexp(files.mat,'ArrowMotivCue_');%keep only files with this in
isOutput = ~cellfun(@isempty,isOutput);
files.mat = files.mat(isOutput);


dataInds = 1:length(files.mat);%change this to only analyse certain files
dataFiles = fullfile(dataFolder, files.mat(dataInds));%make path to files

dataStruct = ArrowDataLoad(dataFiles);%analyse each file
n = length(dataStruct);%num pps
nTrials = dataStruct(1).nTrials;
nTT = dataStruct(1).nTrialTypes;
%% which measures to plot


inds = [1 4];%which measures to plot 1:4 or [1 4], for just error + totalRT
nInds = length(inds);




%% format

% [shortPrelo, longPreLo, shortPostLo, longPostLo, shortPreHi, longPreHi, shortPostHi, longPostHi]
allPrec = permute( abs(rad2deg(reshape([dataStruct.prec],nTrials/nTT,nTT,n))), [3,2,1]) ;
allStartRT = permute( reshape([dataStruct.startRT],nTrials/nTT,nTT,n), [3,2,1]);
allEndRT = permute( reshape([dataStruct.endRT],nTrials/nTT,nTT,n), [3,2,1]);
allRT = allStartRT + allEndRT;

prec = abs(rad2deg(nancat(1, dataStruct.precision)));
startRT = abs(rad2deg(nancat(1, dataStruct.startRTs)));
endRT = abs(rad2deg(nancat(1, dataStruct.endRTs)));
totRT = startRT + endRT;

dvLabels = {'error','startRT','endRT','totalRT'};

%% anova

allMeasures = nancat(3, nanmean(allPrec,3), nanmedian(allStartRT,3), nanmedian(allEndRT,3), nanmedian(allRT,3));%conc
allMeasuresLog = nancat(3, nanmean(allPrec,3), nanmean(log(allStartRT),3), nanmean(log(allEndRT),3), nanmean(log(allRT),3));%conc
allMeasuresTrial = nancat(3, rad2deg(abs([dataStruct.precision]')), ...
    [dataStruct.startRTs]', [dataStruct.endRTs]', [dataStruct.startRTs]' + [dataStruct.endRTs]');

x = zeros(n,1);%zeros
del = [x;x+1;x;x+1;x;x+1;x;x+1];%set to 1 when reward was high, 0 for low
cue = [x;x;x+1;x+1;x;x;x+1;x+1];%pre/post cue
rew = [x;x;x;x;x+1;x+1;x+1;x+1];%short/long delay
ppInd = repmat([1:n]',nTT,1);%participant indices 1:n repeated for each condition
modelMat = [0 0 0 1; 1 0 0 0;0 1 0 0; 0 0 1 0;...
    1 1 0 0; 1 0 1 0; 0 1 1 0;...
    1 1 1 0];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors = {rew, cue, del, ppInd};%combine factors into cells - IVs and participant
factorLabels = {'rewardLevel','cueTime','delayLength','pp'};%labels for factors


for i = 1:nInds
%     DVVec =  reshape(allMeasuresLog(:,:,i),1,nTT*n)';%make a matrix into a vector
%     [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',4,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms
    
    rmStats{i} = rmanova(reshape(allMeasuresLog(:,:,inds(i)), [],2,2,2), {'pp','delay','cue','reward'});
    pVals(:,i) = rmStats{i}.pValue;
end

%% trials


trialRew = [dataStruct.rewardLevel]';
trialCue = [dataStruct.cueTime]';
trialDel = [dataStruct.delayLength]';
trialPP = reshape(repmat([1:n],[nTrials,1]),[],1);
modelMat = [0 0 0 1; 1 0 0 0;0 1 0 0; 0 0 1 0;...
    1 1 0 0; 1 0 1 0; 0 1 1 0;...
    1 1 1 0];
trialFactors = {trialDel, trialCue, trialRew, trialPP};
factorLabels = {'delayLength','cueTime','rewardLevel','pp'};


for i = 1:nInds
    DVVec = [allMeasuresTrial(:,:,inds(i))];
    if regexp(dvLabels{inds(i)},'RT')
        DVVec = log(DVVec);
    end
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



%% main figs
yLabs = {'error (deg)', 'start RT (ms)', 'end RT (ms)', 'total RT (ms)'};
figure();
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 5 3 7; 2 6 4 8];
for i = 1:nInds
    if nInds > 1
        subplot(ceil(nInds/2),2,i)
    end

    set(gca,'ColorOrder',lineColours);
    h = errorBarPlot(nancat(3, allMeasures(:,condInds(1,:),inds(i)), ...
                    allMeasures(:,condInds(2,:),inds(i))),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    hold on
    JoinUpErrorPoints(h, [1 2; 3 4]);
%     for j = 1:length(h)
%         h(j).LineWidth = 2;
%         h(j).Marker = markers{j};
%         h(j).Color = lineColours(j,:);
%         if size(condInds,2) == 4
%             h(j).LineStyle = 'none';
%             plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:));
%             plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:));
%         end
%     end
% 
%     hold on
%     for j = 1:2
%         h(j).LineWidth = 2;
%         h(j).Marker = markers{j};
%         h(j).Color = lineColours(j,:);
%     end
    
    set(gca,'XTick',1:4,'XTickLabel',{'Low','High','Low','High'})
    xlabel('Pre                Post')
    xlim([0.5 4.5])
    ylabel(yLabs{inds(i)})
    if i==1
        legend(fliplr(h), fliplr({'Short','Long'}),'Location','Best')
    end
    box off

end

%%
allMeanConds = NaN(n,2,3,nInds);
cols = {[1 3 5 7];[2 4 6 8];...
    [1 2 5 6]; [3 4 7 8];...
    [1 2 3 4]; [5 6 7 8];};
for i = 1:nInds%for each DV
    for j = 1:3%for each factor
        %get means
        allMeanConds(:,:,j,i) = [nanmean(allMeasures(:,cols{j*2-1},inds(i)),2), nanmean(allMeasures(:,cols{j*2},inds(i)),2)];
    end
    
    %get ylims for each measure
    yLims (i,:) = [min(nanmean(allMeanConds(:,:,:,i),1),[],'all'),max(nanmean(allMeanConds(:,:,:,i),1),[],'all')] .* [.9 1.1];
end


%% plot main effects
xFacLabels = {'short','long';...
    'pre','post';...
    'low','high';};
figure()

for i = 1:nInds
    for j = 1:3
        subplot(nInds,3,(i-1)*3+j)
        errorBarPlot(allMeanConds(:,:,j,i),'type','bar');
        set(gca,'XTick',1:2,'XTickLabel',xFacLabels(j,:))
        if i==1
            title(factorLabels{4-j})
        end
        if j==1
            ylabel(dvLabels{inds(i)})
        end
        ylim(yLims(i,:))
        text(0.5, 0.8, stars(j,i),'Units','normalized','FontSize',14)

    end
end



%% plot interaction of cue*delay
intCols =  { [1 5], [2 6], [3 7], [4 8];  ...%rew*cue
        [1 3], [2 4], [5 7], [6 8];...%rew*delay
        [1 2],[3 4],[5 6],[7 8]; };%cue*delay

allInts = [];
for i = 1:nInds
    for j = 1:3
        allInts(:,:,:,j,i) = nancat(3, [nanmean(allMeasures(:,intCols{j,1},inds(i)),2), nanmean(allMeasures(:,intCols{j,2},inds(i)),2)], [nanmean(allMeasures(:,intCols{j,3},inds(i)),2), nanmean(allMeasures(:,intCols{j,4},inds(i)),2)]);
    end
%     allInts = nancat(4, allInts, precInt); 
    yLims (i,:) = [min(nanmean(allInts(:,:,:,:,i),1),[],'all'),max(nanmean(allInts(:,:,:,:,i),1),[],'all')] .* [.9 1.1];
end
%%
intLabels = {'pre','post';...
    'low','high';...
    'low','high';};
intLegends = {'short','long';...
    'short','long';...
    'pre','post';};
intTitles = {'cue*delay','rew*delay','rew*cue'};

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
            title(factorLabels{4-j})
        end
        ylim(yLims(i,:))
    end
end

%%

save('./Data/ArrowMotivAnalysis.mat')

% end