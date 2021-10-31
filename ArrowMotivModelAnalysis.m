% ArrowsModelAnalysis
% load up ArrowModelFits.mat
% compare models, extract pars, analyse and plot pars

close all
clear
%%

load('./Data/ArrowMotivAnalysis.mat','dataStruct')

n = length(dataStruct);%num pps

modelStruct = ArrowModelLoad(dataStruct); % get data for model fits

splitBy = 'trialTypes'; % fit each condition separately
useMemFit = 0; % use MLE non MemFit

allSwap = ArrowModelCall(modelStruct, SwapModel(), splitBy,useMemFit); % fit the misbinding model
allMix = ArrowModelCall(modelStruct, StandardMixtureModel(), splitBy,useMemFit); % fit the model without misbinding


%%
do180 = 0; % fit model with sep misbinding to 180degrees from target?

if do180
    % set prevResp to 180 (i.e. 180 from target)
    for i = 1:n
        modelStruct(i).data.prevResp = ones(size(modelStruct(i).data.errors))*180;%mod(modelStruct(i).data.targets + 360,360)-180;    
    end
    % fit with extra misbinding
    all180 = ArrowModelCall(modelStruct, SwapModelPrev, splitBy,useMemFit);

end

%% compare overall fits (across all conditions)
models = {SwapModel(), StandardMixtureModel()};
if do180
    models{end+1} = SwapModelPrev();
end

for i=1:length(models)
    allTrialsFit(i) = ArrowModelCall(modelStruct, models{i}, 'allTrials',0); % fit the misbinding model
end
meanBIC = nanmean([allTrialsFit.bic])
[m,i] = min(meanBIC)


%% model comparison

mean(sq([allSwap.bic] < [allMix.bic]))

bestModel = allSwap;
bestModel.model = SwapModel();

%% pars

bestPars = [bestModel.pars];%8 x nPars x n
g = sq(bestPars(:,1,:)); % guessing
sd = sq(bestPars(:,3,:)); % imprecisoin
 
nPars = bestModel.nPars + 1; % add target param
if nPars == 4 % if misbinding model
    b = sq(bestPars(:,2,:));%beta - only in swap
    a = 1 - g - b;
elseif nPars == 5 % if another param
    b = sq(bestPars(:,2,:));%beta - only in swap
    pr = sq(bestPars(:,4,:)); % get extra
    a = 1 - g - b - pr;
else % if mixture model
    b = NaN(n,8);
    a = 1 - g;
end

bic = sq([bestModel.bic]);

%% stats

x = zeros(n,1);
del = [x;x+1;x;x+1;x;x+1;x;x+1];%set to 1 when reward was high, 0 for low
cue = [x;x;x+1;x+1;x;x;x+1;x+1];%pre/post cue
rew = [x;x;x;x;x+1;x+1;x+1;x+1];%short/long delay

ppInd = repmat([1:n]',8,1);
modelMat = [0 0 0 1; 1 0 0 0;0 1 0 0; 0 0 1 0;...
    1 1 0 0; 1 0 1 0; 0 1 1 0;...
    1 1 1 0];
factors = {rew, cue, del, ppInd};
factorLabels = {'rewardLevel','cueTime','delayLength','pp'};

allPars = nancat(3, sd, a, g);%combine
if nPars==4
    allPars = nancat(3, allPars,b);
elseif nPars == 5
    allPars = nancat(3, allPars, b, pr);
end

for i = 1:nPars
%     DVVec =  reshape(allPars(:,:,i),1,8*n)';%make a matrix into a vector
%     [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',4,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms

    rmStats{i} = rmanova(reshape(allPars(:,:,i), [],2,2,2), {'pp','delay','cue','reward'});
    pVals(:,i) = rmStats{i}.pValue(2:end);
end


%%
% pVals = [stats.p];
alphas = [.05, .01, .001, .0001];
for i = 1:size(pVals,2)
    nStars(:,i) = sum(pVals(:,i) < alphas,2);
    for j = 1:7
        stars{j,i} = repmat('*',1,nStars(j,i));
    end
end

%%
parNames = {'SD','target','guess','misbind','180'};

allStats = table();
allStats1 = table();

for i = 1:nPars
    allStats1.measure = repmat(parNames(i),7,1);
    allStats1.effect = rmStats{i}.Term(2:end);
    allStats1.df = rmStats{i}.DF1(2:end);
    allStats1.dfErr = rmStats{i}.DF2(2:end);
    allStats1.F = rmStats{i}.FStat(2:end);
    allStats1.p = rmStats{i}.pValue(2:end);
    allStats1.stars = stars(:,i);
    
    
    allStats = vertcat(allStats, allStats1);
end  
    

disp(allStats(:,[1,2,6,7]))%,11,12]))


%% figs
yLabs = {'imprecision (SD)', 'proportion targets', 'proportion guesses', 'proportion misbinds', 'proportion 180'};
figure()
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 5 3 7; 2 6 4 8];
for i = 1:nPars
    subplot(nPars-2,2,i)
    
    h = errorBarPlot(nancat(3, allPars(:,condInds(1,:),i), ...
                    allPars(:,condInds(2,:),i)),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    hold on
    JoinUpErrorPoints(h, [1 2; 3 4]);
    for j = 1:length(h)
        h(j).Marker = markers{j};
    end


    set(gca,'XTick',1:4,'XTickLabel',{'Low','High','Low','High'})
    xlabel('Pre                Post')
    xlim([0.5 4.5])
    ylabel(yLabs{i})
    if i==1
        legend(fliplr(h), fliplr({'Short','Long'}),'Location','North')
    end

    box off

end
%%

allMeanConds = NaN(n,2,3,nPars);
cols = {[1 3 5 7];[2 4 6 8];...
    [1 2 5 6]; [3 4 7 8];...
    [1 2 3 4]; [5 6 7 8];};
for i = 1:nPars%for each DV
    for j = 1:3%for each factor
        %get means
        allMeanConds(:,:,j,i) = [nanmean(allPars(:,cols{j*2-1},i),2), nanmean(allPars(:,cols{j*2},i),2)];
    end
    
    %get ylims for each measure
    yLims (i,:) = [min(nanmean(allMeanConds(:,:,:,i),1),[],'all'),max(nanmean(allMeanConds(:,:,:,i),1),[],'all')] .* [.9 1.1];
end



%% plot main effects
xFacLabels = {'short','long';...
    'pre','post';...
    'low','high';};
figure()

for i = 1:nPars
    for j = 1:3
        subplot(nPars,3,(i-1)*3+j)
        errorBarPlot(allMeanConds(:,:,j,i),'type','bar');
        set(gca,'XTick',1:2,'XTickLabel',xFacLabels(j,:))
        if i==1
            title(factorLabels{4-j})
        end
        if j==1
            ylabel(parNames{i})
        end
        ylim(yLims(i,:))
        text(0.5, 0.8, stars(j,i),'Units','normalized','FontSize',14)

    end
end


%% plot interaction of cue*delay

intCols =  {   [1 5], [2 6], [3 7], [4 8];%cue*delay
        [1 3], [2 4], [5 7], [6 8];...%rew*delay
        [1 2],[3 4],[5 6],[7 8];};%rew*cue

allInts = [];
for i = 1:nPars
    for j = 1:3
        allInts(:,:,:,j,i) = nancat(3, [nanmean(allPars(:,intCols{j,1},i),2),nanmean(allPars(:,intCols{j,2},i),2)],[nanmean(allPars(:,intCols{j,3},i),2),nanmean(allPars(:,intCols{j,4},i),2)]);
    end
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
intCols = {[1 .3 0], [.6, .2, .6];...
    [1 .8 .2], [0 .5 .8];...
    [.1 .8 .2], [.85, .32, .98]};
figure()


for i = 1:nPars
    
    for j = 1:3
        subplot(nPars,3,(i-1)*3+j)
        h = errorBarPlot(allInts(:,:,:,j,i), 'type', 'line','plotargs',{'LineWidth',2});
        h(1).Color = intCols{j,1};
        h(2).Color = intCols{j,2};
        if i==nPars
            set(gca,'XTick',1:2,'XTickLabel',intLabels(j,:))
        else
            set(gca,'XTickLabel','')
        end
        xlim([.5 2.5])
        ylim([yLims(i,:)])
        text(0.5, 0.8, stars(j+3,i),'Units','normalized','FontSize',14)
        if i==1
            title(intTitles{j})
        end
        if i==4
            legend(intLegends(j,:),'Location','Best')
        end
        if j==1
            ylabel(parNames{i})
        end
        
    end
end

%% effect size for pps

figure()

effects = allMeanConds(:,2,:,:) - allMeanConds(:,1,:,:);%high error - low error
%

for i = 1:nPars
    
   yLims(i,:) = repmat( max(abs([min(effects(:,:,:,i),[],'all'), max(effects(:,:,:,i),[],'all')])),1,2) .* [-1.1 1.1];%get ylimits 

    for j = 1:3
        
        subplot(nPars,3,(i-1)*3+j)
        bar(sort(effects(:,:,j,i),'ascend'));
        if j==1
            ylabel(['\Delta ' parNames{i}])
        end
        if i==1
            title(factorLabels{4-j})
        end
        ylim(yLims(i,:))
    end
end

%%

save('./Data/ArrowMotivModelAnalysis.mat')