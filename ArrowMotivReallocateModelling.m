% ArrowMotivReallocateModelling
% compare models, extract pars, analyse and plot pars

close all
clear
%%
load('./Data/ArrowMotivReallocateAnalysis.mat','dataStruct')

n = length(dataStruct);%num pps

%%
modelStruct = ArrowModelLoad(dataStruct); % get data for model fits

splitBy = 'trialTypes';
useMemFit = 0; % use MLE non MemFit

allSwap = ArrowModelCall(modelStruct, SwapModel(), splitBy,useMemFit); % fit the misbinding model
allMix = ArrowModelCall(modelStruct, StandardMixtureModel(), splitBy,useMemFit); % fit the model without misbinding

allVarPrec = ArrowModelCall(modelStruct, VariablePrecisionModel(), splitBy,useMemFit);

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
nPars = length(bestModel.model.paramNames);
g = sq(bestPars(:,1,:));
sd = sq(bestPars(:,3,:));

nPars = bestModel.nPars + 1;
if nPars == 4
    b = sq(bestPars(:,2,:));%beta - only in swap
    a = 1 - g - b;
elseif nPars == 5
    b = sq(bestPars(:,2,:));%beta - only in swap
    pr = sq(bestPars(:,4,:));
    a = 1 - g - b - pr;
else
    b = NaN(n,8);
    a = 1 - g;
end

bic = sq([bestModel.bic]);

allPars = nancat(3, sd, a, g);%conc
if nPars==4
    allPars = nancat(3, allPars,b);
elseif nPars == 5
    allPars = nancat(3, allPars, b, pr);
end

%% stats
nTT = size(allPars,2);
x = zeros(n,1);%zeros
motiv = [x+1;x;x;x+1;x+1;x;x;x+1];%1 = 
rew = [x;x;x+1;x+1;x;x;x+1;x+1];%1=low, 2=hi
cue = [x;x;x;x;x+1;x+1;x+1;x+1];%1=pre, 2=post
ppInd = repmat([1:n]',nTT,1);%participant indices 1:n repeated for each condition
modelMat = [0 0 0 1; 1 0 0 0;0 1 0 0; 0 0 1 0;...
    1 1 0 0; 1 0 1 0; 0 1 1 0;...
    1 1 1 0];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors = {rew, cue, motiv, ppInd};%combine factors into cells - IVs and participant
factorLabels = {'rewardLevel','cueTime','isEqual','pp'};%labels for factors

allParsRM = allPars(:,[2 3 1 4 6 7 5 8],:);
pVals = [];
for i = 1:nPars
%     DVVec =  col(allPars(:,:,i));%make a matrix into a vector
%     [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',4,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms

    rmStats{i} = rmanova(reshape(allParsRM(:,:,i), [],2,2,2), {'pp','reward','equal','cue'});
    pVals(:,i) = rmStats{i}.pValue;
end


%%
% pVals = [stats.p];
alphas = [.05, .01, .001, .0001];
nStars = zeros(8,8);
stars = cell(7,8);
for i = 1:size(pVals,2)
    nStars(:,i) = sum(pVals(:,i) < alphas,2);
    for j = 1:7
        stars{j,i} = repmat('*',1,nStars(j+1,i));
    end
end

%%
parNames = {'imprecision','target','guess','misbind','m180'};

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

%% also do stats on un/equal separately
%% anova
allParsConds = nancat(4, allPars(:,[1 4 5 8],:), allPars(:,[2 3 6 7],:));

x = zeros(n,1);%zeros
rew = [x;x+1;x;x+1];%1=low, 2=hi
cue = [x;x;x+1;x+1];%1=pre, 2=post
ppInd = repmat([1:n]',nTT./2,1);%participant indices 1:n repeated for each condition
modelMat2 = [0 0 1; 1 0 0;0 1 0; ...
    1 1 0 ; ];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors2 = {rew, cue, ppInd};%combine factors into cells - IVs and participant
factorLabels2 = {'rewardLevel','cueTime','pp'};%labels for factors

for j = 1:2
    for i = 1:nPars
%         DVVec =  reshape(allParsConds(:,:,i,j),1,nTT*n/2)';%make a matrix into a vector
%         [stats2(i,j).p,stats2(i,j).tab,stats2(i,j).s] = anovan(DVVec,factors2,'varnames',factorLabels2,'display','off','model','full','random',3,'model',modelMat2);
        rmStats2{i,j} = rmanova(reshape(sq(allParsConds(:,:,i,j)),[],2,2), {'pp','reward','cue'});
        pVals2(:,i,j) = rmStats2{i,j}.pValue;

    end
end

%% combine stats into one table and display

% pVals = nancat(3, [stats2(:,1).p],[stats2(:,2).p]);
alphas = [.05, .01, .001, .0001];

nStars2 = zeros(3,nPars);
stars2 = cell(3,nPars);
for k = 1:2
    for i = 1:nPars
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
    for i = 1:nPars
        allStats1.cond = repmat(condLabels(1,j),3,1);
        allStats1.measure = repmat(parNames(i),3,1);
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

%% figs

figure()
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 4 2 3; 5 8 6 7];
for i = 1:nPars
    if nPars>1
        subplot(ceil(nPars/2),2,i)
    end

    h = errorBarPlot(nancat(3, allPars(:,condInds(1,:),i), ...
                    allPars(:,condInds(2,:),i)),'type','line');
    hold on
    for j = 1:length(h)
        h(j).LineWidth = 2;
        h(j).Marker = markers{j};
        h(j).Color = lineColours(j,:);
        if size(condInds,2) == 4
            h(j).LineStyle = 'none';
            plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:));
            plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:));
        end
    end

    hold on
    for j = 1:2
        h(j).LineWidth = 2;
        h(j).Marker = markers{j};
        h(j).Color = lineColours(j,:);
    end
    
    set(gca,'XTick',1:4,'XTickLabel',{'Low','High','Low','High'})
    xlabel('Equal                Unequal')
    xlim([0.5 4.5])
    ylabel(parNames{i})
    if i==1
        legend(fliplr(h), fliplr({'Pre','Post'}),'Location','Best')
    end
    box off


end


%%
allMeanConds = NaN(n,2,3,nPars);
cols = {[1 2 5 6]; [3 4 7 8];... % low vs high
    [1 2 3 4]; [5 6 7 8];... % pre vs post
    [1 4 5 8];[2 3 6 7];}; % motiv vs reallocate
for i = 1:nPars%for each DV
    for j = 1:3%for each factor
        %get means
        allMeanConds(:,:,j,i) = [nanmean(allPars(:,cols{j*2-1},i),2), nanmean(allPars(:,cols{j*2},i),2)];
    end
    
    %get ylims for each measure
    yLims (i,:) = [min(nanmean(allMeanConds(:,:,:,i),1),[],'all'),max(nanmean(allMeanConds(:,:,:,i),1),[],'all')] .* [.9 1.1];
end




%% plot main effects
xFacLabels = {'low','high';...
    'pre','post';...
    'motiv','reallo'};
figure()

for i = 1:nPars
    for j = 1:3
        subplot(nPars,3,(i-1)*3+j)
        errorBarPlot(allMeanConds(:,:,j,i),'type','bar');
        set(gca,'XTick',1:2,'XTickLabel',xFacLabels(j,:))
        if i==1
            title(factorLabels{j})
        end
        if j==1
            ylabel(parNames{i})
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
for i = 1:nPars
    for j = 1:3
        allInts(:,:,:,j,i) = nancat(3, [nanmean(allPars(:,intCols{j,1},i),2),nanmean(allPars(:,intCols{j,2},i),2)],...
            [nanmean(allPars(:,intCols{j,3},i),2),nanmean(allPars(:,intCols{j,4},i),2)]);
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

intColours = {[.1 .8 .2], [.85, .32, .98];...
    [1 .8 .2], [0 .5 .8];...
    [1 .3 0], [.6, .2, .6];};
figure()

for i = 1:nPars
    
    for j = 1:3
        subplot(nPars,3,(i-1)*3+j)
%         h = errorBarPlot(permute(allInts(:,:,:,j,i),[1,3,2]), 'type', 'line','plotargs',{'LineWidth',2});
        h = errorBarPlot(allInts(:,:,:,j,i), 'type', 'line','plotargs',{'LineWidth',2});
%         h(1).Color = intColours{j,1};
%         h(2).Color = intColours{,2};
        if i==nPars
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
        if i==nPars && j==3
            legend(intLabels(1,:),'Location','Best')
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
        
%         yLims = repmat( max(abs([min(effects,[],'all'), max(effects,[],'all')])),1,2) .* [-1.1 1.1];%get ylimits
        subplot(nPars,3,(i-1)*3+j)
        bar(sort(effects(:,:,j,i),'ascend'));
        if j==1
            ylabel(['\Delta ' parNames{i}])
        end
        if i==1
            title(factorLabels{j})
        end
        ylim(yLims(i,:))
    end
end



%%

save('./Data/ArrowMotivReallocateModelling.mat')