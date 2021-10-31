% ArrowMotivRetrocueModelling

close all
clear
%%
load('./Data/ArrowMotivRetrocueAnalysis.mat','dataStruct','condLabels')

n = length(dataStruct);%num pps

%% remove invalid/incorrect trials
% find which trials you want to remove
% set those to NaN

% remove all trials where the retrocue response was incorrect
for i = 1:n % for each person
    toRemove = find(dataStruct(i).RetrocueResponseCorrect == 0); % trials where retrocue answered wrong
    if toRemove
        for j = 1:length(toRemove)
            dataStruct(i).result.data(toRemove(j)) = structfun(@(x) NaN(size(x)),...
                dataStruct(i).result.data(toRemove(j)),'UniformOutput',0); % set all info for this trial to NaN
        end
    end
end

%%
modelStruct = ArrowModelLoad(dataStruct); % get data for model fits

%%
splitBy = 'trialTypes';
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
    b = NaN(n,4);
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

x = zeros(size(allPars(:,:,1)));%zeros
rew = x;
rew(:,3:4) = 1;
cueCongr = x;
cueCongr(:,[2 4]) = 1;

ppInd = repmat([1:n]',1,4);%participant indices 1:n repeated for each condition
modelMat = [0 0 1; 1 0 0;0 1 0;...
    1 1 0;];%matrix of effects/interactions to test. has all simple effects, and all interactions (except those with ppInd)
factors = {col(rew), col(cueCongr), col(ppInd)};%combine factors into cells - IVs and participant
factorLabels = {'rewardLevel','cueCongr','pp'};%labels for factors


for i = 1:nPars
    DVVec =  col(allPars(:,:,i));%make a matrix into a vector
    [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',3,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms
    
    rmModelStats{i} = rmanova(reshape(allPars(:,:,i), n,2,2), {'pp','cong','rew'});
    pVals(:,i) = rmModelStats{i}.pValue;

end



%%
% pVals = [stats.p];
alphas = [.05, .01, .001, .0001];
nStars = zeros(4,nPars);
stars = cell(3,nPars);
for i = 1:size(pVals,2)
    nStars(:,i) = sum(pVals(:,i) < alphas,2);
    for j = 1:3
        stars{j,i} = repmat('*',1,nStars(j+1,i));
    end
end

%%
parNames = {'imprecision','target','guess','misbind','m180'};

allStats = table();
allStats1 = table();

for i = 1:nPars
    allStats1.measure = repmat(parNames(i),3,1);
    allStats1.effect = rmModelStats{i}.Term(2:end);
    allStats1.df = rmModelStats{i}.DF1(2:end);
    allStats1.dfErr = rmModelStats{i}.DF2(2:end);
    allStats1.F = rmModelStats{i}.FStat(2:end);
    allStats1.p = rmModelStats{i}.pValue(2:end);
    allStats1.stars = stars(:,i);
    
    
    allStats = vertcat(allStats, allStats1);
end  
    

disp(allStats(:,[1,2,6,7]))%,11,12]))


%% figs

figure()
markers = {'o','^'};
lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
condInds = [1 3; 2 4];
for i = 1:nPars
    subplot(ceil(nPars/2),2,i)
    
%     for k = 1:2
        h = errorBarPlot(nancat(3, allPars(:,condInds(1,:),i), ...
            allPars(:,condInds(2,:),i)),'type','line');%,'xaxisvalues',[k*2-1;k*2]);
%         hold on
        for j = 1:2
            h(j).LineWidth = 2;
            h(j).Marker = markers{j};
            h(j).Color = lineColours(j,:);
        end
%     end
    set(gca,'XTick',1:2,'XTickLabel',condLabels(2,:))
    xlabel('reward')
    xlim([0.5 2.5])
    ylabel(parNames{i})
    if i==1
        legend(condLabels(1,:), 'Location','Best')
    end

    box off

end


%%

save('./Data/ArrowMotivRetrocueModelling.mat')