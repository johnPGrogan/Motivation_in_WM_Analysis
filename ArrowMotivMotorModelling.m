function ArrowMotivMotorModelling(exptVersion)

%%

if ~(strcmp(exptVersion,'A') || strcmp(exptVersion,'B'))
    error('exptVersion must be A or B')
end
load(sprintf('./Data/ArrowMotivMotorAnalysis_%c.mat', exptVersion), 'dataStruct', 'allMeasures')

n = length(dataStruct);%num pps

%%
modelStruct = ArrowModelLoad(dataStruct); % get data for model fits

if exptVersion=='A'
    nTrialsToSkip = sum(isnan(allMeasures(1,:,1))); % skip trials 1:2 if they are set to NaN
    for i = 1:length(modelStruct)
        modelStruct(i).data.errors(1:nTrialsToSkip) = NaN;
    end
else
    modelStruct(6).data.errors(26) = NaN; % outlying error - skip this
end

splitBy = 'trialTypes';
useMemFit = 0; % use MLE non MemFit

allSwap = ArrowModelCall(modelStruct, SwapModel(), splitBy,useMemFit); % fit the misbinding model
allMix = ArrowModelCall(modelStruct, StandardMixtureModel(), splitBy,useMemFit); % fit the model without misbinding
allNoGuess = ArrowModelCall(modelStruct, NoGuessingModel(), splitBy,useMemFit); % fit the model without misbinding

%% model comparison

nanmean(cat(2, allSwap.bic, allMix.bic, allNoGuess.bic))

% the noguessing model will probably be best, due to low guessing and
% misbinding, but this just means those params will be close to zero in
% misbinding model
bestModel = allSwap;

%% pars

bestPars = [bestModel.pars];%8 x nPars x n
g = sq(bestPars(:,1,:));
sd = sq(bestPars(:,end,:));

nPars = bestModel.nPars;
if nPars>1; nPars = nPars+1; end
if nPars == 4
    b = sq(bestPars(:,2,:));%beta - only in swap
    a = 1 - g - b;
elseif nPars == 5
    b = sq(bestPars(:,2,:));%beta - only in swap
    pr = sq(bestPars(:,4,:));
    a = 1 - g - b - pr;
elseif nPars==1
    a = ones(n,2);
    g = NaN(n,2);
else
    b = NaN(n,2);
    a = 1 - g;
end

bic = sq([bestModel.bic]);

allPars = nancat(3, sd, a, g);%conc
if nPars==4
    allPars = nancat(3, allPars,b);
elseif nPars == 5
    allPars = nancat(3, allPars, b, pr);
elseif nPars == 1
    allPars = allPars(:,:,1);
end


%% stats

x = ones(n,2,12);%zeros
rew = x;
rew(:,1,:) = 0;
rewTrial = col(rew);

ppInd = repmat([1:n]', 1, 2);
ppIndTrial = col(ppInd);
% ppInd = repmat([1:n]',nTT,1);%participant indices 1:n repeated for each condition
modelMat = [0 1; 1 0 ;]; % random pp + reward main

factorsTrial = {rewTrial,  ppIndTrial};%combine factors into cells - IVs and participant
factors = {col(nanmean(rew,3)), col(nanmean(ppInd,3))}; % average across trials
factorLabels = {'rewardLevel','pp'};%labels for factors

for i = 1:nPars
    DVVec =  col(allPars(:,:,i));%make a matrix into a vector
    [stats(i).p,stats(i).tab,stats(i).s] = anovan(DVVec,factors,'varnames',factorLabels,'display','off','model','full','random',2,'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms
    
    if max(diff(allPars(:,:,i)),[],'all') < .00001
        pVals(:,i) = 1;
    else
        rmModelStats{i} = rmanova(allPars(:,:,i), {'pp','rew'});
        pVals(:,i) = rmModelStats{i}.pValue;
    end
end



%%
% pVals = [stats.p];
alphas = [.05, .01, .001, .0001];
nStars = zeros(2,nPars);
stars = cell(1,nPars);
for i = 1:size(pVals,2)
    nStars(:,i) = sum(pVals(:,i) < alphas,2);
    for j = 1
        stars{j,i} = repmat('*',1,nStars(j+1,i));
    end
end

%%
parNames = {'imprecision','target','guess','misbind','pr'};

allStats = table();
allStats1 = table();

for i = 1:nPars
    allStats1.measure = repmat(parNames(i),1,1);
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
condInds = [1 2];
for i = 1:nPars
    if nPars>1
        subplot(ceil(nPars/2),2,i)
    end
    
    for j = 1
        h(j) = errorBarPlot(allPars(:,:,i),'type','line',...
            'type','line','plotargs',{'LineWidth',2.5,'Marker',markers{j}});
        hold on
%         h(j).XData(3-j) = NaN;    
    end
%         hold on
%         for j = 1:2
%             h(j).LineWidth = 2;
%             h(j).Marker = markers{j};
%             h(j).Color = lineColours(j,:);
%         end
%     end
    set(gca,'XTick',1:2,'XTickLabel',{'Low','High'})
    xlabel('Reward')
    xlim([0.5 2.5])
    ylabel(parNames{i})
%     if i==1
%         legend({'low','high'},'Location','Best')
%     end

    

end


%%

save(sprintf('./Data/ArrowMotivMotorModelling_%c.mat', exptVersion));