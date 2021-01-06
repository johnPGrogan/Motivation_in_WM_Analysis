% ArrowFitTrajInterp
% fit models to points along response trajectories for all experiments

close all; clear; clc;

%% load trajectory data

if exist('./Data/ArrowFitTrajInterp.mat','file')
    load('./Data/ArrowFitTrajInterp.mat','allParsConds','o','nInterps','parNames');
else
    

    matFiles = {'ArrowMotivAnalysis.mat';
                'ArrowMotivCueBlockedAnalysis.mat';
                'ArrowMotivRetrocueAnalysis.mat';
                'ArrowMotivReallocateAnalysis.mat';
                'ArrowMotivMotorAnalysis_A.mat';
                'ArrowMotivMotorAnalysis_B.mat';
                };

    dataFolder = './Data/';

    nInterps = 100; % number of points to interpolate

    nExpts = 6;%length(matFiles);
    parfor i = 1:nExpts

        d = load(fullfile(dataFolder, matFiles{i}), 'dataFiles', 'dataStruct');

        o{i} = ArrowModelFitsTrajNormalise(d.dataFiles, d.dataStruct, nInterps, i==2);

    end

    
end

%%

% split expt 4 into a+b
o2 = o([1,2,3,4,4,5,6]); % copy, duplicate expt 4
eqConds = [1 0 0 1 1 0 0 1] == 1; % equal conds
fn = {'fits','pars','negLogLike','aic','bic','nTrialsEachType'}; % fields to split
for i = 1:length(fn)
    for j = 1:nInterps
        o2{4}(j).(fn{i}) = o{4}(j).(fn{i})(:,:,eqConds);
        o2{5}(j).(fn{i}) = o{4}(j).(fn{i})(:,:,~eqConds);
    end
end
[o2{4}.nTrialTypes] = deal(4);
[o2{5}.nTrialTypes] = deal(4);

%% do some plotting


condCols = {[1 1 1 1 2 2 2 2; 1 1 2 2 1 1 2 2; 1 2 1 2 1 2 1 2]; % rew, cue, delay
            [1 1 2 2;         1 2 1 2];                          % rew, cue
            [1 1 2 2;         1 2 1 2];                          % rew, congr
            [1 2 1 2; 1 1 2 2];                                  % rew, cue %             [1 1 2 2 1 1 2 2; 1 1 1 1 2 2 2 2; 1 2 2 1 1 2 2 1]; % rew, cue, equal
            [1 2 1 2; 1 1 2 2];                                  % rew, cue
            [1 2];                                               % rew
            [1 2];                                               % rew
            };
condCols2 = {   [1 1 1 1];
                [1 1];
                [1 1]; 
                [1 1]; %[1 3; 2 4];
                [1 1];
                [];
                [];
                };
            
condLabels = { {'1p','50p';     'pre','post';       'short','long'          };
               {'1p','50p';     'pre','post';                               };
               {'1p','50p';     'incongr','congr';                          };
               {'1p','50p';     'pre','post';};%      'equal', 'unequal'      };
               {'1p','50p';     'pre','post';};%      'equal', 'unequal'      };
               {'1p','50p'};
               {'1p','50p'};
               };
        
exptNums = {'1','2','3','4Eq','4Uneq','5a','5b'};

sepConds=0;
ybar = [20 1 0 0]; % for pbar

useClust=1; nPerms=5000;

nExpts = 5; % don't plot 5a & b

% nExpts = length(exptNums);
figure();
for j = 1:nExpts
   
    %[pp, time, cond, par]
    allPars = permute(nancat(4, o2{j}.pars),[1,4,3,2]);
    allPars(:,:,:,4) = 1 - allPars(:,:,:,1) - allPars(:,:,:,2);
    allPars = allPars(:,:,:,[3,4,1,2]);
    parNames = {'imprecision','P(target)','P(guess)','P(misbind)'};
    
    if sepConds

        parsConds = cat(5, allPars(:,:,condCols{j}(1,:)==1,:), allPars(:,:,condCols{j}(1,:)==2,:));
        if condCols2{j}
            parsConds2 = [];
            for i = 1:size(condCols2{j},1)
                parsConds2(:,:,i,:,:) = nanmean(parsConds(:,:,condCols2{j}(i,:),:,:),[3]);
            end
            parsConds = parsConds2;
        end
        parsConds = movmean(diff(parsConds,[],5),5,2); % rew eff, smoothed
    else
        parsConds = [];
        for k = 1:2
            parsConds(:,:,:,k) = nanmean(allPars(:,:,condCols{j}(1,:)==k,:),3);
        end
        parsConds = permute(parsConds,[1,2,4,3]);
    %     smooth
        parsConds = movmean(parsConds, 5, 2);
    end
    
    for iP = 1:4
        subplot(4, nExpts, (iP-1)*nExpts+j)
%         subplot(2,2,i);
        h = errorBarPlot(parsConds(:,:,:,iP),'area',1);

        if j==1; ylabel(parNames{iP}); end
        if iP==4; xlabel('interp(t)'); end

        if iP==1; title(['Expt ' exptNums{j}]);end
        
        % perm test
        if ~sepConds
            y1 = sq(diff(parsConds(:,:,:,iP),[],3));
    %         ybar = max(abs(ylim));
            [~,pp]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
            hold on;
            pbar(pp, 'yVal', ybar(iP), 'plotargs', {'Color', 'k', 'LineWidth',3});
        else
            ybar = [4 0.15 0.04 0.04]; % for pbar
            for k = 1:size(parsConds,3)
                [~,pp]=permutationOLS( parsConds(:,:,k,iP), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);

%             draw pbar
                hold on;
                pbar(pp, 'yVal', ybar(iP).*(1+k/20), 'plotargs', {'Color', h{k,1}.Color, 'LineWidth',3});
            end
            yline(0,':k');
        end
        
        xlim([0.5 nInterps + .5])
    end

%     legend([h{:,1}], condLabels{j}(1,:), 'Location','Best');
    allParsConds{j} = parsConds;
    
end

for iP=1:4
   makeSubplotScalesEqual(4, nExpts, iP*nExpts-(nExpts-1):iP*nExpts)
end

% legend
% h = findobj(gca,'Type','Line');
% legend(h([3,2,1]), [condLabels{j}(1,:), 'p < .05'], 'Location','Best');


%% do rmanova at each time point

figure(); clf
p = NaN(7,100,4,nExpts);

factorNames = { {'pp','rew','cue','delay'}; % names of factors in expts
                {'pp','rew','cue'};
                {'pp','rew','congr'};
                {'pp','rew','cue'};%,'equal'};
                {'pp','rew','cue'};%,'equal'};
                {'pp','rew'};
                {'pp','rew'};
                };
nTerms = [7;3;3;3;3;1;1;]; % number of terms from anovas
% will permute factors so reward always at top
factorPerms = { [1 2 6 5 4 3];  % pp time par rew cue delay
                [1 2 5 4 3];    % pp time par rew cue
                [1 2 5 4 3];    % pp time par rew congr
                [1 2 5 3 4];    % pp time par rew cue [1 2 6 4 5 3]
                [1 2 5 3 4];    % pp time par rew cue 
                [1 2 4 3];      % pp time par rew
                [1 2 4 3];      % pp time par rew
                };

for j = 1:nExpts
    allPars = permute(nancat(4, o2{j}.pars),[1,4,3,2]); % get this expt
    allPars(:,:,:,4) = 1 - allPars(:,:,:,1) - allPars(:,:,:,2);
    allPars = allPars(:,:,:,[3,4,1,2]);
    
    nFac = length(factorNames{j})-1; % number of factors for anova
    
    newDims = [size(allPars,1), 100, ones(1, nFac)*2, 4]; % dimensions to get
    allPars = reshape(allPars, newDims); % reshape into each factor
    allPars = permute(allPars, factorPerms{j}); % permute order
    
    catFacs = (1:nFac) + 1;
    for iPar = 1:4 % for each par
        for iT = 1:100 % for each timestep, do anova
            try
                stats = rmanova(sq(allPars(:,iT,iPar,:,:,:,:)), factorNames{j},'categorical',catFacs, 'DummyVarCoding','effects'); 
                p(1:nTerms(j),iT,iPar,j) = stats.pValue(2:end); % store pvals
            end
        end
        
        % plot pvals
        % make any > .05 = NaN
        p(p > .05) = NaN;
        subplot(4, nExpts, (iPar-1)*nExpts+j)
        imagep(p(:,:,iPar,j), stats.Term(2:end));
        colorbar off;
        if j==1; ylabel(parNames{iPar}); end
        if iPar==4; xlabel('interp(t)'); end

        if iPar==1; title(['Expt ' exptNums{j}]);end
   
    end
end


%% put all params/expts into one big table for regression, column per 
parNames2 = {'imprecision','target','guess','misbind'}; % wihout parentheses

allFactorNames = {'rew','cue','delay','congr'};%,'equal'};

allParsTab = table();
for j = 1:nExpts
    % get param across all expts
    allPars = permute(nancat(4, o2{j}.pars),[1,4,3,2]); % get this expt
    allPars(:,:,:,4) = 1 - allPars(:,:,:,1) - allPars(:,:,:,2);
    allPars = allPars(:,:,:,[3,4,1,2]);
    
    nFac = length(factorNames{j})-1; % number of factors for anova
    
    newDims = [size(allPars,1), 100, ones(1, nFac)*2, 4]; % dimensions to get
    allPars = reshape(allPars, newDims); % reshape into each factor
    allPars = permute(allPars, factorPerms{j}); %[pp time pars factors]
    
    
    % depivot to get factors, repeat to get all params as sep columns
    for iP = 1:4
        t1 = dePivot(sq(allPars(:,:,iP,:,:,:,:)));
        if iP==1
            t = t1;
        else
            t = [t t1(:,end)];
        end
    end
    t = [ones(size(t,1),1) * j, t]; % expt number
    t(:,2) = t(:,2) + (j*100); % make pp nums distinct
    
    % zscore
    t(:,4:end) = nanzscore(t(:,4:end));
    
    % make into table
    t = array2table(t, 'VariableNames', {'exptNum','pp','time',factorNames{j}{2:end},parNames2{:}});
    
    % add in columns of NaNs for missing factors
    missingFactors = allFactorNames(~ismember(allFactorNames, factorNames{j}));
    for i = 1:length(missingFactors)
        t.(missingFactors{i}) = NaN(height(t),1);
    end
    
    
    allParsTab = vertcat(allParsTab, t);
end

% % make factors -1, 1, 0
for iP = 1:length(allFactorNames)
    y = allParsTab.(allFactorNames{iP});
%     y(y==1) = -1;
%     y(y==2) = 1;
    y(isnan(y)) = 0;
    allParsTab.(allFactorNames{iP}) = y;
end


% expt5 = unequal
allParsTab.equal = nanzscore(allParsTab.exptNum ~= 5);
%% now regress all factors except time and reward, then regress across expts

timeResidsP = NaN(4,100,2);
timeResidsData = [];
timeResidsGroup = [];
for iP=1:4
    % now get p-values for each time point?
%     
    allParsTab.(['timeResids_' parNames2{iP}]) = NaN(height(allParsTab),1); % preallocate
    for j = 1:nExpts
        t = allParsTab(allParsTab.exptNum == j,:);

        fac = factorNames{j}(2:end);
        terms = flat([factorNames{j}(2:end); [repmat({'*'},1,length(factorNames{j})-2), {''}]])';
        
        formula = [parNames2{iP} ' ~ 1 + ' terms{:} ' + (1 | pp)']; 

        formula2 = [parNames2{iP} ' ~ 1 + ' terms{3:end} ' + (1 | pp)']; 


%       regress each par and expt, leaving rew+time resids
        timeResids{iP,j} = fitlme(t, formula2,'DummyVarCoding','effects');

        allParsTab.(['timeResids_' parNames2{iP}])(allParsTab.exptNum==j) = timeResids{iP,j}.residuals; % store this expts resids
    end
    
    % regress rew effect across expts per time point
    for iT=1:nInterps
        
        timeResids2 = fitlme(allParsTab(allParsTab.time==iT & allParsTab.exptNum ~= 5,:), ['timeResids_' parNames2{iP} '~ rew + (1 | pp)']);
        timeResidsP(iP,iT,1) = timeResids2.Coefficients.pValue(2:end);
        
        timeResids2 = fitlme(allParsTab(allParsTab.time==iT & allParsTab.exptNum == 5,:), ['timeResids_' parNames2{iP} '~ rew + (1 | pp)']);
        timeResidsP(iP,iT,2) = timeResids2.Coefficients.pValue(2:end);
    end
    
    % split by rew + time
    b = allParsTab(allParsTab.exptNum~=5,:); 
    timeResidsData1 = groupMeans(groupMeans(b.(['timeResids_' parNames2{iP}]),1,b.time,'dim'),2,groupMeans(round(b.rew),1,b.time,'dim'),'dim');
    timeResidsGroup1 = groupMeans(groupMeans(b.pp,1,b.time,'dim'),2,groupMeans(round(b.rew),1,b.time,'dim'),'dim'); % group for permOLS
    
    b2 = allParsTab(allParsTab.exptNum==5,:);  
    timeResidsData2 = groupMeans(groupMeans(b2.(['timeResids_' parNames2{iP}]),1,b2.time,'dim'),2,groupMeans(round(b2.rew),1,b2.time,'dim'),'dim');
    timeResidsGroup2 = groupMeans(groupMeans(b2.pp,1,b2.time,'dim'),2,groupMeans(round(b2.rew),1,b2.time,'dim'),'dim');
    
    timeResidsData = nancat(4, timeResidsData, nancat(5, timeResidsData1, timeResidsData2));
    timeResidsGroup = nancat(4, timeResidsGroup, nancat(5, timeResidsGroup1, timeResidsGroup2));
    
    % also do a rew*time regression
    f = fitlme(b, ['timeResids_' parNames2{iP} ' ~ 1 + rew*time + (1 | pp)'],'DummyVarCoding','effects');
    timeRewIntP(:,iP,1) = f.Coefficients.pValue(2:end);
    
    % for unequal
    f = fitlme(b2, ['timeResids_' parNames2{iP} ' ~ 1 + rew*time + (1 | pp)'],'DummyVarCoding','effects');
    timeRewIntP(:,iP,2) = f.Coefficients.pValue(2:end);
    
    % include equality as a factor, to see if rew*time differs between
    f = fitlme(allParsTab, ['timeResids_' parNames2{iP} ' ~ 1 + rew*time*equal + (1 | pp)'],'DummyVarCoding','effects');
    timeRewIntP2(:,iP,1) = f.Coefficients.pValue(2:end);
end

timeResidsData = permute(timeResidsData,[3,1,2,4,5]); %[tr time rew par expt]
timeResidsGroup = permute(timeResidsGroup,[3,1,2,4,5]);

%% do perm tests
doPerm=0;
if doPerm
    permArgs = {'cluster',0,'clustermethod','mean'};
    for iP=1:4
        for j=1:2
            Y = [timeResidsData(:,:,1,iP,j);timeResidsData(:,:,2,iP,j)]; % data
            X = [zeros(size(timeResidsData,1),1);ones(size(timeResidsData,1),1)]; % design [low; high]
            G = [timeResidsGroup(:,1,1,iP,j);timeResidsGroup(:,1,1,iP,j) ]; % group = pp
            [~,permP(iP,:,j)] = permutationOLS( Y , X , [] , G, permargs{:});
        end
    end
end
%% plot

exptLabels = {'equal','unequal'};

stars = p2stars(timeRewIntP);
for j = 1:2
    figure();
    cols = get(gca,'ColorOrder');
    cols = cols(4:5,:);
    
    for iP = 1:4
        subplot(2,2,iP)
        set(gca,'ColorOrder',cols);
        h(:,:,iP) = errorBarPlot(timeResidsData(:,:,:,iP,j),'area',1);
        ylabel([parNames(iP)]);
        xlabel('% of movement time');
        hold on;
        yline(0,':k');
        yl = [-1.2 1.2];
        ylim(yl);
        box off;
        if doPerm
            h1 = pbar(sq(permP(iP,:,j)),'yVal',-1,'alpha',.05); % p values per time
        else
            h1 = pbar(sq(timeResidsP(iP,:,j)),'yVal',min(yl),'alpha',.05); % p values per time
        end
        
        % rew*time int
        text(50 - length(stars{3,iP,j})*3., max(yl)*.8, stars(3,iP,j),  'fontsize', 20);
        
        xlim([1 100]);
        xticks([1 50 100]);
        xticklabels([0 50 100]);
    end

    subplot(2,2,1)
    legend([h{:,1},h1],{'1p','50p','p<.05'},'Location',[0.3172 0.7731 0.1411 0.1119]);
    makeSubplotScalesEqual(2,2);

%     SuperTitle(['residuals averaged by (' upper(exptLabels{j}) ') reward level']);

end

%% show reward effects over time

figure();
for iP=1:4
    subplot(2,2,iP)
    h = errorBarPlot(movmean(sq(diff(timeResidsData(:,:,:,iP,:),[],3)),5,2),'area',1); 
    yline(0,':k'); 
    ylabel(['reward effect on ' parNames{iP}]);
    xlabel('% of interpolated time');
end
makeSubplotScalesEqual(2,2);
legend([h{:,1}], {'equal','unequal'},'Location','Best');


%% save

save('./Data/ArrowFitTrajInterp.mat','allParsConds','o','o2','condLabels',...
    'nInterps','nExpts','parNames','sepConds','timeResidsData');
