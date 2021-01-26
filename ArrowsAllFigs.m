function ArrowsAllFigs
% Draw figs for WM Motiv paper
% Requires you have run all main analyses and saved outputs (e.g. files
% ending in *Analysis.m, *Modelling.m, and ArrowMotivMotorCorrels.m)
% It will draw subplots for each experiment (combining 5a and 5b), and then
% the motor-wm correlations
% Optional saving

close all; clc; % clear
set(0,'DefaultAxesFontSize',16)

measureNames = {'abs angular error (deg)', 'initiation RT (ms)',...
    'completion RT (ms)', 'total RT (ms)', 'imprecision', 'P(target)',...
    'P(guess)', 'P(misbind)'};

% legOpts = {'Location', [0.4536 0.4826 0.1161 0.0380], 'Orientation', 'Horizontal'};
legOpts = {'Location', 'North'};

% set line spec
markers = {'o','^'};
lineStyles = {'-','--','-.'};
lineColours = [0 .447 .741; .85 .325 .098; .7 0 .7; 0 .7 .2; .929 .694 .125; .635 .078 .184];

% % % each cell in array is a subplot() call. first cell is the line plots,
% % % then trajectory plot, then param scatter plot, then histograms. empty = ignore
% subInds = {[2 3 1; 2 3 2; 2 3 4]; [2 3 3]; [2 3 5]; [4 3 9; 4 3 12]}; % shape, inds for: myFig, traj, cross, rew eff
% subInds = {[2 3 1; 2 3 2; 2 3 4]; [2 3 3]; [2 3 5]; [2 3 6]}; % shape, inds for: myFig, traj, cross, rew eff
subInds = {[2 3 1; 2 3 2; 2 3 4; 2 3 5; 2 3 6]; [2 3 3]; []; []}; % shape, inds for: myFig, traj, cross, rew eff

makeSame = {[5 6]}; % make these subplot inds the same scales

saveNames = cellfun(@(x) sprintf('./Figs/WmMotivExpt%s.svg',x), {'1','2','3','4','5a','5b'},'UniformOutput',0)';
saveNames{7} = 'WmMotivCorrels.svg';
saveNames{8} = 'WmMotivModelTraj.svg';
saveNames = cell(8,1); % uncomment to prevent saving


%% expt 1

% load data for plotting
load('./Data/ArrowMotivAnalysis.mat', 'allMeasures');
load('./Data/ArrowMotivModelAnalysis.mat', 'allPars');

condInds = [1 5 4 8; 3 7 2 6 ]; % [low reward, high reward]

% get these into dimensions for plotting: [pp xticks legend subplot]
allMeasures = permute(nancat(4, allMeasures(:,condInds(1,:),:), allMeasures(:,condInds(2,:),:)),[1,2,4,3]);
allPars = permute(nancat(4, allPars(:,condInds(1,:),:), allPars(:,condInds(2,:),:)),[1,2,4,3]);

allMeasures = nancat(4, allMeasures, allPars); % combine beh + pars

% factor names
xTickNames = {'1p', '50p', '1p', '50p'};
xLabel = {'Short', '                  ', 'Long'};
leg = {'1p',' 50p','Pre', 'Post', 'Short', 'Long'};

% trim and reorder
allMeasures = allMeasures(:,:,:,[1, 4, 5, 7, 8]); % 4th dim = [error, endRT, SD, g, B]

% load trajectories during responses
traj = load('./Data/ArrowTrajAll.mat');
iTraj = 1; % expt num

modelTraj = load('./Data/ArrowFitTrajInterp.mat'); % or [] if not plotting these

% call plotting function
mySubPlots(allMeasures, subInds, markers, lineColours, measureNames([1, 4, 5, 7, 8]),...
    xTickNames, xLabel, leg, legOpts, saveNames{1}, traj, iTraj,1, lineStyles, modelTraj, makeSame);


%% expt 2

% load data for plotting
load('./Data/ArrowMotivCueBlockedAnalysis.mat','allMeasures');
load('./Data/ArrowMotivCueBlockedModelling.mat', 'allPars');

condInds = [1 3; 2 4];

% [pp xticks legend subplot]
allMeasures = permute(nancat(4, allMeasures(:,condInds(1,:),:), allMeasures(:,condInds(2,:),:)),[1,2,4,3]);
allPars = permute(nancat(4, allPars(:,condInds(1,:),:), allPars(:,condInds(2,:),:)),[1,2,4,3]);

allMeasures = nancat(4, allMeasures, allPars); % combine beh + pars

% trim and reorder
allMeasures = allMeasures(:,:,:,[1, 4, 5, 7, 8]);

% factor names
xTickNames = {'1p', '50p'};
xLabel = {'Reward'};
leg = {'1p', '50p','Pre', 'Post'};

mySubPlots(allMeasures, subInds, markers, lineColours, measureNames([1, 4, 5, 7, 8]),...
    xTickNames, xLabel, leg, legOpts, saveNames{2}, traj, 2 , 1, lineStyles, modelTraj, makeSame);
%% expt 3

% load data for plotting
load('./Data/ArrowMotivRetrocueAnalysis.mat','allMeasures');
load('./Data/ArrowMotivRetrocueModelling.mat', 'allPars');

condInds = [1 3; 2 4];

% [pp xticks legend subplot]
allMeasures = permute(nancat(4, allMeasures(:,condInds(1,:),:), allMeasures(:,condInds(2,:),:)),[1,2,4,3]);
allPars = permute(nancat(4, allPars(:,condInds(1,:),:), allPars(:,condInds(2,:),:)),[1,2,4,3]);

allMeasures = nancat(4, allMeasures(:,:,:,1:4), allPars); % combine beh + pars

% trim and reorder
allMeasures = allMeasures(:,:,:,[1, 4, 5, 7, 8]);

% factor names
xTickNames = {'1p', '50p'};
xLabel = {'Reward'};
leg = {'1p', '50p','Incongruent', 'Congruent'};
cols = lineColours([5,6,3,4],:);

mySubPlots(allMeasures, subInds, markers, cols, measureNames([1, 4, 5, 7, 8]),...
    xTickNames, xLabel, leg, legOpts, saveNames{3}, traj, 3, 1, lineStyles, modelTraj, makeSame);

%% expt 4

% load data for plotting
load('./Data/ArrowMotivReallocateAnalysis.mat','allMeasures');
load('./Data/ArrowMotivReallocateModelling.mat', 'allPars');

condInds = [1 4 2 3; 5 8 6 7];

% [pp xticks legend subplot]
allMeasures = permute(nancat(4, allMeasures(:,condInds(1,:),:), allMeasures(:,condInds(2,:),:)),[1,2,4,3]);
allPars = permute(nancat(4, allPars(:,condInds(1,:),:), allPars(:,condInds(2,:),:)),[1,2,4,3]);

allMeasures = nancat(4, allMeasures, allPars); % combine beh + pars

% trim and reorder
allMeasures = allMeasures(:,:,:,[1, 4, 5, 7, 8]);

% factor names
xTickNames = {'1p', '50p', '1p', '50p'};
xLabel = {'Equal','              ','Unqeual'};
leg = {'1p','50p','Pre', 'Post','Equal','Unequal'};

mySubPlots(allMeasures, subInds, markers, lineColours, measureNames([1, 4, 5, 7, 8]),...
    xTickNames, xLabel, leg, legOpts, saveNames{4}, traj, 4, 1, lineStyles([1 3]), modelTraj, makeSame);

%% expt 5a

% load data for plotting
load('./Data/ArrowMotivMotorAnalysis_A.mat','allRewMeasures');
load('./Data/ArrowMotivMotorModelling_A.mat', 'allPars');

%%%%% maybe exclude some trials?

% [pp rew measure]
allMeasures = nanmean(allRewMeasures,3);
allPars = permute(allPars, [1 2 4 3]);
allMeasures = nancat(4, allMeasures, allPars); % combine beh + pars

% trim
allMeasures = allMeasures(:,:,:,[1,4]);

subInds2 = {[2 3 1; 2 3 2]; [2 3 3]; []; []}; % shape, inds for: myFig, traj, cross, rew eff

% factor names
xTickNames = {'1p', '50p'};
xLabel = {'Reward'};

mySubPlots(allMeasures, subInds2, markers, lineColours, measureNames([1, 4]),...
    xTickNames, xLabel, [], [], [], traj, 5, 1, lineStyles , modelTraj);


%% expt 5b

% load data for plotting
load('./Data/ArrowMotivMotorAnalysis_B.mat','allRewMeasures');
load('./Data/ArrowMotivMotorModelling_B.mat', 'allPars');

% [pp rew measure]
allMeasures = nanmean(allRewMeasures,3);
allPars = permute(allPars, [1 2 4 3]);
allMeasures = nancat(4, allMeasures, allPars); % combine beh + pars

% trim
allMeasures = allMeasures(:,:,:,1:4);

subInds2 = {[2 3 4; 2 3 5]; [2 3 6]; []; []}; % shape, inds for: myFig, traj, cross, rew eff

makeSame2 = {[1 4]; [2 5]};

% factor names
xTickNames = {'1p', '50p'};
xLabel = {'Reward'};

mySubPlots(allMeasures, subInds2, markers, lineColours, measureNames([1, 4]),...
    xTickNames, xLabel, [], [], saveNames{6}, traj, 6, 0, lineStyles , modelTraj, makeSame2);


%% correls

load('./Data/ArrowMotivMotorCorrels.mat','motorRew','wmRew','fields');
plotargs = {'pearson',0,'plot_ci',1,'text',0, 'plotline', 2, 'showzero',1}; % use spearman correlation?

varNames = {'Error (deg)', 'initRT', 'complRT', 'RT (ms)'};
corInds = [1 4]; % only plot error + RT
f = figure();
set(f, 'DefaultAxesFontSize',12);
for i = 1:length(corInds)
    subplot(length(corInds)/2,2,i)
    
    x = col(diff(motorRew.(fields{corInds(i)})));
    y = col(diff(wmRew.(fields{corInds(i)})));
    toRemove = isnan(x) | isnan(y);
    x(toRemove) = [];
    y(toRemove) = [];
    scatterRegress(x, y, plotargs{:});
    %     lsline;
    title(varNames{corInds(i)});
    xlabel({'motor reward effect'});% varNames{corInds(i)}});
    ylabel({'WM reward effect'});% varNames{corInds(i)}});
    
    %     makeSubplotScalesEqual(2,3,[i i+3]);
    
    set(gca,'XTickLabel', cellfun(@num2str, num2cell(xticks),'UniformOutput',0));
    set(gca,'YTick', sort([0 minMax(yticks,2)]));

    axis('square');
end

subplot(length(corInds)/2,2,1);
h = findobj(gca,'Type','Patch');
h.FaceColor = [ .2 .2 .2];
if ~isempty(saveNames{7})
    saveas(f, saveNames{7});
end

%% trajectory modelling

load('./Data/ArrowFitTrajInterp.mat','timeRewIntP','timeResidsData','parNames','timeResidsP')

exptLabels = {'equal','unequal'};
f = figure();
clear h;
stars = p2stars(timeRewIntP);

cols = get(gca,'ColorOrder');
cols = cols(4:5,:);

for iP = 1:4
    subplot(2,2,iP)
    set(gca,'ColorOrder',cols);
    h(:,:,iP) = errorBarPlot(timeResidsData(:,:,:,iP,1),'area',1);
    ylabel([parNames(iP)]);
    xlabel('% of movement time');
    hold on;
    yline(0,':k');
    yl = [-1.2 1.2];
    ylim(yl);
    box off;
    h1 = pbar(sq(timeResidsP(iP,:,1)),'yVal',min(yl),'alpha',.05); % p values per time

    % rew*time int
    text(50 - length(stars{3,iP,1})*3., max(yl)*.8, stars(3,iP,1),  'fontsize', 20);

    xlim([1 100]);
    xticks([1 50 100]);
    xticklabels([0 50 100]);
end

subplot(2,2,1)
legend([h{:,1},h1],{'1p','50p','p<.05'},'Location',[0.3172 0.7731 0.1411 0.1119]);
makeSubplotScalesEqual(2,2);


if ~isempty(saveNames{8})
    saveas(f, saveNames{8});
end
end

function mySubPlots(allMeasures, subInds, markers, lineColours, yLabs, xTickNames, xLabs, leg, legOpts, saveName, traj, iTraj, noNew, lineStyles, modelTraj, makeSame)
% calls different plotting functions to make subplots for one experiment.
% 

if ~exist('noNew','var') || noNew % make new figure
    f = figure();
    set(f,'Position', get(0,'Screensize'));
else
    f = gcf; % reuse current;
end

%% loop through line plots for each measure
if ~isempty(subInds{1})
    n = size(subInds{1},1); % number of measures to plot
    
    for i = 1:n
        subplot(subInds{1}(i,1), subInds{1}(i,2), subInds{1}(i,3)); % current subplot
        if i > 1 % legend on first one
            myFig(allMeasures(:,:,:,i), markers, lineColours, yLabs(i), xTickNames, xLabs, [], [], lineStyles);
        else
            myFig(allMeasures(:,:,:,i), markers, lineColours, yLabs(i), xTickNames, xLabs, leg, legOpts, lineStyles); % do legend once
        end
    end
    
end
%% traj plot - averages for each reward level
if ~isempty(subInds{2})
    subplot(subInds{2}(1,1), subInds{2}(1,2), subInds{2}(1,3));
    myTrajPlot(traj, iTraj, 1, lineColours, lineStyles);
end

%% cross plot - params: guess vs misbinding
if ~isempty(subInds{3})
    subplot(subInds{3}(1,1), subInds{3}(1,2), subInds{3}(1,3));
    myErrorXY(allMeasures(:,:,:,4:5), 0, leg, xLabs, lineColours, yLabs(end-1:end), lineStyles);
end

%% individual data histograms
if ~isempty(subInds{4})
    for i = 1:size(subInds{4},1)
        subplot(subInds{4}(i,1), subInds{4}(i,2), subInds{4}(i,3));
        indivPlot(allMeasures(:,:,:,1), xLabs, leg(3:end), lineColours, lineStyles);
%         modelTrajPlot(modelTraj, iTraj, xLabs(1:2:end), lineColours, i, lineStyles);
    end
    if subInds{4}(i,3) == (subInds{4}(i,1) * subInds{4}(i,2))
%         xlabel('timepoint in movement');
    end
end

%% make axis scales same 

if exist('makeSame','var') && ~isempty(makeSame)
    for i = 1:length(makeSame)
        makeSubplotScalesEqual(subInds{1}(1,1), subInds{1}(1,2), makeSame{i});
    end
end

drawnow;

%% save figure?
if exist('saveName','var') && ~isempty(saveName)
    keyboard;
    saveas(f, saveName);
end

end



function h = myFig(allMeasures, markers, lineColours, yLabs, xTickNames, xLabs, leg, legOpts, lineStyles)
% draw errorbarplot line figs, 2*2*2 design, for a single measure

% plot
h = errorBarPlot(allMeasures, 'type', 'line');
hold on
% set line style
for j = 1:length(h)
    
    h(j).LineWidth = 2;
    %     h(j).Marker = markers{j};
    h(j).Color = lineColours(j,:);
    if size(allMeasures,2) == 4 % break the line into pairs
        h(j).LineStyle = 'none';
        plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:), 'LineStyle',lineStyles{1});
        plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:), 'LineStyle',lineStyles{2});
    end
    
end

% rew markers
x = [h(1,:).XData];
y = [h(1,:).YData];
for i = 1:2
    
    h1(i) = plot(x(i:2:end), y(i:2:end),'Marker',markers{mod(i-1,2)+1},...
        'MarkerFaceColor',lineColours(mod(i-1,2)+3,:),'LineStyle','none',...
        'LineWidth',2,'Color',lineColours(mod(i-1,2)+3,:));
end

set(gca,'XTick',1:length(xTickNames),'XTickLabel',xTickNames)
xlabel([xLabs{:}])
xlim([0.5 length(xTickNames)+.5])
ylabel(yLabs)
%         yticks( prctile(yticks, [0 50 100]) )
box off

if exist('leg','var') && ~isempty(leg)
      
    opts = [markers,lineStyles(1), lineStyles(1), lineStyles(1:2)];
    cols = {lineColours(3,:),lineColours(4,:),lineColours(1,:),lineColours(2,:), [0 0 0], [0 0 0]};
    if size(allMeasures,2)==4, n=6; else n=4;end
    emptyLegend(n, opts, cols, leg, legOpts);
    
end


end


function myTrajPlot(traj, iTraj, leg, lineColours, lineStyles)
% trajectories conditional plot

% reward colours
cols = [0 0 0; .5 .5 .5];
set(gca,'ColorOrder', flipud(lineColours(3:4,:)));

[x,y,~,~,h] = conditionalPlot(traj.o2{iTraj}.colTrajTimes, traj.o2{iTraj}.colTargDist);
hold on;

if traj.sepConds
    y1 = y;
    yline(0,'k:');
    if size(h,1) > 1
        h(2,1).LineStyle = lineStyles{2};
        h(2,2).LineStyle = lineStyles{2};
        h(2,2).EdgeColor = cols(2,:);
    end
else
    y1 = diff(y,[],2);
end
%         ybar = max(abs(ylim));

useClust=0; nPerms=5000;
hold on;
yMin = min([col([h(:,2).YData]); min(ylim)]); % get highest point
yBit = 0;%abs(diff(ylim))/25; % 20% of ylims
for i = 1:size(y1,2)
    [~,pp(i,:)]=permutationOLS( sq(y1(:,i,:)), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
    pbar(pp(i,:), 'yVal', yMin + yBit*i, 'xVal', sq(x(1,1,:)), 'plotargs', {'Color', cols(i,:), 'LineWidth',5});
end
    

% x axis
h = findobj(gca,'Type','Line');
x = sq(x(1,1,[1 end]))';
set(gca, 'XTick', sort([mean(x) x]), 'XTickLabel', [0 50 100])
xlabel('% time through movement')
xlim(x)

% xlabel('timepoint in movement');
ylabel(traj.ylab)
box off

% legend
if exist('leg','var') && ~isempty(leg)
    h = [findobj(gca,'Type','Patch'); h];
    legend(h,[traj.condLabels{1}(1,:), 'p < .05'],'Location','Best');
end

end



function myErrorXY(X, semBars, leg, xLabs, lineColours, parNames, lineStyles, markArgs)
% xy error bar cross plot

%%

if ~exist('lineStyles','var') || isempty(lineStyles)
    lineStyles = {'k-','k--'};
    %     cols = get(gca,'ColorOrder');
end
if ~exist('markArgs','var') || isempty(lineStyles)
    markArgs = {'o','^'};
end

%% transform data

% it is currently [pp, rew(1,2,1,2), factors]
[pp,c,p,n] = size(X);
X = reshape(X, pp, 2, c/2, p, n);
% [pp, rew, cue, delay/etc, param]

%%
means = nanmean(X);
[~,r,c,p,~] = size(means);
sems = nanstd(X) ./sqrt(pp);

% for i = 1:c % over cond1
%     for j = 1:p  % cond2
%         plot(means(1,:,i,j), means(1,:,i,j,2), lineArgs{1,i}, 'Color', lineArgs{2,j}, 'LineWidth', 2);
%         hold on;
%     end
% end


for i = 1:c % draw crosses
    for j = 1:p
        plot(means(1,:,i,j), means(1,:,i,j,2), 'Color', lineColours(j,:), 'LineWidth', 2, 'LineStyle', lineStyles{1,i});
        
        for k = 1:r
            hold on;
            if semBars
                errorbarxy(means(1,k,i,j), means(1,k,i,j,2), sems(1,k,i,j), sems(1,k,i,j,2), sems(1,k,i,j), sems(1,k,i,j,2), {'.'}, {'color', lineColours(k+2,:),'LineWidth',2});
            end
            hold on;
            
            % markers
            plot(means(1,k,i,j), means(1,k,i,j,2), markArgs{k}, 'Color', lineColours(k+2,:), 'MarkerFaceColor', lineColours(k+2,:),'MarkerSize',8);
            
        end
    end
end
xlabel(parNames(1))
ylabel(parNames(2))
box off

% equal axes
ax = [xlim; ylim];
ax = [min(ax(:,1)) max(ax(:,2))];
axis([ax ax]);

% equal ticks
xt = xticks;
yt = yticks;
if length(xt) < length(yt)
    set(gca,'YTick',xt);
else
    set(gca,'XTick',yt);
end
%% legend
if exist('leg', 'var') && ~isempty(leg)
    if (r+c+p) == 6, n=6; else n = 4; end % num factors/points
    opts = [markArgs,lineStyles(1), lineStyles(1), lineStyles(1:2)];
    cols = {lineColours(3,:),lineColours(4,:),lineColours(1,:),lineColours(2,:), [0 0 0], [0 0 0]};
    emptyLegend(n, opts, cols, leg);
end
end



function indivPlot(measure,xlab, leg, lineColours, lineStyles)
% ksdensity for histograms plot

% get low and hig rew conds
low = measure(:,1:2:end,:);
high = measure(:,2:2:end,:);

rewEff = (high - low); % reward effect

[~,c,p] = size(rewEff);

cols = lineColours;

for i = 1:p
    for j = 1:c
        hold on;
        [y,x] = ksdensity(rewEff(:,j,i));
        plot(x,y, lineStyles{j}, 'Color', cols(i,:), 'LineWidth',2);
        xline(nanmedian(rewEff(:,j,i)),lineStyles{j}, 'Color', cols(i,:), 'LineWidth',1.5);
    end
end
% xline(0,'-.k','LineWidth', 2);

xlabel('reward effect on error');
ylabel('frequency');
box off

xl = max(abs(minMax(xlim,2)));
xlim([-xl xl]); 
%% legend stuff
if exist('leg', 'var') && ~isempty(leg)
    n = p*c;
    opts = [lineStyles(1), lineStyles(1), lineStyles(1,:)];
    cols = {lineColours(1,:),lineColours(2,:), [0 0 0], [0 0 0]};
    emptyLegend(n, opts, cols, leg);

end

end


function modelTrajPlot(traj, iTraj, leg, lineColours, iPar, lineStyles)
% trajectories conditional plot - modelling stuff

% reward colours
cols=[0 0 0; 0.5 0.5 0.5];
set(gca,'ColorOrder', cols);


parsConds = traj.allParsConds{iTraj}(:,:,:,iPar);

h = errorBarPlot(parsConds,'area',1);


% perm test
if traj.sepConds
    y1 = sq(parsConds);
    yline(0,'k:');
    if size(h,1) > 1
        h{2,1}.LineStyle = lineStyles{2};
        h{2,2}.LineStyle = lineStyles{2};
        h{2,2}.EdgeColor = cols(2,:);
    end
else
    y1 = sq(diff(parsConds,[],3));
end
%         ybar = max(abs(ylim));

useClust=1; nPerms=5000;
hold on;
yMax = max([col(cellfun(@(x) max(x.YData), h)); max(ylim)]); % get highest point
yBit = abs(diff(ylim))/25; % 20% of ylims
for i = 1:size(y1,3)
    [~,pp(i,:)]=permutationOLS( y1(:,:,i), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
    pbar(pp(i,:), 'yVal', yMax + yBit*i, 'plotargs', {'Color', cols(i,:), 'LineWidth',3});
end

%         for k = 1:size(parsConds,3)
%             [~,pp]=permutationOLS( parsConds(:,:,k,i), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);

% draw pbar
%             hold on;
%             pbar(pp, 'yVal', ybar(i).*(1+k/20), 'plotargs', {'Color', h{k,1}.Color, 'LineWidth',3});
%         end
% yline(0,':k');
% xlim([0.5 traj.nInterps + .5])

x = [1 traj.nInterps];
set(gca, 'XTick', x, 'XTickLabel', {'Start', 'Finish'})
xlim(x + [-1 1])

box off

ylabel(['reward effect: ' traj.parNames{iPar}]);
% if iPar==2; xlabel('timepoint in movement'); end
if iPar==1 && size(h,1)>1; legend([h{:,1}], leg, 'Location','Best');end

end