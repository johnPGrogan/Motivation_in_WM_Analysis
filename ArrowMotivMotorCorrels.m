% ArrowMotivMotorCorrels

% load up all motor control data
% load up data from those people on WM tasks
% correlate motor and WM errors and RTs, to see how much of WM is due to
% motor effects

%%
clear
close all

dataFolder = '../../Experiments/Results/ArrowMotivCue';%folder with data in
files1 = what(dataFolder);%find files there

cellregexpi = @(cellArray, pattern) ~cellfun(@isempty, regexpi(cellArray, pattern));

motorFiles = files1.mat(cellregexpi(files1.mat,'ArrowMotivMotorControl_'));%keep only files with this in
nPP = length(motorFiles);
%% get ppIDs
cellRegexpiInds = @(cellArray, pattern) cell2mat(regexpi(cellArray, pattern));
x = cellRegexpiInds(motorFiles, 'Control_');
y = cellRegexpiInds(motorFiles, '.mat');
for i = 1:nPP
    motorIDs{i,1} = motorFiles{i}(x(i)+8:y(i)-1);
end

%% find other files with those IDs

for i = 1:nPP
    wmFiles(i,1) = files1.mat(cellregexpi(files1.mat, motorIDs{i}) & ~cellregexpi(files1.mat, 'ArrowMotivMotorControl_'));
    wmExpt{i,1} = wmFiles{i}(1: cellRegexpiInds(wmFiles(i), '_')-1); % get expt type
end

% get number for groupMeans
wmExptNames = {'ArrowMotivCueBlocked', 'ArrowMotivReallocate', 'ArrowMotivRetrocue'};
nExpts = length(wmExptNames);
for j = 1:nExpts
    wmExptNum(:,j) = strcmp(wmExpt, wmExptNames{j});
end
wmExptNum = wmExptNum * [1:length(wmExptNames)]';

%% split wm and motor by that

motorFiles1 = groupMeans(motorFiles, 1, wmExptNum, 'dim')';
wmFiles1 = groupMeans(wmFiles, 1, wmExptNum, 'dim')';

nPerExpt = size(wmFiles1, 1);

%% load each expt type and get rew effect on error and RT

missing = cellfun(@isempty, motorFiles1);
for i = 1:nExpts
    motorStruct{i} = ArrowDataLoad(fullfile(dataFolder, motorFiles1(~missing(:,i),i)));%analyse each file
    wmStruct{i} = ArrowDataLoad(fullfile(dataFolder, wmFiles1(~missing(:,i),i)));
end

%% get prec, RTs for each col

[motor, wm] = deal(struct());

fields = {'prec','startRT','endRT','totalRT'};
nFields = length(fields);
for j = 1:nFields-1
    motor.(fields{j}) = [];
    wm.(fields{j}) = [];
    for i = 1:nExpts
        motor.(fields{j}) = nancat(4,  motor.(fields{j}), nancat(3, motorStruct{i}.(fields{j})));
        wm.(fields{j}) = nancat(4,  wm.(fields{j}), nancat(3, wmStruct{i}.(fields{j})));
        % [trials, conds, pp, expt]
    end
end

motor.totalRT = motor.startRT + motor.endRT;
wm.totalRT = wm.startRT + wm.endRT;

% convert to degrees
motor.prec = rad2deg(motor.prec);
wm.prec = rad2deg(wm.prec);

% % make sure they are mod()
% motor.prec = mod(motor.prec+180,360)-180;
% wm.prec = mod(wm.prec+180,360)-180;

%% average over diff conds
motorConds = {'low','high'; 'low','high'};
wmConds = nancat(1, {'lowPre','lowPreCatch','highPre','highPreCatch','lowPost','lowPostCatch','highPost','highPostCatch'},...
    {'lowMotivPre','lowRealloPre','highRealloPre','highMotivPre','lowMotivPost','lowRealloPost','highRealloPost','highMotivPost'},...
    {'lowIncon','lowCon','highIncon','highCon';}...
    )';

motorRewInds = repmat(permute([0 1; 0 1; 0 1]', [4,1,3,2]), [size(motor.prec,1) 1 size(motor.prec,3),1]); % 1 = high rew
wmRewInds = repmat(permute([0 0 1 1 0 0 1 1; 0 NaN NaN 1 0 NaN NaN 1; 0 0 1 1 NaN NaN NaN NaN]', [4,1,3,2]), [size(wm.prec,1),1,size(wm.prec,3),1]);

motorRew = structfun(@(x) permute(nanmean(x,1), [3,2,4,1]), motor,'UniformOutput',0);
wmRew = structfun(@(x) permute(nanmean(groupMeans(x, 2, wmRewInds,'dim'),[1 5]), [3,2,4,1]), wm,'UniformOutput',0);
% [rew, pp, expt, other conds]


%% scatter regress

plotargs = {'pearson',0,'plot_ci',2,'text',2, 'plotline', 2, 'showzero',1}; % use spearman correlation?

figure();
c = get(gca,'ColorOrder');
for i = 1:nFields
    subplot(2,nFields,i)
    
    for j = 1:3
        [~,~,~,~,~,h(j,:)] = scatterRegress(diff(motorRew.(fields{i})(:,:,j),[],2), [ diff(wmRew.(fields{i})(:,:,j),[],2)], plotargs{:});
        h(j,1).CData = c(j,:);
        hold on;
    end
    title(fields{i});
    xlabel('motor rew effect');
    ylabel('WM rew effect');
end


% combine across expts


for i = 1:nFields
    subplot(2,nFields,i+nFields)
    
    [~,~,~,~,stats(i,:),h(i,:)] = scatterRegress(col(diff(motorRew.(fields{i}))), col(diff(wmRew.(fields{i}))), plotargs{:});
%     lsline;
    title(fields{i});
    xlabel('motor rew effect');
    ylabel('WM rew effect');
%     makeSubplotScalesEqual(2,3,[i i+3]);
end


% legend(h(:,1), wmExptNames, 'Location',[0.5731 0.3183 0.2696 0.1048])

%% exclude an outlier
% one person (TRI_015) was super slow in endRT on motor

replaceNaN = ones(size(motorRew.prec));
replaceNaN(29,:,3) = NaN;

motorRew = structfun(@(x) x.*replaceNaN, motorRew,'UniformOutput',0);

% now redo figs
figure();
for i = 1:nFields
    subplot(2,nFields,i)
    
    for j = 1:3
        % stats is [R^2, F, p, error var est]
        [~,~,~,~,~,h(j,:)] = scatterRegress(diff(motorRew.(fields{i})(:,:,j),[],2), [ diff(wmRew.(fields{i})(:,:,j),[],2)], plotargs{:});
        h(j,1).CData = c(j,:);
        hold on;
    end
    title(fields{i});
    xlabel('motor rew effect');
    ylabel('WM rew effect');
end

% combine across expts

for i = 1:nFields
    subplot(2,nFields,i+nFields)
    
    [~,~,~,~, stats1(i,:)] = scatterRegress(col(diff(motorRew.(fields{i}))), col(diff(wmRew.(fields{i}))), plotargs{:});
%     lsline;
    title(fields{i});
    xlabel('motor rew effect');
    ylabel('WM rew effect');
%     makeSubplotScalesEqual(2,3,[i i+3]);

end

% legend(h(:,1), wmExptNames, 'Location',[0.5731 0.3183 0.2696 0.1048])

%%
plotargs = {'pearson',0,'plot_ci',1,'text',0, 'plotline', 2, 'showzero',1}; % use spearman correlation?

varNames = {'error', 'initRT', 'complRT', 'totalRT'};
f = figure();
set(f, 'DefaultAxesFontSize',12);
for i = 1:nFields
    subplot(2,2,i)
    
    x = col(diff(motorRew.(fields{i})));
    y = col(diff(wmRew.(fields{i}))); 
    toRemove = isnan(x) | isnan(y);
    x(toRemove) = [];
    y(toRemove) = [];
    [~,~,~,~, stats2(i,:)] = scatterRegress(x, y, plotargs{:});
%     lsline;
%     title(fields{i});
    xlabel(['motor rew effect: ' varNames{i}]);
    ylabel(['WM rew effect: ' varNames{i}]);
%     makeSubplotScalesEqual(2,3,[i i+3]);

end

subplot(2,2,1);
h = findobj(gca,'Type','Patch');
h.FaceColor = [ .2 .2 .2];


%%

save('./Data/ArrowMotivMotorCorrels.mat');