function output = ArrowRaceAnalysis(output)

% close all

%%
struct2workspace(output);


pars2 = cell2mat(permute(pars,[1,4,2,3]));

RT2 = cell2mat(permute(RT,[4,1,2,3]));
Corr2 = cell2mat(permute(Corr,[4,1,2,3]));

nRTs=length(rtCols);

nFactors = size(factorNames,2)-1;

rtLabels = rtLabels(rtCols);

%% draw figure

% rtCols = 1;
% rtLabel = {'startRT','endRT','totalRT'};

% condCols = {[1 2; 3 4]; [1 3; 2 4]}; %%%%%%% use this to split conds
if drawKs
    c = get(gca,'ColorOrder');
    vals = [1 50];

    rt = permute(groupMeans(RT2,1,Corr2,'dim'), [5,2,3,1,4]);
    for iRT = 1:nRTs
        figure();clf
        for j = 1:nFactors
            for i = 1:2
                subplot(2,nFactors,j)

                [f,x] = ksdensity(col(rt(:,:,condCols{j}(i,:),2,iRT)),[0:6000]);
                fill([x, fliplr(x)], [f, zeros(size(f))],c(i,:),'FaceAlpha',0.5);
                hold on

                [f,x] = ksdensity(col(rt(:,:,condCols{j}(i,:),1,iRT)),[0:6000]);
                fill([x, fliplr(x)], -[f, zeros(size(f))],c(i,:),'FaceAlpha',0.5);
                hold on
                if j==1
                    ylabel('model')
                end

                % now plot real data
                a = [data{:,condCols{j}(i,:),iRT}];
                r1  = cat(1,a.rt_correct);
                rt1 = nancat(1, r1{:});
                r2  = cat(1,a.rt_error);
                rt2 = nancat(1, r2{:});

                title(factorNames{j+1})

                subplot(2,nFactors,j+nFactors)
                [f,x] = ksdensity(rt1,[0:6000]);
                fill([x, fliplr(x)], [f, zeros(size(f))],c(i,:),'FaceAlpha',0.5);
                hold on
                
                if ~isempty(rt2) % only plot if there are incorrect resps
                    [f,x] = ksdensity(rt2,[0:6000]);
                    fill([x, fliplr(x)], -[f, zeros(size(f))],c(i,:),'FaceAlpha',0.5);
                    hold on
                end
                if j==1
                    ylabel('data')
                end
                xlabel('RT (ms)')
            end
        end
        SuperTitle(rtLabels{iRT});
%         saveas(figure(iRT),sprintf('%s_%d.jpg',saveName,rtCols(iRT)))
    end

end
%% stats on pars

% rew = permute(repmat([0 0 1 1],nPP,1,1),[1,3,2]);
% cue = permute(repmat([0 1 0 1],nPP,1,1),[1,3,2]);
% pp  = permute(repmat(1:nPP, 1,1,nConds), [2,1,3]);
% 
% modelMat = [0 0 1; 1 0 0; 0 1 0;
%     1 1 0;];
% factorNames = {'rew','cue','pp'};
% factors = {col(rew),col(cue),col(pp)};

for iRT = 1:nRTs
    for i = 1:4
%         [stats(i,iRT).p,stats(i,iRT).tab,...
%             stats(i,iRT).s] = anovan(col(pars2(:,i,:,iRT)),factors,...
%             'varnames',factorNames,'display','off',...
%             'model','full','random',length(factorNames),'model',modelMat);%anova, with factors and labels, no figure, a random effect of ppInd, and the modelMat terms
        
        % rmanova
        x = reshape( sq(pars2(:,i,:,iRT)), reshapeDims);
        rmModelStats{i,iRT} = rmanova(x, factorNames);
        pVals(:,i,iRT) = rmModelStats{i,iRT}.pValue;
    end
end

% pVals = [];
% for iRT=1:nRTs
%     pVals = nancat(3, pVals, [stats(:,iRT).p]);
% end
alphas = [.05, .01, .001, .0001];

nEffects = size(pVals,1) - 1;

nStars = zeros(nEffects,4,nRTs);
for iRT = 1:nRTs
    for i = 1:4
        nStars(:,i,iRT) = sum(pVals(2:end,i,iRT) < alphas,2);
        for j = 1:nEffects
            stars{j,i,iRT} = repmat('*',1,nStars(j,i,iRT));
        end
    end
end

allStats = table();
allStats1 = table();
parNames = {'tA','dA','ndt','thresh'};

for iRT = 1:nRTs
    for i = 1:4
        allStats1.measure = repmat(rtLabels(iRT),nEffects,1);
        allStats1.param = repmat(parNames(i),nEffects,1);
        allStats1.effect = rmModelStats{i,iRT}.Term(2:end);
        allStats1.df = rmModelStats{i,iRT}.DF1(2:end);
        allStats1.dfErr = rmModelStats{i,iRT}.DF2(2:end);
        allStats1.F = rmModelStats{i,iRT}.FStat(2:end);
        allStats1.p = rmModelStats{i,iRT}.pValue(2:end);
        allStats1.stars = stars(:,i,iRT);

        allStats = vertcat(allStats, allStats1);
    end  
end
disp(allStats(:,[1 2 3 7 8]))

%% motiv/realloc split

if exist('splitAnova','var')
    pVals2 = [];
    for j = 1:2
        for iRT = 1:nRTs
            for i = 1:4
                % rmanova
                x = reshape( sq(pars2(:,i,condCols{3}(j,:),iRT)), reshapeDims(1:end-1));
                rmModelStats2{j,i,iRT} = rmanova(x, factorNames);
                pVals2(:,j,i,iRT) = rmModelStats2{j,i,iRT}.pValue;
            end
        end
    end

    alphas = [.05, .01, .001, .0001];

    nStars2 = zeros(3,2,4,4);
    stars2 = cell(3,2,4,4);
    for j = 1:2
        for iRT = 1:nRTs
            for i = 1:4
                nStars2(:,j,i,iRT) = sum(pVals2(2:end,j,i,iRT) < alphas,2);
                for k = 1:3
                    stars2{k,j,i,iRT} = repmat('*',1,nStars2(k,j,i,iRT));
                end
            end
        end
    end

    allStats2 = table();
    allStats1 = table();
    nEffects = size(pVals2,1)-1;
    
    condLabels = {'motivate','reallocate'; 'low','high';'pre','post'};
    for j = 1:2
        for iRT = 1:nRTs
            for i = 1:4
                allStats1.cond = repmat(condLabels(1,j),3,1);
                allStats1.measure = repmat(rtLabels(iRT),nEffects,1);
                allStats1.param = repmat(parNames(i),nEffects,1);
                allStats1.effect = rmModelStats2{j,i,iRT}.Term(2:end);
                allStats1.df = rmModelStats2{j,i,iRT}.DF1(2:end);
                allStats1.dfErr = rmModelStats2{j,i,iRT}.DF2(2:end);
                allStats1.F = rmModelStats2{j,i,iRT}.FStat(2:end);
                allStats1.p = rmModelStats2{j,i,iRT}.pValue(2:end);
                allStats1.stars = stars2(:,j,i,iRT);

                allStats2 = vertcat(allStats2, allStats1);
            end
        end
    end

    disp(allStats2(:,[1,2,3,4,8,9]))%,11,12]))
    output.allStats2 = allStats2;

end
%% mean pars

if drawPars
    figure()
    clf;
    markers = {'o','^'};
    lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
    % condInds = [1 3 2 4];

    for iRT = 1:nRTs
        for i = 1:4
            subplot(4,nRTs,(i-1)*nRTs+iRT)


    %         for k = 1:size(condInds,1)
    %             h = errorBarPlot(nancat(3, sq(pars2(:,i,condInds{k,1},iRT)), ...
    %                 sq(pars2(:,i,condInds{k,2},iRT))),'type','line','xaxisvalues',[k*2-1;k*2]);
                h = errorBarPlot( nancat(3, sq(pars2(:,i,condInds{1},iRT)), sq(pars2(:,i,condInds{2},iRT))),'type','line');
                hold on
                for j = 1:length(h)
                    h(j).LineWidth = 2;
                    h(j).Marker = markers{j};
                    h(j).Color = lineColours(j,:);
                    if length(condInds{1}) == 4
                        h(j).LineStyle = 'none';
                        plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:));
                        plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:));
                    end
                end
    %         end
            set(gca,'XTick',1:length(xTickLabels),'XTickLabel',xTickLabels)
            xlabel(xLabel)
            xlim([0.5 length(xTickLabels) + .5])
            ylabel(parNames{i})

            box off
            if i==1
                title(rtLabels{iRT})
            end

            if length(factorOrder)>0
                text(0.5,0.8,stars{factorOrder(1),i,iRT},'Units','Normalized','FontSize',12,'Color',[.85 .325 .098]);
                if length(factorOrder)>1
                    text(0.25,0.7,stars{factorOrder(2),i,iRT},'Units','Normalized','FontSize',12);
                    if length(factorOrder)>2
                        text(0.5,0.6,stars{factorOrder(3),i,iRT},'Units','Normalized','FontSize',12);
                    end
                end
            end
        end
    end
    if exist('legendLabels', 'var')
        legend(legendLabels,'Location',[.007 .264 .125 .1])
    end
    for i = 1:4
        makeSubplotScalesEqual(4,nRTs,i*(nRTs)-(nRTs-1):i*(nRTs))
    end
    % makeSubplotScalesEqual(4,nRTs,nRTs+1:nRTs*2)

end

%% plot median RTs for misbinds

if drawRTs
    figure()
    clf;
    markers = {'o','^'};
    lineColours = [0 0.447 0.741; 0.85 0.325 0.098];
    % condInds = [1 3 2 4];
    paramNames = {'target RT','guess RT', 'misbind RT'};
    
    for iRT = 1:nRTs
        for i = 1:3
            subplot(3,nRTs,(i-1)*nRTs+iRT)
            
            d = nancat(3, rtSplit2(:,condInds{1},i,:,iRT), rtSplit2(:,condInds{2},i,:,iRT));
            d = permute(nanmedian(d,1),[4,2,3,1]);

                h = errorBarPlot( d,'type','line');
                hold on
                for j = 1:length(h)
                    h(j).LineWidth = 2;
                    h(j).Marker = markers{j};
                    h(j).Color = lineColours(j,:);
                    if length(condInds{1}) == 4
                        h(j).LineStyle = 'none';
                        plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:));
                        plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:));
                    end
                end
    %         end
            set(gca,'XTick',1:length(xTickLabels),'XTickLabel',xTickLabels)
            xlabel(xLabel)
            xlim([0.5 length(xTickLabels) + .5])
            ylabel(paramNames{i})

            box off
            if i==1
                title(rtLabels{iRT})
            end
% 
%             if length(factorOrder)>0
%                 text(0.5,0.8,stars{factorOrder(1),i,iRT},'Units','Normalized','FontSize',12,'Color',[.85 .325 .098]);
%                 if length(factorOrder)>1
%                     text(0.25,0.7,stars{factorOrder(2),i,iRT},'Units','Normalized','FontSize',12);
%                     if length(factorOrder)>2
%                         text(0.5,0.6,stars{factorOrder(3),i,iRT},'Units','Normalized','FontSize',12);
%                     end
%                 end
%             end
        end
    end
    if exist('legendLabels', 'var')
        legend(legendLabels,'Location',[.007 .264 .125 .1])
    end
    for iRT = 1:nRTs
        makeSubplotScalesEqual(3,nRTs,iRT:nRTs:(3*nRTs))
    end
    % makeSubplotScalesEqual(4,nRTs,nRTs+1:nRTs*2)

end

%%

output.allStats = allStats;
% save('../Data/ArrowMotivCueBlockedRace.mat','output')

end