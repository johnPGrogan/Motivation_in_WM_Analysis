function o = PlotTrajPrecVsTimeInterp(o, sepConds)
% Subplots conditionalPlots split by several conditions, plotting
% interpolated trajectory time vs absolute error
% o = struct containing o.allTrajTimes, o.splitBy, o.allTargDists,
% o.condCols, o.condLabels, 
% This can be obtained from ArrowTrajAnalysis()
% 
% returns o = workspace2struct();

warning('off', 'NANCAT:emptyCells')
%% plot targDist vs time

o.splitBy = nancat(1, o.dataStruct.trialTypes);

% uncomment this to plot interp points (i.e. same time scale)
o.interpTimePoints = permute(repmat([1:100],o.nPP, 1, o.nTrials),[1 3 2]);

% o.interpTimePoints = o.interpTimePoints * 1000;

% split by trial type
o.allTrajTimesRew = groupMeans(o.interpTimePoints, 2, repmat(o.splitBy,1,1,size(o.interpTargDists,3)),'dim');
o.allTargDistRew = abs(rad2deg(groupMeans(o.interpTargDists, 2, repmat(o.splitBy,1,1,size(o.interpTimePoints,3)),'dim')));

% uncomment this to normalise each start error to 100
% o.allTargDistRew = abs(o.allTargDistRew) ./ nanmean(abs(o.allTargDistRew(:,:,:,1)),[1 2]) * 100;


n = size(o.condCols,1);

% for each cond split
for i = 1:n
    
    
    % split by cond, and reshape and permute into [trajT x trial x type, pp, cond]
%     for j = 1:2
        o.allTargDistRewEff = cat(5, o.allTargDistRew(:,o.condCols(i,:)==2,:,:),  o.allTargDistRew(:,o.condCols(i,:)==1,:,:));
        o.allTrajTimesRewEff = cat(5, o.allTrajTimesRew(:,o.condCols(i,:)==2,:,:),  o.allTrajTimesRew(:,o.condCols(i,:)==1,:,:));
%     end
    
    if exist('sepConds','var') && sepConds
        nConds = size(o.allTargDistRewEff,2);

        % plot each cond
        if any(o.condCols2)
            x = [];
            for j = 1:size(o.condCols2,1)
                x = cat(6, x, o.allTargDistRewEff(:,o.condCols2(j,:),:,:,:));
            end
            o.allTargDistRewEff = x;
            
            nConds = size(o.condCols2,1);

        end
            
            
        
        o.colTargDist = reshape(permute(diff(o.allTargDistRewEff,[],5),[4,3,1,2,5,6]),[],o.nPP,nConds);
        o.colTrajTimes = reshape(permute(o.allTrajTimesRewEff(:,:,:,:,1),[4,3,1,2]),[],o.nPP,nConds);
        [x,y] = conditionalPlot(o.colTrajTimes, o.colTargDist);
        h = findobj(gca, 'Type', 'Patch'); % the patches are upside down so flip legend order

        useClust=0; nPerms=5000;
        x1 = sq(nanmean(x,[1,2])); % get xvalues
        hold on;
        yMax = max([col([h(:).YData]); max(ylim)]); % get highest point
        yBit = abs(diff(ylim))/25; % 20% of ylims
        for j=1:nConds
            [~,pp(j,:)]=permutationOLS( sq(y(:,j,:)), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
            pbar(pp(j,:), 'yVal',yMax + yBit*j,'xVals',x1,'plotargs',{'Color', h(nConds-j+1).FaceColor,'LineWidth',3});
        end
        yline(0,':k');
        if nConds>1
            legend(flipud(h), o.legends, 'Location', 'Best'); 
        end

        
    else
        o.colTargDist = reshape(permute(o.allTargDistRewEff,[4,3,2,1,5]),[],o.nPP,2);
        o.colTrajTimes =  reshape(permute(o.allTrajTimesRewEff,[4,3,2,1,5]),[],o.nPP,2);


        [x,y] = conditionalPlot(o.colTrajTimes, o.colTargDist);
%         [x,y] = conditionalPlot(o.colTrajTimes(:,:,1), diff(o.colTargDist,[],3));
        h = findobj(gca, 'Type', 'Patch'); % the patches are upside down so flip legend order
    %     legend(flipud(h(1:2)), [o.condLabels(i,:) 'p < .05'],'Location','Best'); 

        % use cluster permutation testing
        y1 = sq(diff(y,[],2));
        useClust=0; nPerms=5000;
        [~,pp]=permutationOLS( y1, [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);

        hold on;
        x1 = sq(nanmean(x(:,1,:),[1,2])); % get xvalues
        pbar(pp, 'yVal', min(ylim),'xVals',x1);
    end
%     xlabel('time from resp start (s)')
%     ylabel('abs(error)')
    box off
    xticks(prctile(xticks,[0 50 100]));
    yticks(prctile(yticks,[0 50 100]));

        
end

warning('on')
end