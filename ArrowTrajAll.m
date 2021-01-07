% ArrowTrajAll

close all; clear; clc;
%%

matFiles = {'ArrowMotivAnalysis.mat';
            'ArrowMotivCueBlockedAnalysis.mat';
            'ArrowMotivRetrocueAnalysis.mat';
            'ArrowMotivReallocateAnalysis.mat';
            'ArrowMotivMotorAnalysis_A.mat';
            'ArrowMotivMotorAnalysis_B.mat';
            };
        
dataFolder = './Data/';

n = 6; % length(matFiles);

condCols = {[1 1 1 1 2 2 2 2; 1 1 2 2 1 1 2 2; 1 2 1 2 1 2 1 2]; % rew, cue, delay
            [1 1 2 2 1 1 2 2; 1 1 1 1 2 2 2 2; 1 2 1 2 1 2 1 2]; % rew, cue, catch
            [1 1 2 2;         1 2 1 2];                          % rew, congr
            [1 1 2 2 1 1 2 2; 1 1 1 1 2 2 2 2; 1 2 2 1 1 2 2 1]; % rew, cue, equal
            [1 2];                                               % rew
            [1 2];                                               % rew
            };

condCols2 = {   [1 1 1 1]; % for averaging across some conds
                [1 1 1 1];
                [1 1];
                [1 3; 2 4];
                [];
                [];
                };
            
condLabels = { {'1p','50p';     'pre','post';       'short','long'          };
               {'1p','50p';     'pre','post';       'notCatch', 'isCatch'   };
               {'1p','50p';     'incongr','congr';                          };
               {'1p','50p';     'pre','post';       'equal', 'unequal'      };
               {'1p','50p'};
               {'1p','50p'};
               };
        
legends = { {'preShort','preLong','postShort','postLong';}
            {'pre','preCatch','post','postCatch';}
            {'incon','con';}
            {'pre-eq','pre-uneq','post-eq','post-uneq';}
            {''};
            {''};
            };

exptNums = {'1','2','3','4','5a','5b'};
doPlots = 0;
for i = 1:n
    
    d = load(fullfile(dataFolder, matFiles{i}), 'dataFiles', 'dataStruct');
    
    o{i} = ArrowTrajAnalysis(d.dataFiles, d.dataStruct, doPlots);
    
    o{i}.condCols = condCols{i,:}(1,:);
    o{i}.condCols2 = condCols2{i,:};

    o{i}.condLabels = condLabels{i,:}(1,:);
    o{i}.legends = legends{i,:};
    
    
end

% remove initial two trials for expt 5a
o{5}.interpTargDists(:,1:2,:) = NaN;


%%
sepConds=0;

xlab = 'time point since initiation';
ylab = 'angular error (deg)';
figure();
% n = length(o);
for i = 1:n
    subplot(ceil(n/2),2,i)
    

    % get o.condCols, o.labels
%     o1{i} = PlotTrajPrecVsTime(o{i});
    o{i} = PlotTrajPrecVsTimeInterp(o{i}, sepConds);
    
    title(['Experiment ', exptNums{i}]);
    drawnow;
    
    if i > n-2
        xlabel(xlab);
    end
    if mod(i,2)
        ylabel(ylab)
    end
    
    o2{i}.colTrajTimes = o{i}.colTrajTimes;
    o2{i}.colTargDist = o{i}.colTargDist;
end

% h = findobj(gca,'Type','Patch');
% legend(flipud(h(1:2)),condLabels{1}(1,:),'Location','Best');

%%

save('./Data/ArrowTrajAll.mat','o2', 'condCols','condLabels','exptNums', 'xlab', 'ylab','sepConds','-v7.3');