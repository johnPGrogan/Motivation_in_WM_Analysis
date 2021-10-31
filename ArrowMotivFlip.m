% ArrowMotivFlip
% Look at the 180 flipped model fittings

%% 
dataFolder = './Data';

loadFiles ={'ArrowMotivModelAnalysis180.mat';
            'ArrowMotivCueBlockedModelling180.mat';
            'ArrowMotivRetrocueModelling180.mat';
            'ArrowMotivReallocateModelling180.mat';
            };
        
nExpts = 4;

clear d;
for i = 1:nExpts
    d1 = load(fullfile(dataFolder, loadFiles{i}),'meanBIC', 'allPars', 'n', 'allStats', 'allStats2');
    if i> 1
        [d, d1] = ensureStructsAssignable(d, d1);
    end
    d(i) = d1;
end

% d = nancat(d);

%% extract

bic = nancat(1, d.meanBIC);
n = nancat(1, d.n);
pars = nancat(4, d.allPars);


%% hists
parNames = {'imprecision','target','guessing','misbinding','flips'};
nPars = length(parNames);

figure();
for iP = 1:nPars
    
    subplot(3,2,iP);
    hist(col(pars(:,:,iP,:)), 100);
    title(parNames{iP});
end
makeSubplotScalesEqual(3,2,2:5);