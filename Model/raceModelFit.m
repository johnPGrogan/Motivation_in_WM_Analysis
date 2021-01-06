function [pars,nll,Corr,RT] = raceModelFit(data,useParfor,fitWithin)
% use gridsearch and then fminsearch to fit LBA race model to data split by
% condition

if ~exist('useParfor','var')
    useParfor=0;
end
if ~exist('fitWithin','var')
    fitWithin=0;
end

pars1 = gridSearch(data,useParfor,fitWithin);

%%
[nll,pars, Corr, RT] = minSearch(pars1, data, fitWithin);




end

function pars = gridSearch(data,useParfor,fitWithin)

if fitWithin
    wrapFunc = @raceModelWrapperWithin;
else
    wrapFunc = @raceModelWrapper;
end

%% grid search pars
ta = 0:6;
da = -2:4;
% b = 0;
ndt = 0:.5:2;
thresh = 0:8;

parsMat = CombVec(ta,da,ndt)';

parsMat(parsMat(:,2)>parsMat(:,1),:) = [];
nSteps = length(thresh);

%% search

%split into chunks and parfor
like = cell(nSteps,1);
if useParfor
    parfor j=1:nSteps
        disp(j)
        for i = 1:length(parsMat)
            like{j}(i,1) = wrapFunc( [parsMat(i,:) thresh(j)] ,data);
        end
    end
    
    nll = cat(1, like{:});

    parsMat = CombVec(parsMat',thresh)';
    
else % or just run through all
    parsMat = CombVec(parsMat',thresh)';
    for i = 1:length(parsMat)
%         disp(i)
        nll(i,1) = wrapFunc( repmat(parsMat(i,:),8,1) ,data);
    end
    
end
%% find min



[m,i] = min(nll);

pars = parsMat(i,:);

end


function [nll,x, Corr, RT] = minSearch(pars,data,fitWithin)
%% fit
if fitWithin
    wrapFunc = @raceModelWrapperWithin;
else
    wrapFunc = @raceModelWrapper;
end

[x,fval] = fminsearch(wrapFunc, pars,[],data);

%%


[nll,~, Corr, RT] = wrapFunc(x,data);

end

