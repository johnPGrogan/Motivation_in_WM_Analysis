function output = ArrowRaceModelFitting(dataStruct)
% output = ArrowRaceModelFitting(rt, modelStruct, pars, factorLabels, condCols)
% fit the LBA race model to data (pooled across participants, splitting by
% condition), and draw some figures

struct2workspace(dataStruct);

if ~exist('splitBy','var')
    splitBy = 'trialTypes';
end

[nPP, nConds, nTrials, nRTs] = size(rt);

%% use model to classify each trial as target/non

if ~exist('model','var')
    model = SwapModel();
end

pSep = [];
for i = 1:nPP
    pSep2 = [];
    modelStruct(i).data.(splitBy) = modelStruct(i).(splitBy); % add splitby to data
    datasets = SplitDataByField(modelStruct(i).data, splitBy); % split by that 
    for j = 1:nConds
        [~,pSep1] = model.pdf(datasets{j}, pars(i,1,j), pars(i,2,j), pars(i,3,j));
        pSep2 = nancat(3, pSep2, pSep1);
    end
    pSep = nancat(4, pSep, pSep2);
end

pSep = permute(pSep, [4,3,2,1]);



%% split RTs by target/non - into correct/incorrect rts

[~,maxInd] = max(pSep,[],3);

respType = sq(maxInd);
respType(maxInd>2) = 3;% store resp as 1=targ, 2=guess, 3=misb


maxInd = (sq(maxInd)==1) + 1;

if size(rt,2) > 2
    rtSplit = groupMeans(rt, 3, repmat(maxInd,1,1,1,nRTs), 'dim'); %[pp conds correct iRT trial]
    rtSplit2 = groupMeans(rt, 3, repmat(respType,1,1,1,nRTs), 'dim'); %[pp conds respType iRT trial]
    if nRTs > 1
        rtSplit = permute(rtSplit, [5 2 3 1 4]);%[trial cond correct pp iRT]
        rtSplit2 = permute(rtSplit2, [5 2 3 1 4]);%[trial cond respType pp iRT]
    else
        rtSplit = permute(rtSplit, [4 1 2 3]);
        rtSplit2 = permute(rtSplit2, [4 1 2 3]);
    end
else
    rtSplit = groupMeans(rt, 1, repmat(maxInd,1,1,1,nRTs), 'dim');%[corr cond tr iRT pp]
    rtSplit2 = groupMeans(rt, 1, repmat(respType,1,1,1,nRTs), 'dim');%[respType cond tr iRT pp]
    if nRTs > 1
        rtSplit = permute(rtSplit, [3 2 1 5 4]);%[trial cond correct pp iRT]
        rtSplit2 = permute(rtSplit2, [3 2 1 5 4]);%[trial cond respType pp iRT]
    else
        rtSplit = permute(rtSplit, [3 2 4 1]);
        rtSplit2 = permute(rtSplit2, [3 2 4 1]);
    end
end

for iRT = 1:nRTs
    for iCond = 1:nConds
        for iPP = 1:nPP
            data{iPP,iCond,iRT}.rt_correct{1} = rmmissing(col(rtSplit(:,iCond,2,iPP,iRT)));
            data{iPP,iCond,iRT}.rt_error{1} = rmmissing(col(rtSplit(:,iCond,1,iPP,iRT)));   
        end
    end
end
%% fit model

useParfor = 0;
[fitPars,nll,Corr,RT] = deal(cell(nPP,nConds));

tic
for iPP = 1:nPP
    disp(iPP)
    for iRT = 1:nRTs
        for iCond = 1:nConds
    %         disp(iCond)
            [fitPars{iPP,iCond,iRT},nll{iPP,iCond,iRT},Corr{iPP,iCond,iRT},RT{iPP,iCond,iRT}] = ...
                raceModelFit(data{iPP,iCond,iRT},useParfor);
        end
    end
end
t=toc

%%

output.RT = RT;
output.Corr = Corr;
output.pars = fitPars;
output.data = data;
output.t = t;
output.nPP = nPP;
output.nConds = nConds;
output.nTrials = nTrials;
output.nll = nll;

output.rtSplit2 = rtSplit2;
end