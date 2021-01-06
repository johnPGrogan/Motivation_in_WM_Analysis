function [likelihoods] = get_neg_log_like(data, corr, rt, num_conds, limit)
% calculate negative log likelihoods for data using the LBA model's
% predictions
% 
% Robert Udale & Sanjay G Manohar 2019
% Adapted by John Grogan
%% Put predicted data into structure and censor long RTs.
for i = 1:num_conds
    % Add predictions to a structure.
    preds.rt_correct{i}   = rt(corr(:,i)==1, i);
    preds.rt_error{i}     = rt(corr(:,i)==0, i);
    
    % Remove predicted trials which are nan (i.e. when RT = inf).
    preds.rt_correct{i} (isnan(preds.rt_correct{i})) = [];
    preds.rt_error{i}   (isnan(preds.rt_error{i}))   = [];
    
    % censor long RTs from both SIM and DATA
    data.rt_correct{   i}(data.rt_correct{   i}>limit) = nan;
    data.rt_error{ i}(data.rt_error{ i}>limit) = nan;
    preds.rt_correct{  i}(preds.rt_correct{  i}>limit) = nan;
    preds.rt_error{i}(preds.rt_error{i}>limit) = nan;
end

%% Initialise the likelihoods as nans.
likelihoods_correct = nan(1, num_conds);
likelihoods_error = nan(1, num_conds);

for i = 1 : num_conds
    
    %% Calculate probabilities of data being above or below censor threshold.
    
    n_data_rt_correct_toolong = sum(isnan(data.rt_correct{i})); % number of data trials that are too long
    n_data_rt_error_toolong   = sum(isnan(data.rt_error{i})); % number of data trials that are too long
    
    data.rt_correct{i}(isnan(data.rt_correct{i})) = [];         % Remove too-long correct trials.
    data.rt_error{i}(isnan(data.rt_error{i})) = [];             % Remove too-long incorrect trials.
    
    p_model_rt_correct_toolong = mean( isnan(preds.rt_correct{  i}) ); % proportion of model trials that are too long
    p_model_rt_error_toolong = mean( isnan(preds.rt_error{  i}) ); % proportion of model trials that are too long
    
    p_model_rt_correct = nanmean(corr(:,i)); % modelled probability of any response being correct
    p_model_rt_error = 1-nanmean(corr(:,i)); % modelled probability of any response being correct
    
    %% Calculate likelihoods.
    
    % Only get correct likelihoods if there is correct rts in the data (which there will be).
    if ~isempty (data.rt_correct{i})
        
        % If the model predicts no data (but there is real data), then model is implausible.
        if nansum(preds.rt_correct{i}) == 0
            likelihoods_correct(i) = -inf;
%             correct_likelihoods = inf; 
            correct_likelihoods = nan; 
        else
            
            % Get the correct likelihoods.
            correct_likelihoods = ...
                p_model_rt_correct * ...        % normalise the pdf for correct trials, by the modelled probability of us being in the "correct" distribution
                (1-p_model_rt_correct_toolong) * ... % and by the modelled probability of a correct RT being less than 5
                ksdensity (preds.rt_correct{  i}, data.rt_correct{  i}) ...% and take pdf conditioned on being correct and less than 5.
                +eps;                 % Add a very small number.
            
            % Get the total likelihoods for all correct trials.
            likelihoods_correct (i) =  nansum(log(correct_likelihoods)) ... % combine the log likelihoods
                +   (n_data_rt_correct_toolong+eps) * log(p_model_rt_correct_toolong+eps); % add on the modelled log likelihoods for the trials that were longer than 5.

        end
    else
        correct_likelihoods = nan;
    end
    
    % Only get error likelihoods if there are error RTs in the data (which sometimes happens).
    if ~isempty(data.rt_error{i})
        
        % If the model predicts no data (but there is real data), then model is implausible.
        if nansum(preds.rt_error{i}) == 0
            likelihoods_error(i) = -inf;
            %             error_likelihoods = inf;
            error_likelihoods = nan;
        else
            
            % Get the likelihoods for each of the error trials.
            error_likelihoods = ...
                p_model_rt_error * ...        % normalise the pdf for correct trials, by the modelled probability of us being in the "correct" distribution
                (1-p_model_rt_error_toolong) * ... % and by the modelled probability of an RT being less than 5
                ksdensity (preds.rt_error{  i}, data.rt_error{  i})  ...% and take pdf conditioned on being correct and less than 5.
                + eps ;
            
            % Get the total likelihood for all error trials.
            likelihoods_error(i) = nansum(log(error_likelihoods)) ... % combine the log likelihoods
                +  (n_data_rt_error_toolong+eps) * log(p_model_rt_error_toolong+eps); % add on the modelled log likelihoods for the trials that were longer than 5.
            
        end
       
    else
        error_likelihoods = nan;
    end
    
end

%% Get neg log like.

% Add up all the likelihoods and make them negative (This is the negative log likelihood).
likelihoods.neg_log_like = -nansum([likelihoods_correct, likelihoods_error]);
likelihoods.correct = correct_likelihoods;
likelihoods.error = error_likelihoods;

end

