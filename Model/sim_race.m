function [Corr, RT] = sim_race(conditions, N)
% Simulate race model for multiple conditions
% 
% Robert Udale & Sanjay G Manohar, 2019
% Adapted by John Grogan


for c = 1:size(conditions,1)
    
    rng(1);
    thr = conditions{c, 7};  % Threshold.
    base_mu   = conditions{c,2}; % randomise rate of accumulation across trials
    bias = conditions{c, 5};
    nd_time = conditions {c, 6}; % Non-decision time.
    
    mu = repmat(bias + base_mu, N, 1) + randn(N,size(base_mu, 2));
    %     mu = exp(repmat(base_mu, N, 1) + randn(N,size(base_mu, 2)));
%     mu = repmat(bias + base_mu, N, 1) + exp(randn(N,size(base_mu, 2)));
    
    
    
    time      = (thr./mu) ;    % when does each unit reach threshold?
    time(time<0) = inf;
    [rt_sort,ord]   = sort(time,2); % put events in order
    % ord now contains the event index happening at each time
    and_group  = conditions{c,3}';
    NG         = max(and_group);  % number of groups
    and_values = false(N,NG); % matrix showing if each group is satisfied
    % get number of triggers required for each group
    clear group_size
    for g=1:NG
        group_size(g) = sum(and_group==g);
    end
    winning_group = nan(N,1); % for each trial, find the winning group.
    rt            = nan(N,1); % and RT
    % step through each event occurring (in temporal order)
    for i=1:size(ord,2)  % events
        % which and-groups have been triggered so far?
        groups_count_so_far = and_group( ord(:,1:i) );
        % this will contain which group has been triggered
        group_triggered = nan(N,1);
        
        for g=1:NG % for each group,
            % how many triggers for this group so far?
            n_triggers_for_group = sum( groups_count_so_far == g, 2);
            % has the group therefore triggered? if so, put the group number into
            % 'group_triggered'.
            group_triggered( n_triggers_for_group >= group_size(g) ) = g;
        end
        % vector of trials that have finished at this event, i.e.
        % does not yet have a winning_group, but now has a group_triggered
        % value indicating which group.
        finished_now = isnan(winning_group) & ~isnan(group_triggered);
        % set the winning group for those trials, to the group number that
        % finished.
        winning_group(finished_now) = group_triggered(finished_now);
        % select the current event's column, and put the times for the finished
        % trials into winning_rt.
        rt(   finished_now) = rt_sort(finished_now, i) ;
    end
    resp = winning_group;
    if false
        [rt,resp] = min(time,[],2); % rt = find the smallest time.
    end
    isCorr    = conditions{c,4}(resp)'; % was the winning response 'correct'?
    
    rt = rt + nd_time; 
    rt = rt * 1000;
    
    isCorr(isinf(rt)) = nan;
    rt(isinf(rt)) = nan;
        
    RT(:,c)   = rt;
    Corr(:,c) = isCorr;
    
end



end

