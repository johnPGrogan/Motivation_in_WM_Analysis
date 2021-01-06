function [like, likelihoods, Corr, RT] = raceModelWrapper(pars,data)

if any(pars([3 4]) < 0)
    like = 1000000*(1-min(pars(3:4)));
    return
end

ta = pars(1,1);
da = pars(1,2);
% b = pars(1,3);
ndt = pars(1,3);
thresh = pars(1,4);



conditions = {
    %  name, rates, and_group, isCorrect(for group), Bias, Non-decision time, Threshold
    'low', [ ta,da ], [ 1, 2 ], [ 1, 0 ], 0, ndt, thresh;
    };

N = 10000;

[Corr, RT] = sim_race(conditions, N);


num_conds = 1;
limit = 10000;
[likelihoods] = get_neg_log_like(data, Corr, RT, num_conds, limit);

like = likelihoods.neg_log_like;

end