function stars = p2stars(p, criteria)
% converts p values into cell array of strings of stars. default criteria =
% * = p < .05
% ** = p < .01
% *** = p < .001
% **** = p < .0001
% 
% inputs: 
%   p = matrix of p values, any size
%   criteria = vector of decreasing p value criteria. default is [.05, .01, .001, .0001]
% 
% outputs:
%   stars = cell array matching size of p, with one '*' per criteria passed
% 


if ~exist('criteria','var') || isempty(criteria)
    criteria = [.05, .01, .001, .0001];
end

stars = cell(size(p)); % preset

for i = 1:length(criteria)
    stars(p <= criteria(i)) = {repmat('*', 1, i)};
end
