function SuperTitle(titleText)
% puts one title at the top of a subplot

g = gca; % store current axes

set(gcf,'NextPlot','add');
axes; % make new

h = title(titleText); % put title
set(gca,'Visible','off'); % hide current
set(h,'Visible','on'); % show title

axes(g); % restore previous axes
end