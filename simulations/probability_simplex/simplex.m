function [x,y] = simplex(h1,data,gibbs_dist)

%%% Simplex (Ternary) plot
%%% Author: Didier Gonze
%%% Created 10/12/2015
%%% Updated 10/12/2015

data = bsxfun(@rdivide,data,sum(data,2));
a = data(:,1); b = data(:,2); c = data(:,3);

colors = distinguishable_colors(size(gibbs_dist,2));
colors = mix_colors(gibbs_dist,colors);

%clf

%%% Axes

xx=[0 1 1/2 0];
yy=[0 0 sqrt(3)/2 0];

plot(h1,xx,yy,'k')

% text(-0.05,0,'A','fontsize',18,'HorizontalAlignment','center')
% text(1.05,0,'B','fontsize',18,'HorizontalAlignment','center')
% text(0.5,0.9,'C','fontsize',18,'HorizontalAlignment','center')


%%% Data

x=(1/2)*(2*b+c)./(a+b+c);
y=(sqrt(3)/2)*c./(a+b+c);

hold on;
%subplot(1,3,1)
scatter(h1,x,y,200,colors,'.')
%scatter(x,y)

h1.XAxis.Visible = 'off';
h1.YAxis.Visible = 'off';
%set(gcf,'Color','w')
h1.Color = 'w';
