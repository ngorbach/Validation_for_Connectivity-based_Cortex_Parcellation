function h = setup_plots

figure(1); clf
set(1, 'Position',  [0, 200, 1600, 800]);
axis tight
h{1} = subplot(2,4,1); h{1}.FontSize = 15; h{1}.Title.String = 'GCM Gibbs Distribution'; hold on

h{4} = subplot(2,4,4); h{4}.FontSize = 15; h{4}.Title.String = 'ERM consistency';
h{4}.XTick = []; h{4}.YTick = []; h{4}.Box = 'on';

h{5} = subplot(2,4,5); h{5}.FontSize = 15; h{5}.Title.String = 'GCM consistency';
h{5}.XTick = []; h{5}.YTick = []; h{5}.Box = 'on';

h{2} = subplot(2,4,2); h{2}.FontSize = 15; h{2}.Title.String = 'Generalization capacity'; 
h{2}.XLabel.String = 'inverse temperature (complexity)'; h{2}.YLabel.String = 'bits'; grid on; hold on

h{3} = subplot(2,4,3); h{3}.FontSize = 15; h{3}.Title.String = 'Information content'; 
h{3}.XLabel.String = 'number of potential clusters'; h{3}.YLabel.String = 'bits'; grid on; hold on
% h{4} = subplot(2,3,4); h{4}.FontSize = 15; h{4}.Title.String = 'Centroids dendrogram'; 
% h{4}.XLabel.String = 'potential cluster index'; hold on

h{7} = subplot(2,4,7); h{7}.FontSize = 15; h{7}.Title.String = 'BIC and AIC'; h{7}.XLabel.String = 'number of clusters'; grid on; hold on
%h{7}.YLabel.String = 'log Bayesian evidence 1';

h{8} = subplot(2,4,8); h{8}.FontSize = 15; h{8}.Title.String = 'Log Bayesian evidence'; h{8}.XLabel.String = 'number of clusters'; 
grid on; hold on
hold on;
