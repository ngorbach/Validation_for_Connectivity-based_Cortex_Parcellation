function h = setup_plots

figure(1); clf
set(1, 'Position', [0, 200, 1700, 1000]);
axis tight
h{1} = subplot(3,4,1); h{1}.FontSize = 15; h{1}.Title.String = 'GCM cortex parcellation'; hold on
h{1}.XTick = []; h{1}.YTick = []; h{1}.Box = 'on';

h{2} = subplot(3,4,2); h{2}.FontSize = 15; h{2}.Title.String = 'Minimal projection of ERM onto dissim. matrix';
h{2}.XTick = []; h{2}.YTick = []; h{2}.Box = 'on'; h{2}.XLabel.String = 'seed voxels'; h{2}.YLabel.String = 'seed voxels';

h{5} = subplot(3,4,5); h{5}.FontSize = 15; h{5}.Title.String = 'Distance between centroids';
h{5}.XTick = []; h{5}.YTick = []; h{5}.Box = 'on';

h{3} = subplot(3,4,3); h{3}.FontSize = 15; h{3}.Title.String = 'Generalization capacity'; h{3}.XLabel.String = 'inverse temperature (complexity)'; grid on; hold on
h{3}.YLabel.String = 'bits';

h{4} = subplot(3,4,4); h{4}.FontSize = 15; h{4}.Title.String = 'Information content'; h{4}.XLabel.String = 'number of potential clusters'; grid on; hold on
h{4}.YLabel.String = 'bits';

h{7} = subplot(3,4,7); h{7}.FontSize = 15; h{7}.Title.String = 'ERM consistency';
h{7}.XTick = []; h{7}.YTick = []; h{7}.Box = 'on'; h{7}.XLabel.String = 'seed voxels'; h{7}.YLabel.String = 'seed voxels';

h{8} = subplot(3,4,8); h{8}.FontSize = 15; h{8}.Title.String = 'GCM consistency';
h{8}.XTick = []; h{8}.YTick = []; h{8}.Box = 'on'; h{8}.XLabel.String = 'seed voxels'; h{8}.YLabel.String = 'seed voxels';

h{9} = subplot(3,4,9); h{9}.FontSize = 15; h{9}.Title.String = 'BIC'; h{9}.XLabel.String = 'number of potential clusters'; grid on; hold on

h{10} = subplot(3,4,10); h{10}.FontSize = 15; h{10}.Title.String = 'Log Bayesian evidence'; h{10}.XLabel.String = 'number of potential clusters'; grid on; hold on

hold on;
