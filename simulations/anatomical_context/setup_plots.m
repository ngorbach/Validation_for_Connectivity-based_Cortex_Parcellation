function h = setup_plots

%try; close(1); end;
figure(1); clf
set(1, 'Position', [0, 200, 1600, 1000]);
axis tight

h{1} = subplot(3,4,1); h{1}.FontSize = 15; h{1}.Title.String = 'True cortex parcellation'; hold on
h{1}.XTick = []; h{1}.YTick = []; h{1}.Box = 'on';

h{2} = subplot(3,4,2); h{2}.FontSize = 15; h{2}.Title.String = 'GCM cortex parcellation'; hold on
h{2}.XTick = []; h{2}.YTick = []; h{2}.Box = 'on';

h{3} = subplot(3,4,3); h{3}.FontSize = 15; h{3}.Title.String = 'Connectivity matrix';
h{3}.XTick = []; h{3}.YTick = []; h{4}.Box = 'on';

h{4} = subplot(3,4,4); h{4}.FontSize = 15; h{4}.Title.String = 'Dissimilarity matrix';
h{4}.XTick = []; h{4}.YTick = []; h{4}.Box = 'on';

h{8} = subplot(3,4,8); h{8}.FontSize = 15; h{8}.Title.String = 'ERM consistency';
h{8}.XTick = []; h{8}.YTick = []; h{8}.Box = 'on';

h{9} = subplot(3,4,9); h{9}.FontSize = 15; h{9}.Title.String = 'GCM consistency';
h{9}.XTick = []; h{9}.YTick = []; h{9}.Box = 'on';

h{7} = subplot(3,4,7); h{7}.FontSize = 15; h{7}.Title.String = 'Centroid dendrogram';
h{7}.XTick = []; h{7}.YTick = []; h{7}.Box = 'on';

h{5} = subplot(3,4,5); h{5}.FontSize = 15; h{5}.Title.String = 'Generalization capacity'; h{5}.XLabel.String = 'inverse temperature (complexity)'; grid on; hold on
h{5}.YLabel.String = 'bits';

h{6} = subplot(3,4,6); h{6}.FontSize = 15; h{6}.Title.String = 'Information content'; h{6}.XLabel.String = 'number of potential clusters'; grid on; hold on
h{6}.YLabel.String = 'bits';

h{10} = subplot(3,4,10); h{10}.FontSize = 15; h{10}.Title.String = 'BIC and AIC'; h{10}.XLabel.String = 'number of potential clusters'; grid on; hold on

h{11} = subplot(3,4,11); h{11}.FontSize = 15; h{11}.Title.String = 'Log Bayesian evidence'; h{11}.XLabel.String = 'number of potential clusters'; grid on; hold on

hold on;
