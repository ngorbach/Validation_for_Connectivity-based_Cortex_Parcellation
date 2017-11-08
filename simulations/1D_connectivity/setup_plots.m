function h = setup_plots(diffusion)

figure(1); clf
set(1, 'Position', [0, 200, 1200, 600]);
axis tight

h{1} = subplot(2,2,1); h{1}.FontSize = 15; h{1}.Title.String = 'Generalization capacity'; h{1}.YLabel.String = 'bits';
h{1}.XLabel.String = 'inverse temperature (complexity)'; grid on; hold on

for i = 2:4
    h{i} = subplot(2,2,i); h{i}.FontSize = 15; h{i}.Title.String = ['GCM Gibbs distribution, diffusion = ' num2str(diffusion(i-1))]; hold on
    hold on;
    
    h{i}.YTick = []; h{i}.YLabel.String = 'seed space'; h{i}.YLabel.FontSize = 20;
    h{i}.XTick = [];
    h{i}.Box = 'on';
    h{i}.XLabel.String = 'time'; h{i}.XLabel.FontSize = 20;
end

% yyaxis right; h{2} = gca; h{2}.YTick = []; h{2}.YLabel.String = 'target space'; h{2}.YLabel.FontSize = 20;
% h{2}.YLabel.Color = [0,0,0]; 
% yyaxis left; h{1}.YLabel.Color = [0,0,0];
