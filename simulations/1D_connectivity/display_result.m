function display_result(gc,inv_temp,gibbs_dist_packed1,gibbs_dist_packed2,diffusion,...
    seed_coords,data1,trajectory_all,time,connectivity_matrix,rem_idx)

h = setup_plots(diffusion);
gc{1}(gc{1}<0) = 0; gc{2}(gc{1}<0) = 0; gc{3}(gc{3}<0) = 0;
plot(h{1},inv_temp,gc{1},'k-','LineWidth',2); plot(h{1},inv_temp,gc{2},'k-.','LineWidth',2); plot(h{1},inv_temp,gc{3},'k.','LineWidth',2);
for n = 1:length(diffusion)
    gc{n}(gc{n}<0) = 0;
    %plot(h{1},inv_temp,gc{n},'LineWidth',2);
    
    [~,max_gc_idx] = max(gc{n});
    gibbs_dist1 = gibbs_dist_packed1{n}(:,:,max_gc_idx-1);
    
%     lin_idx = sub2ind(size(gibbs_dist1),[1:size(gibbs_dist1,1)]',true_cluster_labels(:,n));
%     gibbs_dist1 = zeros(size(gibbs_dist1));
%     gibbs_dist1(lin_idx) = 1;
 
    % cluster coordinates
    cluster_coords = gibbs_dist1' * seed_coords' ...
        ./ sum(gibbs_dist1,1)';
    [~,sort_idx] = sort(cluster_coords);
    gibbs_dist1 = gibbs_dist1(:,sort_idx);

    centroids{n} = gibbs_dist1'*data1;
    centroids{n} = bsxfun(@rdivide,centroids{n},sum(centroids{n},2));
    
    colors = distinguishable_colors(size(gibbs_dist1,2));
    colors = mix_colors(gibbs_dist1,colors);
    
    plot_idx = randperm(size(trajectory_all{n},1)); plot_idx = plot_idx(1:50);
    p = plot(h{n+1},time,trajectory_all{n}(plot_idx,2:end),'LineWidth',1);
    y = reshape([p.YData],length(p(1).YData),[]); y = y(1,:);
    for i = 1:length(p)
        p(i).Color = colors(y(i),:);
    end
    h{n+1}.YLim = h{2}.YLim;
    subplot(2,2,n+1);
    yyaxis right; hyy = gca; hyy.YTick = []; hyy.YLabel.String = 'target space'; hyy.YLabel.FontSize = 20;
    hyy.YLabel.Color = [0,0,0];
    yyaxis left; h{1}.YLabel.Color = [0,0,0];
end
h{1}.XLim(2) = 3;
%legend(h{1},['SDE diffusion = ' num2str(diffusion(1))],['SDE diffusion = ' num2str(diffusion(2))],['SDE diffusion = ' num2str(diffusion(3))],'Location','SouthEast','Box','off')
legend(h{1},'low stochasticity','medium stochasticity','high stochasticity','Location','SouthEast','Box','off')

figure(3); clf 
set(3, 'Position', [0, 200, 1000, 400]); axis tight
for n = 1:3
    c = zeros(size(connectivity_matrix{n}{1},1),length(rem_idx(n,:)));
    c(:,~rem_idx(n,:)) = connectivity_matrix{n}{1};
    c = sum(c,1);
    s{n} = subplot(1,3,n);
    c2 = zeros(1,length(c)); c2(1:end) = c(1:end); a = area(c2); a.FaceColor =[117,112,179]./255; hold on;
%     c2 = zeros(1,length(c)); c2(12:35) = c(12:35); a = area(c2); a.FaceColor = [117,112,179]./255; hold on;
%     c2 = zeros(1,length(c)); c2(36:end) = c(36:end); a = area(c2); a.FaceColor = [117,112,179]./255; hold on;
   %plot(c,'LineWidth',2); hold on;
    s{n}.XLim = [0,60]; s{n}.XTick = []; s{n}.YTick = []; s{n}.Title.String = ['SDE Diffusion = ' num2str(diffusion(n))]; s{n}.Title.FontSize = 15;
end
s{1}.Title.String = 'low stochasticity';
s{2}.Title.String = 'med stochasticity';
s{3}.Title.String = 'high stochasticity';

%% Centroid dendrogram

figure(4); clf 
set(4, 'Position', [0, 200, 1000, 400]); axis tight
for n = 1:length(diffusion)
    s2{n} = subplot(1,3,n);
    %dsim_centroids = pdist2(centroids{n},centroids{n});
    for i = 1:size(centroids{n},1)
        for j = i+1:size(centroids{n},1)
            dsim_centroids(i,j) = JSDiv(centroids{n}(i,:),centroids{n}(j,:));
            dsim_centroids(j,i) = dsim_centroids(i,j);
        end
    end
    d = dendrogram(linkage(dsim_centroids));
    hold on; plot(s2{n},[0,size(dsim_centroids,1)+0.5],[0.08,0.08],'--','LineWidth',2)
    s2{n}.FontSize = 15;
    s2{n}.YLim(1) = 0; 
    s2{n}.Title.String = ['SDE Diffusion = ' num2str(diffusion(n))];
    centroid_labels{n} = cluster(linkage(dsim_centroids),'cutoff',0.08,'criterion','distance');
    h{n+1}.Title.String = [h{n+1}.Title.String ', k_{effective} = ' num2str(length(unique(centroid_labels{n})))];
    
    indexmax = find(max(gc{n}) == gc{n});
    xmax = inv_temp(indexmax);
    ymax = gc{n}(indexmax);
    strmax = ['k_{effective} = ' num2str(length(unique(centroid_labels{n})))];
    text(h{1},xmax,ymax,strmax,'HorizontalAlignment','left','FontSize',15);
end
s2{1}.Title.String = 'low stochasticity';
s2{2}.Title.String = 'med stochasticity';
s2{3}.Title.String = 'high stochasticity';
%% ERM and GCM consistency
figure(5); clf 
set(5, 'Position', [0, 200, 1200, 800]); axis tight

for n = 1:length(diffusion)
    
    % ERM clustering consistency
    h{n}=subplot(2,3,n);
    [~,seed_labels] = max(gibbs_dist_packed1{n}(:,:,end),[],2);
    lin_idx = sub2ind(size(gibbs_dist1),[1:size(gibbs_dist1,1)]',seed_labels);
    gibbs_dist1 = zeros(size(gibbs_dist1));
    gibbs_dist1(lin_idx) = 1;
    
    mask = round(gibbs_dist1*gibbs_dist1' ...
        +gibbs_dist_packed2{n}(:,:,end)*gibbs_dist_packed2{n}(:,:,end)');
    img = cat(3,zeros(size(mask)),zeros(size(mask)),zeros(size(mask)));
    
    img(logical(cat(3,mask==0,zeros(size(mask)),zeros(size(mask))))) = 55/255;
    img(logical(cat(3,zeros(size(mask)),mask==0,zeros(size(mask))))) = 126/255;
    img(logical(cat(3,zeros(size(mask)),zeros(size(mask)),mask==0))) = 184/255;
    
    img(logical(cat(3,mask==1,zeros(size(mask)),zeros(size(mask))))) = 228/255;
    img(logical(cat(3,zeros(size(mask)),mask==1,zeros(size(mask))))) = 26/255;
    img(logical(cat(3,zeros(size(mask)),zeros(size(mask)),mask==1))) = 28/255;
    
    img(logical(cat(3,mask==2,zeros(size(mask)),zeros(size(mask))))) = 1;
    img(logical(cat(3,zeros(size(mask)),mask==2,zeros(size(mask))))) = 1;
    img(logical(cat(3,zeros(size(mask)),zeros(size(mask)),mask==2))) = 51/255;
    
    imagesc(h{n},img)
    
%     imagesc(gibbs_dist_packed1{n}(:,:,end)*gibbs_dist_packed1{n}(:,:,end)' ...
%         +gibbs_dist_packed2{n}(:,:,end)*gibbs_dist_packed2{n}(:,:,end)');
    h{n}.FontSize = 15; h{n}.XTick = []; h{n}.YTick = []; h{n}.Box = 'on';
    h{n}.Title.String = ['ERM disagreement, k_{potential} = ' num2str(length(unique(seed_labels)))];
    h{n}.XLabel.String = 'seed voxels'; h{n}.YLabel.String = 'seed voxels';
    
    % GCM clustering consistency
    
    if length(unique(centroid_labels{n})) < length(unique(seed_labels))
        seed_labels_gcm = zeros(length(seed_labels),1);
        for i = centroid_labels{n}'
            idx = find(centroid_labels{n}==i);
            for j = 1:length(idx)
                seed_labels_gcm(seed_labels==idx(j)) = i;
            end
        end
    else
        seed_labels_gcm = seed_labels;
    end
    h{n+3} = subplot(2,3,n+3);
    [~,sort_idx] = sort(seed_labels_gcm);
    imagesc(h{n+3},gibbs_dist_packed1{n}(sort_idx,:,max_gc_idx)*gibbs_dist_packed1{n}(sort_idx,:,max_gc_idx)' ...
        +gibbs_dist_packed2{n}(sort_idx,:,max_gc_idx)*gibbs_dist_packed2{n}(sort_idx,:,max_gc_idx)');
    h{n+3}.FontSize = 15; h{n+3}.XTick = []; h{n+3}.YTick = []; h{n+3}.Box = 'on';
    h{n+3}.Title.String = ['GCM consistency, k_{effective} = ' num2str(length(unique(centroid_labels{n})))];
    h{n+3}.XLabel.String = 'seed voxels'; h{n+3}.YLabel.String = 'seed voxels';
end