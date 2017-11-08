function display_result(gc,info_content,gibbs_dist_packed1,gibbs_dist_packed2,inv_temp,log_bayes_evidence,BIC,AIC,K,dsim,centroids,number_misclustered_objects,seed_coord)

 disp(['Number of potential clusters: ' num2str(K(end))])
 
h = setup_plots;

%~,max_k_idx] = max(info_content);
max_k_idx = length(info_content);
[~,max_gc_idx] = max(gc{max_k_idx});
gibbs_dist1 = gibbs_dist_packed1{max_k_idx}(:,:,max_gc_idx-1);
colors = distinguishable_colors(size(gibbs_dist1,2));
colors = mix_colors(gibbs_dist1,colors);
scatter3(h{1},seed_coord(:,1),seed_coord(:,2),seed_coord(:,3),5,colors,'filled');


[~,seed_labels] = max(gibbs_dist1,[],2); 

%cluster seed labels
if size(centroids,1)<20
    Kmax = size(centroids,1);
else
    Kmax = 20;
end


plot(h{5},[1,K],number_misclustered_objects([1,K])./length(seed_labels),'LineWidth',2)
h{5}.Box = 'on';  h{5}.FontSize = 15; h{5}.YGrid = 'on'; h{5}.XGrid = 'on';
h{5}.Title.String = '% of misclustered seed voxels by the ERM'; h{5}.XLabel.String = 'number of potential clusters';

% dsim_centroids = -centroids * log(centroids)';
% dsim_centroids = (dsim_centroids + dsim_centroids') / 2;
% centroid_coords = bsxfun(@rdivide,gibbs_dist_packed1{k}(:,:,max_gc_idx)' * seed_coord,...
%     sum(gibbs_dist_packed1{k}(:,:,max_gc_idx),1)');
%     
% % sort_idx_gcm = symrcm(dsim_centroids);
% % dsim_centroids = dsim_centroids(sort_idx_gcm,sort_idx_gcm);
% % imagesc(h{3},dsim_centroids);
% % h{3}.Title.String = 'Centroids dissimilarity matrix';
% 
% [~,sort_idx_gcm] = sort(centroid_coords(:,1));
% dsim_centroids = dsim_centroids(sort_idx_gcm,sort_idx_gcm);
% tree = linkage(dsim_centroids);
% subplot(2,3,3); dendrogram(tree);
% h{3}.Title.String = 'Centroids tree';

[~,max_gc_idx] = max(gc{max_k_idx});
gc_plot = gc{max_k_idx}; gc_plot(max_gc_idx+find(gc_plot(max_gc_idx:end)<-0.1,1,'first')+1:end) = NaN;
hold on; plot(h{3},inv_temp,gc_plot,'LineWidth',2);
h{3}.YLim(1) = 0;

plot(h{4},[1,K],info_content([1,K]),'LineWidth',2);
h{4}.XLim(1) = 1;


plot(h{9},K,BIC(K),'LineWidth',2); 
hold on; plot(h{9},K,AIC(K),'LineWidth',2); 
subplot(3,4,9);
%legend('BIC','AIC','Location','southwest'); legend('boxoff')

subplot(3,4,10);
plot(h{10},K,log_bayes_evidence{1}(K),'LineWidth',2);
hold on; plot(h{10},K,log_bayes_evidence{2}(K),'LineWidth',2);
%legend('hyperparameter 1','hyperparameter 2','Location','southwest'); legend('boxoff')

figure(1)
%dsim_centroids = pdist2(centroids,centroids);
for i = 1:size(centroids,1)
    for j = i+1:size(centroids,1)
        dsim_centroids(i,j) = JSDiv(centroids(i,:),centroids(j,:));
        dsim_centroids(j,i) = dsim_centroids(i,j);
    end
end

h{6} = subplot(3,4,6); d = dendrogram(linkage(dsim_centroids));
hold on; plot(h{6},[0,size(dsim_centroids,1)+0.5],[0.1,0.1],'--','LineWidth',2)
h{6}.FontSize = 15; h{6}.Title.String = 'Centroids dendrogram'; 
h{6}.XTick = [];%h{4}.XLabel.String = 'potential cluster index'; hold on
h{6}.YLim(1) = 0;

centroid_labels = cluster(linkage(dsim_centroids),'cutoff',0.1,'criterion','distance');
centroid_labels = kmeans(centroids,length(unique(centroid_labels)));
h{1}.Title.String = [h{1}.Title.String ', k_{effective} = ' num2str(length(unique(centroid_labels)))];



seed_labels_gcm = zeros(length(seed_labels),1);
for i = centroid_labels'
    idx = find(centroid_labels==i);
    for j = 1:length(idx)
        seed_labels_gcm(seed_labels==idx(j)) = i;
    end
end

[~,sort_idx_gcm] = sort(seed_labels_gcm);
imagesc(h{2},dsim(sort_idx_gcm,sort_idx_gcm));

t0 = 0;
for k = unique(seed_labels_gcm)'
    t1 = find(seed_labels_gcm(sort_idx_gcm)==k,1,'last')-t0-1;
    hold on; rectangle(h{2},'Position',[t0+1,t0+1,t1,t1],'LineWidth',1,'EdgeColor',[0.2,0.2,0.2]);
    t0 = t0+t1;
end
h{2}.XTick = []; h{2}.YTick = []; h{2}.Box = 'on'; 
h{2}.Title.String = ['Min proj. of GCM onto dissim. matrix'];
h{2}.FontSize = 15;

% ERM clustering consistency
%[~,sort_idx_gcm] = sort(seed_labels_gcm);
[~,sort_idx_gcm] = sort(seed_labels);

mask = round(gibbs_dist_packed1{max_k_idx}(sort_idx_gcm,:,end)*gibbs_dist_packed1{max_k_idx}(sort_idx_gcm,:,end)' ...
    +gibbs_dist_packed2{max_k_idx}(sort_idx_gcm,:,end)*gibbs_dist_packed2{max_k_idx}(sort_idx_gcm,:,end)');
img = cat(3,zeros(size(mask),'single'),zeros(size(mask),'single'),zeros(size(mask),'single'));

zeros_mask = logical(zeros(size(mask)));
img(logical(cat(3,mask==0,zeros_mask,zeros_mask))) = 55/255;
img(logical(cat(3,zeros_mask,mask==0,zeros_mask))) = 126/255;
img(logical(cat(3,zeros_mask,zeros_mask,mask==0))) = 184/255;

img(logical(cat(3,mask==1,zeros_mask,zeros_mask))) = 228/255;
img(logical(cat(3,zeros_mask,mask==1,zeros_mask))) = 26/255;
img(logical(cat(3,zeros_mask,zeros_mask,mask==1))) = 28/255;

img(logical(cat(3,mask==2,zeros_mask,zeros_mask))) = 1;
img(logical(cat(3,zeros_mask,mask==2,zeros_mask))) = 1;
img(logical(cat(3,zeros_mask,zeros_mask,mask==2))) = 51/255;

imagesc(h{7},img);
h{7}.FontSize = 15; h{7}.XTick = []; h{7}.YTick = []; h{7}.Box = 'on';
h{7}.Title.String = ['ERM disagreement, k_{potential} = ' num2str(K(end))];
h{7}.XLabel.String = 'seed voxels'; h{7}.YLabel.String = 'seed voxels';

% t0 = 0;
% for k = unique(seed_labels_gcm)'
%     t1 = find(seed_labels_gcm(sort_idx_gcm)==k,1,'last')-t0-1;
%     hold on; rectangle(h{7},'Position',[t0+1,t0+1,t1,t1],'LineWidth',2,'EdgeColor',[1,1,1]);
%     t0 = t0+t1;
% end

% GCM clustering consistency
imagesc(h{8},gibbs_dist_packed1{max_k_idx}(sort_idx_gcm,:,max_gc_idx)*gibbs_dist_packed1{max_k_idx}(sort_idx_gcm,:,max_gc_idx)' ...
    +gibbs_dist_packed2{max_k_idx}(sort_idx_gcm,:,max_gc_idx)*gibbs_dist_packed2{max_k_idx}(sort_idx_gcm,:,max_gc_idx)');
h{8}.FontSize = 15; h{8}.XTick = []; h{8}.YTick = []; h{8}.Box = 'on';
h{8}.Title.String = ['GCM consistency, k_{effective} = ' num2str(length(unique(centroid_labels)))];
h{8}.XLabel.String = 'seed voxels'; h{8}.YLabel.String = 'seed voxels';

% t0 = 0;
% for k = unique(seed_labels_gcm)'
%     t1 = find(seed_labels_gcm(sort_idx_gcm)==k,1,'last')-t0-1;
%     hold on; rectangle(h{8},'Position',[t0+1,t0+1,t1,t1],'LineWidth',2,'EdgeColor',[1,1,1]);
%     t0 = t0+t1;
% end
%     
drawnow