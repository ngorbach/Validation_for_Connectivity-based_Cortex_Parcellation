function display_result(info_content,gc,max_gc_idx,inv_temp,BIC,AIC,log_bayes_evidence,...
    gibbs_dist_packed1,gibbs_dist_packed2,data1,k,K)

h = setup_plots;
    
[~,opt_k_idx.PA.hc] = max(info_content.hc);
[~,opt_k_idx.PA.kmeans] = max(info_content.kmeans);
simplex(h{1},data1,gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(:,:,max_gc_idx.hc(opt_k_idx.PA.hc)));
gc_plot = gc.hc{opt_k_idx.PA.hc}; gc_plot(max_gc_idx.hc(opt_k_idx.PA.hc)+find(gc_plot(max_gc_idx.hc(opt_k_idx.PA.hc):end)<-0.1,1,'first')+1:end) = NaN;
hold on; plot(h{2},inv_temp.hc,gc_plot,'LineWidth',2); h{2}.YLim(1) = 0;
plot(h{3},[1:k],info_content.hc(1:k),'LineWidth',2);
hold on; plot(h{3},[1:k],info_content.kmeans(1:k),'LineWidth',2,'Color',[0.8,0.8,0.8]);
h{3}.XTick = [1:2:max(K)]; h{3}.XLim(2) = max(K);

figure(1)
centroids = gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(:,:,max_gc_idx.hc(opt_k_idx.PA.hc))'*data1;
centroids = bsxfun(@rdivide,centroids,sum(centroids,2));
dsim_centroids = pdist2(centroids,centroids);
h{6} = subplot(2,4,6); d = dendrogram(linkage(dsim_centroids));
hold on; plot(h{6},[0,size(dsim_centroids,1)+0.5],[0.1,0.1],'--','LineWidth',2)
h{6}.FontSize = 15; h{6}.Title.String = 'Centroids dendrogram'; h{6}.YLim(1) = 0;
h{6}.XLabel.String = 'potential cluster index'; hold on
centroid_labels = cluster(linkage(dsim_centroids),'cutoff',0.1,'criterion','distance');
h{1}.Title.String = [h{1}.Title.String ', k_{effective} = ' num2str(length(unique(centroid_labels)))];

% ERM clustering consistency for histogram clustering
[~,labels] = max(gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(:,:,end),[],2); [~,sort_idx] = sort(labels);
mask = round(gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(sort_idx,:,end)*gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(sort_idx,:,end)' ...
    +gibbs_dist_packed2.hc{opt_k_idx.PA.hc}(sort_idx,:,end)*gibbs_dist_packed2.hc{opt_k_idx.PA.hc}(sort_idx,:,end)');
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

imagesc(h{4},img);
h{4}.FontSize = 15; h{4}.XTick = []; h{4}.YTick = []; h{4}.Box = 'on';
h{4}.Title.String = ['ERM disagreement, k_{potential} = ' num2str(opt_k_idx.PA.hc)];
h{4}.XLabel.String = 'objects'; h{4}.YLabel.String = 'objects';

% GCM clustering consistency
labels_gcm = zeros(length(labels),1);
for i = centroid_labels'
    idx = find(centroid_labels==i);
    for j = 1:length(idx)
        labels_gcm(labels==idx(j)) = i;
    end
end
[~,sort_idx_gcm] = sort(labels_gcm);
imagesc(h{5},gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(sort_idx_gcm,:,max_gc_idx.hc(opt_k_idx.PA.hc))*gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(sort_idx,:,max_gc_idx.hc(opt_k_idx.PA.hc))' ...
    +gibbs_dist_packed2.hc{opt_k_idx.PA.hc}(sort_idx_gcm,:,max_gc_idx.hc(opt_k_idx.PA.hc))*gibbs_dist_packed2.hc{opt_k_idx.PA.hc}(sort_idx_gcm,:,max_gc_idx.hc(opt_k_idx.PA.hc))');
h{5}.FontSize = 15; h{5}.XTick = []; h{5}.YTick = []; h{5}.Box = 'on';
h{5}.Title.String = ['GCM consistency, k_{effective} = ' num2str(length(unique(centroid_labels)))];
h{5}.XLabel.String = 'objects'; h{5}.YLabel.String = 'objects';

t0 = 0;
for k = unique(labels_gcm)'
    t1 = find(labels_gcm(sort_idx_gcm)==k,1,'last')-t0-1;
    hold on; rectangle(h{5},'Position',[t0+1,t0+1,t1,t1],'LineWidth',2,'EdgeColor',[1,1,1]);
    t0 = t0+t1;
end

plot(h{7},K,BIC(K),'LineWidth',2); hold on;
h{7}.XTick = [2:2:max(K)]; h{7}.XLim = [2,max(K)];
hold on; plot(h{7},K,AIC(K),'LineWidth',2);
subplot(2,4,7)
%legend('BIC','AIC','Location','southwest')
%legend('boxoff')
 
subplot(2,4,8);
plot(h{8},K,log_bayes_evidence{1}(K),'LineWidth',2);
yyaxis right
hold on; plot(h{8},K,log_bayes_evidence{2}(K),'LineWidth',2);
h{8}.XTick = [2:2:max(K)]; h{8}.XLim = [2,max(K)];
%legend('BIC','hyperparam 1','hyperparam 2','Location','southwest')
%legend('hyperparameter 1','hyperparameter 2','Location','southeast')
%legend('boxoff')
drawnow


indexmax = find(max(info_content.hc) == info_content.hc);
x = 1:max(K);
xmax = x(indexmax);
ymax = info_content.hc(indexmax);
strmax = ['k_{effective} = ' num2str(length(unique(centroid_labels)))];
%text(h{3},xmax,ymax-0.15,strmax,'HorizontalAlignment','left','FontSize',15);
        
drawnow

% display runtime
disp(['runtime = ' num2str(toc/60) ' minutes']);

figure(2); set(2, 'Position',  [0, 200, 1600, 300]);
h2 = subplot(1,6,1);
simplex(h2,data1,gibbs_dist_packed1.hc{opt_k_idx.PA.hc}(:,:,max_gc_idx.hc(opt_k_idx.PA.hc))); h2.FontSize = 15; h2.Title.String = 'PA, HC';
h2 = subplot(1,6,2);
simplex(h2,data1,gibbs_dist_packed1.kmeans{opt_k_idx.PA.kmeans}(:,:,max_gc_idx.kmeans(opt_k_idx.PA.kmeans))); h2.FontSize = 15; h2.Title.String = 'PA, kmeans';
h2 = subplot(1,6,3);
BIC(1) = NaN; [~,opt_k_idx.BIC] = min(BIC);
simplex(h2,data1,round(gibbs_dist_packed1.hc{opt_k_idx.BIC}(:,:,end))); h2.FontSize = 15; h2.Title.String = 'BIC';
h2 = subplot(1,6,4);
AIC(1) = NaN; [~,opt_k_idx.AIC] = min(AIC);
simplex(h2,data1,round(gibbs_dist_packed1.hc{opt_k_idx.AIC}(:,:,end))); h2.FontSize = 15; h2.Title.String = 'AIC';
h2 = subplot(1,6,5);
log_bayes_evidence{1}(1) = NaN;
[~,opt_k_idx.log_bayes_evidence1] = max(log_bayes_evidence{1});
simplex(h2,data1,round(gibbs_dist_packed1.hc{opt_k_idx.log_bayes_evidence1}(:,:,end)));  h2.FontSize = 15; h2.Title.String = 'Bayes evidence 1';
h2 = subplot(1,6,6);
log_bayes_evidence{2}(1) = NaN;
[~,opt_k_idx.log_bayes_evidence2] = max(log_bayes_evidence{2});
simplex(h2,data1,round(gibbs_dist_packed1.hc{opt_k_idx.log_bayes_evidence2}(:,:,end))); h2.FontSize = 15; h2.Title.String = 'Bayes evidence 2';

% ERM clustering consistency for kmeans
[~,labels] = max(gibbs_dist_packed1.kmeans{opt_k_idx.PA.hc}(:,:,end),[],2); [~,sort_idx] = sort(labels);
mask = round(gibbs_dist_packed1.kmeans{opt_k_idx.PA.hc}(sort_idx,:,end)*gibbs_dist_packed1.kmeans{opt_k_idx.PA.hc}(sort_idx,:,end)' ...
    +gibbs_dist_packed2.kmeans{opt_k_idx.PA.hc}(sort_idx,:,end)*gibbs_dist_packed2.kmeans{opt_k_idx.PA.hc}(sort_idx,:,end)');
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


figure; imagesc(img); h2 = gca;
h2.FontSize = 15; h2.XTick = []; h2.YTick = []; h2.Box = 'on';
h2.Title.String = ['kmeans ERM disagreement, k_{potential} = ' num2str(opt_k_idx.PA.hc)];
h2.XLabel.String = 'objects'; h2.YLabel.String = 'objects';
snapnow