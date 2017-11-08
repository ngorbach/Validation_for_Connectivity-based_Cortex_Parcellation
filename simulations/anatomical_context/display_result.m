function [ktrue2,centroid_labels] = display_result(info_content,gc,inv_temp,BIC,AIC,log_bayes_evidence,...
    gibbs_dist_packed1,gibbs_dist_packed2,K,p2,h,seed_coords,data,fiber_assignment) 

[~,max_k_idx] = max(info_content);
    [~,max_gc_idx] = max(gc{max_k_idx});
    gibbs_dist1 = gibbs_dist_packed1{max_k_idx}(:,:,max_gc_idx-1);
    %centroid_coords = princomp((gibbs_dist1' * seed_coords)');
    %[~,sort_idx] = sort(centroid_coords(:,1),'descend');
    %gibbs_dist1 = gibbs_dist1(:,sort_idx);
    colors = distinguishable_colors(size(gibbs_dist1,2));
    gibbs_dist_colors = mix_colors(gibbs_dist1,colors);
    for i = 1:length(p2)
        p2{i}.Color = gibbs_dist_colors(i,:);
        p2{i}.LineWidth = 1;
    end
    scatter(h{2},seed_coords(:,1),seed_coords(:,2),50,gibbs_dist_colors,'filled')
    
       
    [~,max_gc_idx] = max(gc{max_k_idx});
    gc_plot = gc{max_k_idx}; gc_plot(max_gc_idx+find(gc_plot(max_gc_idx:end)<-0.1,1,'first')+1:end) = NaN;
    plot(h{5},inv_temp,gc_plot,'LineWidth',2); h{5}.YLim = [0,ceil(max(gc_plot))];
    
    plot(h{6},1:max(K),info_content,'LineWidth',2); h{6}.XLim(1) = 1;
    h{6}.XTick = [1:max(K)]; h{6}.XLim(2) = max(K); h{6}.YLim(1) = 0;
    
    subplot(3,4,10);
    plot(h{10},K,BIC(K),'LineWidth',2);
    hold on; plot(h{10},K,AIC(K),'LineWidth',2);
    h{10}.XTick = [2:max(K)]; h{10}.XLim = [2,max(K)];
    %legend('BIC','AIC','Location','north')
    %legend('boxoff')
    
    subplot(3,4,11);
    plot(h{11},K,log_bayes_evidence{1}(K),'LineWidth',2);
    hold on; plot(h{11},K,log_bayes_evidence{2}(K),'LineWidth',2);
    h{11}.XTick = [2:max(K)]; h{11}.XLim = [2,max(K)];
    %legend('hyperparameter 1','hyperparameter 2','Location','north')
    %legend('boxoff')
    
%     centroid_coords = gibbs_dist1' * seed_coords;
%     [~,sort_idx] = sort(centroid_coords(:,1),'ascend');
%     gibbs_dist1 = gibbs_dist1(:,sort_idx);
    
    centroids = gibbs_dist1'*data{1};
    centroids = bsxfun(@rdivide,centroids,sum(centroids,2));
   
%     dsim_centroids = zeros(size(centroids,1));
%     for i = 1:size(centroids,1)
%         for j = i+1:size(centroids,1)
%             centroids1 = centroids(i,:); centroids2 = centroids(j,:);
%             idx = centroids1 + centroids2 == 0;
%             centroids1(idx) = []; centroids2(idx) = [];
%             centroids1(centroids1==0) = eps; centroids2(centroids2==0) = eps;
%             dsim_centroids(i,j) = JSDiv(centroids1,centroids2);
%             dsim_centroids(j,i) = dsim_centroids(i,j);
%         end
%     end

    % cluster centroids
    %eval = evalclusters(centroids,'linkage','CalinskiHarabasz','KList',[2:c]);
    %disp(['number of effective clusters: ' num2str(eval.OptimalK)]);
    %centroid_labels = kmeans(centroids,centroids,eval.OptimalK);
    dsim_centroids = pdist2(centroids,centroids);
%     for i = 1:size(centroids,1)
%         for j = i+1:size(centroids,1)
%             dsim_centroids(i,j) = JSDiv(centroids(i,:),centroids(j,:));
%             dsim_centroids(j,i) = dsim_centroids(i,j);
%         end
%     end
    centroid_labels = cluster(linkage(dsim_centroids),'cutoff',0.01,'criterion','distance');
    centroid_labels = kmeans(centroids,length(unique(centroid_labels)));
    %disp(['number of effective clusters: ' num2str(length(unique(centroid_labels)))]);
    
    %[~,sort_idx] = sort(centroid_labels);
    %imagesc(h{5},dsim_centroids(sort_idx,sort_idx)); 
    subplot(3,4,7); d = dendrogram(linkage(dsim_centroids));
    h{7}.FontSize = 15; h{7}.XTick = []; h{7}.YLim(1) = 0;
    hold on; plot(h{7},[0,size(dsim_centroids,1)+0.5],[0.01,0.01],'--','LineWidth',2)
    %h{5}.XTick = [1:size(dsim_centroids,1)]; h{5}.YTick = [1:size(dsim_centroids,1)];
    %h{5}.XTick = []; h{5}.YTick = []; h{5}.Box = 'on';
    h{7}.Title.String = 'Centroids dendrogram'; %colormap(h{4},gray)
    %h{7}.XLabel.String = 'potential cluster index';
    
    
    [~,seed_labels] = max(gibbs_dist_packed1{max_k_idx}(:,:,end),[],2);
    if length(unique(centroid_labels)) < length(unique(seed_labels))
        seed_labels_gcm = zeros(length(seed_labels),1);
        for i = centroid_labels'
            idx = find(centroid_labels==i);
            for j = 1:length(idx)
                seed_labels_gcm(seed_labels==idx(j)) = i;
            end
        end
    else
        seed_labels_gcm = seed_labels;
    end
   
    
     % ERM clustering consistency
     mask = round(gibbs_dist_packed1{max_k_idx}(:,:,end)*gibbs_dist_packed1{max_k_idx}(:,:,end)' ...
         +gibbs_dist_packed2{max_k_idx}(:,:,end)*gibbs_dist_packed2{max_k_idx}(:,:,end)');
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
     
     imagesc(h{8},flipdim(flipdim(img,1),2))
    h{8}.FontSize = 15; h{8}.XTick = []; h{8}.YTick = []; h{8}.Box = 'on';
    h{8}.Title.String = ['ERM disagreement, k_{potential} = ' num2str(K(max_k_idx-1))];
    h{8}.XLabel.String = 'seed voxels'; h{8}.YLabel.String = 'seed voxels';
    
    % GCM clustering consistency
    [~,sort_idx] = sort(seed_labels_gcm);
    imagesc(h{9},gibbs_dist_packed1{max_k_idx}(sort_idx,:,max_gc_idx)*gibbs_dist_packed1{max_k_idx}(sort_idx,:,max_gc_idx)' ...
        +gibbs_dist_packed2{max_k_idx}(sort_idx,:,max_gc_idx)*gibbs_dist_packed2{max_k_idx}(sort_idx,:,max_gc_idx)');
    h{9}.FontSize = 15; h{9}.XTick = []; h{9}.YTick = []; h{9}.Box = 'on';
    h{9}.Title.String = ['GCM consistency, k_{effective} = ' num2str(length(unique(centroid_labels)))];
    h{9}.XLabel.String = 'seed voxels'; h{9}.YLabel.String = 'seed voxels';
    
    
    indexmax = find(max(info_content) == info_content);
    x = 1:max(K);
    xmax = x(indexmax);
    ymax = info_content(indexmax);
    strmax = ['k_{effective} = ' num2str(length(unique(centroid_labels)))];
    %text(h{6},xmax,ymax-0.15,strmax,'HorizontalAlignment','left','FontSize',15);
    
    ktrue2 = sum(cellfun(@(x) ~isempty(x),fiber_assignment));
    if sum(cellfun(@(x) (ismember(6,x) || ismember(7,x)) && length(x)==1,fiber_assignment)) == 2
        ktrue2 = ktrue2-1;
    end
    text(h{1},150,450,['k_{true} = ' num2str(ktrue2)],'HorizontalAlignment','left','FontSize',15)
    text(h{2},150,450,['k_{effective} = ' num2str(length(unique(centroid_labels)))],'HorizontalAlignment','left','FontSize',15)