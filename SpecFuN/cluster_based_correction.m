% cluster-based correction 
% within-subject


n_permutes = 1000;
num_1=size(X,2);
num_2=size(X,3);

voxel_pval   = 0.001;
cluster_pval = 0.005;

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_1);
permuted_vals    = zeros(n_permutes,num_1,num_2);
max_clust_info   = zeros(n_permutes,1);
nsub=size(X,1);


% compute actual t-test of difference
[h,p,ci,stat]=ttest(X,Y);
real_t=squeeze(stat.tstat);
realmean=mean(X)-mean(Y);

eegpower=[X;Y];
% generate pixel-specific null hypothesis parameter distributions
clear permuted_tvals max_pixel_pvals
for permi = 1:n_permutes
    
    
    fake_condition_mapping = randsample(1:nsub*2,nsub*2);
    
    
    % compute t-map of null hypothesis
    [h,p,ci,stat]=ttest(eegpower(fake_condition_mapping(1:nsub),:,:),eegpower(fake_condition_mapping(nsub+1:end),:,:));
    tmap   = squeeze(stat.tstat);
    
    % save all permuted values
    permuted_tvals(permi,:,:) = tmap;
    
    % save maximum pixel values
    max_pixel_pvals(permi,:) = [ min(tmap(:)) max(tmap(:)) ];
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,nsub-1))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(tmap);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % notes: cellfun is superfast, and the zero accounts for empty maps
end


% now compute Z-map
zmap = (real_t-squeeze(mean(permuted_tvals,1)))./squeeze(std(permuted_tvals));

figure
subplot(231)
contourf(xaxislabel,yaxislabel,squeeze(mean(X)),40,'linecolor','none')
axis square
set(gca,'clim',[-0.05 0.1],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Successful Retrieval')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')


subplot(232)
contourf(xaxislabel,yaxislabel,squeeze(mean(Y)),40,'linecolor','none')
axis square
set(gca,'clim',[-0.05 0.1],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Failed Retrieval')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')


subplot(233)
contourf(xaxislabel,yaxislabel,squeeze(mean(X))-squeeze(mean(Y)),40,'linecolor','none')
hold on;
plot([-medianRTs -medianRTs],[min(yaxislabel) max(yaxislabel)],'k-')
axis square
set(gca,'clim',[-0.05 0.05],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Succesfful-Failed')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')


subplot(234)
contourf(xaxislabel,yaxislabel,zmap,40,'linecolor','none')
axis square
set(gca,'clim',[-3 3],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Unthresholded Z map')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')

% apply uncorrected threshold
subplot(235)
contourf(xaxislabel,yaxislabel,zmap,40,'linecolor','none')
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false;
zmapthresh=logical(zmapthresh);
hold on
contour(xaxislabel,yaxislabel,zmapthresh,1,'linecolor','k')

axis square
set(gca,'clim',[-3 3],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Unthresholded Z map')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')


% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

subplot(236)
contourf(xaxislabel,yaxislabel,zmapthresh,40,'linecolor','none')
hold on;
plot([-medianRTs -medianRTs],[min(yaxislabel) max(yaxislabel)],'k-')
axis square
set(gca,'clim',[-3 3],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Cluster-corrected Z map')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')

signYtime=yaxislabel(mean(zmapthresh,2)>0);

hold on;
plot([min(xaxislabel) max(xaxislabel)],[min(signYtime) min(signYtime)],'r-')
plot([min(xaxislabel) max(xaxislabel)],[max(signYtime) max(signYtime)],'r-')


