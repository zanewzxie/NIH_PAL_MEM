% cluster-based correction 
% within-subject
function zmapthresh = ClusterCorrect_WithinSub(X,Y,xaxislabel,yaxislabel,medianRTs)

n_permutes = 1000;
num_1=size(X,2);
num_2=size(X,3);

voxel_pval   = 0.005;
cluster_pval = 0.01;

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_1);
permuted_vals    = zeros(n_permutes,num_1,num_2);
max_clust_info   = zeros(n_permutes,1);
nsub=size(X,1);


baseidx(1)=dsearchn(yaxislabel,-1000);
baseidx(2)=dsearchn(yaxislabel,-500);

baseidx(3)=dsearchn(xaxislabel,-4000);
baseidx(4)=dsearchn(xaxislabel,-3500);

% compute actual t-test of difference
realbaselines = X(:,baseidx(1):baseidx(2),baseidx(3):baseidx(4));


[h,p,ci,stat]=ttest(X-mean(mean(realbaselines,2),3),0);
real_t=squeeze(stat.tstat);
realmean=mean(X-mean(mean(realbaselines,2),3));

eegpower=[X;Y];
triallength=[size(X,1) size(Y,1)];
ntrials=size(eegpower,1);

% nTimepoints=num_1;
% 
% cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1);
nTimepoints=numel(yaxislabel);

    
% generate pixel-specific null hypothesis parameter distributions
clear permuted_tvals max_pixel_pvals
for permi = 1:n_permutes
    
%     fake_condition_mapping = randsample(1:ntrials,ntrials);
    
    cutpoint = randsample(2:nTimepoints-diff(baseidx(1:2))-2,1); 
    
    % compute t-map of null hypothesis
    [h,p,ci,stat]=ttest(X(:,[cutpoint:end 1:cutpoint-1],[cutpoint:end 1:cutpoint-1])-mean(mean(realbaselines,2),3),0);
    tmap   = squeeze(stat.tstat);
    
    % save all permuted values
    permuted_tvals(permi,:,:) = tmap;
    
    % save maximum pixel values
    max_pixel_pvals(permi,:) = [ min(tmap(:)) max(tmap(:)) ];
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,ntrials-2))=0;
    
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
zmapthresh(zmapthresh<0)=0;

subplot(236)
contourf(xaxislabel,yaxislabel,zmapthresh,40,'linecolor','none')
hold on;
plot([-medianRTs -medianRTs],[min(yaxislabel) max(yaxislabel)],'k-')
plot([-0 -0],[min(yaxislabel) max(yaxislabel)],'k-')
axis square
set(gca,'clim',[-3 3],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Cluster-corrected Z map')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')



