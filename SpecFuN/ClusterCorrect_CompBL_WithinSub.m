% cluster-based correction 
% within-subject
function zmapthresh = ClusterCorrect_CompBL_WithinSub(X,Y,xaxislabel,yaxislabel,medianRTs)

n_permutes = 1000;
num_1=size(X,2);
num_2=size(X,3);

voxel_pval   = 0.05;
cluster_pval = 0.10;

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_1);
permuted_vals    = zeros(n_permutes,num_1,num_2);
max_clust_info   = zeros(n_permutes,1);
nsub=size(X,1);


baseidx(1)=dsearchn(yaxislabel,-1000);
baseidx(2)=dsearchn(yaxislabel,-500);

baseidx(3)=dsearchn(xaxislabel,-4000);
baseidx(4)=dsearchn(xaxislabel,-3500);

X2=X-mean(Y);
% compute actual t-test of difference
realbaselines = X2(:,baseidx(1):baseidx(2),baseidx(3):baseidx(4));
realmean      = squeeze(mean(X2-mean(mean(realbaselines,2),3)));
nTimepoints=numel(yaxislabel);


for permi=1:n_permutes
    cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1); 
    
    permuted_vals(permi,:,:) =     squeeze(mean(X2(:,[cutpoint:end 1:cutpoint-1],[cutpoint:end 1:cutpoint-1])-mean(mean(realbaselines,2),3)));
    % btw, using bsxfun instead of repmat increases the speed 
    % of this loop by a factor of ~5.
end

  
zmap = (realmean-squeeze(mean(permuted_vals))) ./ squeeze(std(permuted_vals));
threshmean = realmean;
threshmean(abs(zmap)<norminv(1-voxel_pval))=0;

% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    fakecorrsz = squeeze((permuted_vals(permi,:,:)-mean(permuted_vals,1)) ./ std(permuted_vals,[],1) );
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(fakecorrsz);
    max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps
    % using cellfun here eliminates the need for a slower loop over cells
end

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
plot([-0 -0],[min(yaxislabel) max(yaxislabel)],'k-')
axis square
set(gca,'clim',[-3 3],'xlim',[min(xaxislabel) max(xaxislabel)])
title('Cluster-corrected Z map')
xlabel('Retrieval (ms)'), ylabel('Encoding (ms)')



