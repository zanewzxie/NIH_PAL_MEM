function mi_par_fun_by_pair(pat_s,p,pars,pair_to_inds,lags)
    pat_dir_out = [pars.froot pat_s];
    load([pat_dir_out '/' pat_s '_info.mat'])
    num_events = sum(num_each);
    powDir = [pat_dir_out '/powDir'];
    
    elnums = pair_to_inds{p};
    e1 = elnums(1);
    e2 = elnums(2);
    nlags = numel(lags);
    timesteps = 10;
    points_per_win = 500;
    points_per_slide = 500;
    
    maxlag = max(lags);
    alllags = min(lags):max(lags);
    left_inds = pars.left_inds;
    right_inds = pars.right_inds;
    
    mi_vals_by_event = nan(100,2,timesteps,nlags);
    
    lag_inds = ismember(alllags,lags);
    
    for e = 1:100
        eeg_fname = [powDir '/event' num2str(e) 'eegData.mat'];
        load(eeg_fname,'avg_ref_eeg')
        %display(e)
        for type = 1:2
            data1buf = nan(points_per_win,nlags,timesteps);
            data2buf = nan(points_per_win,nlags,timesteps);
            
            for i = 1:timesteps
                ind1 = 6001 + (i-1)*points_per_slide;%2000ms buffer, 5000msbefore, 10000ms after
                data1 = squeeze(avg_ref_eeg(type,e1,ind1-maxlag:ind1+points_per_win+maxlag-1));
                data1 = buffer(data1,points_per_win,points_per_win-1,'nodelay');
                data2 = squeeze(avg_ref_eeg(type,e2,ind1:ind1+points_per_win-1));
                data1buf(:,:,i) = data1(:,lag_inds);
                data2buf(:,:,i) = repmat(data2,1,size(data1buf,2));
            end
            mi_vals_by_event(e,type,:,:) = abs(squeeze(custom_corr_3d(data1buf,data2buf)))';
        end
    end
    
    %z-score by baseline lags
    
    baseline_vals = mi_vals_by_event(:,:,:,[left_inds right_inds]);%trials X type X time X lag
    baseline_vals = single(permute(baseline_vals,[4 1 2 3]));%lag X trials X type X time
    mean_vals_to_use = squeeze(nanmean(baseline_vals));%trials X type X time
    std_vals_to_use = squeeze(nanstd(baseline_vals));%trials X type X time
    
    normalized_mi_vals = single(bsxfun(@minus,mi_vals_by_event,mean_vals_to_use));
    normalized_mi_vals = single(bsxfun(@rdivide,normalized_mi_vals,std_vals_to_use));%trials X type X time X lag
    
    %Find peak tau
    
    mean_norm_vals = squeeze(nanmean(normalized_mi_vals));%type X time X lag
    mean_norm_vals = squeeze(nanmean(mean_norm_vals));%time X lag
    mean_norm_vals = squeeze(nanmean(mean_norm_vals));%1 X lag
    
    
    
    [~,peak_ind] = max(mean_norm_vals);%scalar value
    
    %Get MI vs time at peak tau
    mi_at_peaks = normalized_mi_vals(:,:,:,peak_ind);%trials X type X time
    
    %Normalize within session
    
    sessionInds = getStructField(all_events.STUDY_PAIR,'session');
    sessionNums = unique(sessionInds);
    
    mi_at_peaks_norm = nan(size(mi_at_peaks));%trials X type X time

    for s = sessionNums
        inds = sessionInds==s;
        inds = inds(1:100);
        sessionData = mi_at_peaks(inds,:,:);
        mean_data = nanmean(sessionData(:));
        std_data = nanstd(sessionData(:));
        mi_at_peaks_norm(inds,:,:) = (mi_at_peaks(inds,:,:) - mean_data)./std_data;
    end
    
    % do all t-tests to determine if dynamic
    mi_at_peaks_to_comp1 = nan(size(mi_at_peaks,1),size(mi_at_peaks,2),45);
    mi_at_peaks_to_comp2 = nan(size(mi_at_peaks,1),size(mi_at_peaks,2),45);
    count = 1;
    for t1 = 1:10
        for t2 = t1+1:10
            mi_at_peaks_to_comp1(:,:,count) = mi_at_peaks_norm(:,:,t1);
            mi_at_peaks_to_comp2(:,:,count) = mi_at_peaks_norm(:,:,t2);
            count = count + 1;
        end
    end
    %second_time_bin = repmat(mi_at_peaks_norm(:,:,2),[1 1 10]);
    [~,pvals] = ttest(mi_at_peaks_to_comp1,mi_at_peaks_to_comp2);
    pvals = squeeze(pvals);%type X comparisons
    peak_val = lags(peak_ind);
    
    
    
    save([powDir '/pair' num2str(p) '_mi_data.mat'],'pvals','peak_val')
    if any(pvals(:)<pars.pval_thresh) && peak_val~=0
        %Compute and save MI vs time at peak tau for this pair
        timesteps = 46;
        points_per_win = 500;
        points_per_slide = 100;
        
        mi_vals_by_event = nan(num_events,2,timesteps,nlags);

        for e = 1:num_events
            eeg_fname = [powDir '/event' num2str(e) 'eegData.mat'];
            load(eeg_fname,'avg_ref_eeg')
            for type = 1:2
                data1buf = nan(points_per_win,nlags,timesteps);
                data2buf = nan(points_per_win,nlags,timesteps);

                for i = 1:timesteps
                    ind1 = 6001 + (i-1)*points_per_slide;%2000ms buffer, 5000msbefore, 10000ms after
                    data1 = squeeze(avg_ref_eeg(type,e1,ind1-maxlag:ind1+points_per_win+maxlag-1));
                    data1 = buffer(data1,points_per_win,points_per_win-1,'nodelay');
                    data2 = squeeze(avg_ref_eeg(type,e2,ind1:ind1+points_per_win-1));
                    data1buf(:,:,i) = data1(:,lag_inds);
                    data2buf(:,:,i) = repmat(data2,1,size(data1buf,2));
                end
                mi_vals_by_event(e,type,:,:) = abs(squeeze(custom_corr_3d(data1buf,data2buf)))';
            end
        end

        %z-score by baseline lags

        baseline_vals = mi_vals_by_event(:,:,:,[left_inds right_inds]);%trials X type X time X lag
        baseline_vals = single(permute(baseline_vals,[4 1 2 3]));%lag X trials X type X time
        mean_vals_to_use = squeeze(nanmean(baseline_vals));%trials X type X time
        std_vals_to_use = squeeze(nanstd(baseline_vals));%trials X type X time

        normalized_mi_vals = single(bsxfun(@minus,mi_vals_by_event,mean_vals_to_use));
        normalized_mi_vals = single(bsxfun(@rdivide,normalized_mi_vals,std_vals_to_use));%trials X type X time X lag
        
        %Get MI vs time at peak tau
        mi_at_peaks = normalized_mi_vals(:,:,:,peak_ind);%trials X type X time

        %Normalize within session

        sessionInds = getStructField(all_events.STUDY_PAIR,'session');
        sessionNums = unique(sessionInds);

        mi_at_peaks_norm = nan(size(mi_at_peaks));%trials X type X time

        for s = sessionNums
            inds = sessionInds==s;
            sessionData = mi_at_peaks(inds,:,:);
            mean_data = nanmean(sessionData(:));
            std_data = nanstd(sessionData(:));
            mi_at_peaks_norm(inds,:,:) = (mi_at_peaks(inds,:,:) - mean_data)./std_data;
        end
        save([powDir '/pair' num2str(p) '_mi_data.mat'],'mi_at_peaks_norm','normalized_mi_vals','-append')
    end
end

function [MI,Hx,Hy,Hxy] = simple_mi(X,Y)
    nbins=30;
    %X = data1;
    %Y = data2;
    %bins={[min([X;Y])+(max([X;Y])-min([X;Y]))/nbins/2:(max([X;Y])-min([X;Y]))/nbins:max([X;Y])-(max([X;Y])-min([X;Y]))/nbins/2],[min([X;Y])+(max([X;Y])-min([X;Y]))/nbins/2:(max([X;Y])-min([X;Y]))/nbins:max([X;Y])-(max([X;Y])-min([X;Y]))/nbins/2]};
    
    bins = linspace(-4.8333,4.8333,nbins);
    bins = {bins;bins};
    [joint]=hist3([X,Y],bins);

    joint=joint/sum(joint(:));
    margX=sum(joint,2);
    margY=sum(joint,1);

    joint(joint==0)=eps;
    margX(margX==0)=eps;
    margY(margY==0)=eps;
    
    Hx = -sum(margX.*log2(margX));
    Hy = -sum(margY.*log2(margY));
    Hxy = -sum(sum(joint.*log2(joint)));
    MI = Hx + Hy - Hxy;
    %MI=sum(sum(joint.* bsxfun(@minus,bsxfun(@minus,log2(joint),log2(margX)),log2(margY))));%-(Mxy-Mx-My-1)/(2*length(X));

end