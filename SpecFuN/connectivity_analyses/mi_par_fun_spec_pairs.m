function mi_par_fun_spec_pairs(eeg_fname,pairs_to_use,lags,pair_to_inds)
    load(eeg_fname,'avg_ref_eeg')
    nlags = numel(lags);
    timesteps = 46;
    points_per_win = 500;
    points_per_slide = 100;
    mi_vals_by_lag_spec_pairs = nan(2,numel(pairs_to_use),timesteps,nlags);
    
    maxlag = max(lags);
    alllags = min(lags):max(lags);
    
    
    for type = 1:2
        pair_count = 1;
        for pair_ind = pairs_to_use
            els = pair_to_inds{pair_ind};
            e1 = els(1);
            e2 = els(2);
            mi_vals = nan(timesteps,nlags);
            %display(pair_ind)
            for i = 1:timesteps
                ind1 = 6001 + (i-1)*points_per_slide;%2000ms buffer, 5000msbefore, 10000ms after
                data1 = squeeze(avg_ref_eeg(type,e1,ind1-maxlag:ind1+points_per_win+maxlag-1));
                data1buf = buffer(data1,points_per_win,points_per_win-1,'nodelay');
                data2 = squeeze(avg_ref_eeg(type,e2,ind1:ind1+points_per_win-1));

                for lag_ind = 1:nlags
                    lag_val = lags(lag_ind);
                    lag = alllags==lag_val;
                    data1 = zscore(data1buf(:,lag));
                    data2 = zscore(data2);
                    mi_vals(i,lag_ind) = simple_mi(data1,data2);
                end
            end
            mi_vals_by_lag_spec_pairs(type,pair_count,:,:) = mi_vals;
            pair_count = pair_count + 1;
        end
    end
    save(eeg_fname,'mi_vals_by_lag_spec_pairs','-append')
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