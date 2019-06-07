clear all

load mi_pat_ids.mat

npats = numel(pat_ids);

pars = set_params9();
pars.bp_flag = 0;

froot = pars.froot;

freqs = pars.freqs;
nFreqBands = pars.nFreqBands;
freqBands = pars.freqBands;

lags = [-250:10:-100 -95:5:-25 -20:20 25:5:95 100:10:250];
lag_inds = lags>=-20 & lags<=20;
mid_lags = -20:20;
left_inds = 1:6;
%mid_inds = 7:47;
right_inds = 98:103;
pars.left_inds = left_inds;
pars.right_inds = right_inds;
pars.pval_thresh = .00001;
tvals = -.5:.1:4;
type_strs = {'enc','ret'};

plots = plot3brains_base;

for i = 1:3
    axes(plots(i))
    alpha(.3)
end

for n = 1:npats
    tic;
    pat_s = pat_ids{n};
    fprintf('%s\n',pat_s)
    
    pat_dir_out = [froot pat_s];
    load([pat_dir_out '/' pat_s '_info.mat'])
    
    %set path to subject directory for output
    pat_dir_out = [froot pat_s];

    %Create powDir
    powDir = [pat_dir_out '/powDir'];
    pars.powDir = powDir;
    
    num_events = sum(num_each);
    
    load([pat_dir_out '/bad_chans.mat'])
    els = load_electrode_info(pat_s,0);
    nelecs = numel(electrodes);
    npairs = nelecs*(nelecs-1)/2;
    pair_to_inds = cell(npairs,1);
    is_good_pair = true(npairs,1);
    pair_ind = 1;
    for e1 = 1:nelecs
        for e2 = e1+1:nelecs
            pair_to_inds{pair_ind} = [e1 e2];
            if ismember(e1,bad_chans) || ismember(e2,bad_chans)
                is_good_pair(pair_ind)=false;
            end
            pair_ind=pair_ind+1;
        end
    end
    parfor p = 1:npairs
    %for p = 107
        if is_good_pair(p)
            mi_par_fun_by_pair(pat_s,p,pars,pair_to_inds,lags);
        end
    end
%    keyboard
    mi_at_peaks_sig_pairs = [];%will be pairs X trials X type X time
    peak_vals_sig_pairs = [];%will be pairs X 1
    mi_vs_tau_sig_pairs = [];%will be pairs X trials X type X time X lag
    pair_inds = [];%will be pairs X 1
    pair_count = 1;
    
    mkdir(['~/Desktop/figures/' pat_s '/corr_task_effects_4'])
    
    for p = 1:npairs
        if is_good_pair(p)
            load([powDir '/pair' num2str(p) '_mi_data.mat'],'pvals','peak_val')
            if any(pvals(:)<pars.pval_thresh) && peak_val~=0
                load([powDir '/pair' num2str(p) '_mi_data.mat'],'mi_at_peaks_norm','normalized_mi_vals')
                %mi_at_peaks_norm: trials X type X time,
                %normalized_mi_vals: trials X type X time X lag
                mi_at_peaks_sig_pairs(pair_count,:,:,:) = mi_at_peaks_norm;
                mi_vs_tau_sig_pairs(pair_count,:,:,:,:) = normalized_mi_vals;
                peak_vals_sig_pairs(pair_count) = peak_val;
                pair_inds(pair_count) = p;
                pair_count = pair_count + 1;
                
                %Make MI vs Time figure
                figure(2)
                set(gcf,'Position',[496 767 1048 489])
                clf
                ylims = [0 0];
                for type = 1:2
                    subplot(1,2,type)
                    hold on
                    mean_vals = squeeze(nanmean(squeeze(mi_at_peaks_sig_pairs(pair_count-1,:,type,:))));
                    std_err_vals = squeeze(std_err(squeeze(mi_at_peaks_sig_pairs(pair_count-1,:,type,:))));
                    shadedErrorBar(tvals,mean_vals,std_err_vals,'b')

                    thisylim = get(gca,'ylim');
                    ylims(1) = min(ylims(1),thisylim(1));
                    ylims(2) = max(ylims(2),thisylim(2));
                end
                for type = 1:2
                    subplot(1,2,type)
                    set(gca,'ylim',ylims)
                    title(type_strs{type})
                    xlabel('Time (s)')
                    ylabel('z-scored |correlation| normalized within session')
                end
                elnums = pair_to_inds{p};
                el1 = els(elnums(1));
                el2 = els(elnums(2));

                suptitle([num2str(elnums(1)) ': ' el1.Loc1(1) ' ' el1.Loc3   '  ' num2str(elnums(2)) ': ' el2.Loc1(1) ' ' el2.Loc3 ' Tau: ' num2str(peak_val)]);
                el1toplot = el1;
                el2toplot = el2;
                print(['~/Desktop/figures/' pat_s '/corr_task_effects_4/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) '.eps'],'-depsc','-loose');
                pause(.1)
                
                %Make MI vs Tau movies
                for type = 1:2
                    vidObj = VideoWriter(['~/Desktop/figures/' pat_s '/corr_task_effects_4/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) type_strs{type} '.avi']);
                    set(vidObj,'FrameRate',4);
                    open(vidObj);
                    
                    all_mean_vals = squeeze(nanmean(squeeze(mi_vs_tau_sig_pairs(pair_count-1,:,:,:,:))));%type X time X lag
                    minval = min(all_mean_vals(:));
                    maxval = max(all_mean_vals(:));
                    ylims = [minval-1 maxval+1];
                    
                    for t = 1:46
                        figure(3)
                        clf
                        hold on
                        mean_vals = squeeze(nanmean(squeeze(mi_vs_tau_sig_pairs(pair_count-1,:,type,t,:))));
                        std_err_vals = squeeze(std_err(squeeze(mi_vs_tau_sig_pairs(pair_count-1,:,type,t,:))));

                        shadedErrorBar(lags,mean_vals,std_err_vals,'b');

                        ylim(ylims);
                        title([num2str(tvals(t)) ' sec'])
                        xlabel('Tau (ms)')
                        ylabel('Correlation z-scored to baseline lags')
                        M = getframe(gcf); 
                        writeVideo(vidObj,M);
                    end
                    close(vidObj);
                end
                
                %Make brain plots
                figure(1)
                try
                    remove_electrode(h1);
                    remove_electrode(h2);
                catch

                end
                h1 = plot_electrode_on_brains(plots,el1toplot,3,'r');
                h2 = plot_electrode_on_brains(plots,el2toplot,3,'b');
                print(['~/Desktop/figures/' pat_s '/corr_task_effects_4/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) 'brain.eps'],'-depsc','-loose');
            end
        end
    end
    save([pat_dir_out '/' pat_s '_mi_vals_sig_pairs.mat'],'mi_at_peaks_sig_pairs','peak_vals_sig_pairs','-v7.3')
    %Add computing reinstatement and making plots (brain, networks,
    %time series, mi vs tau movies, and reinstatement
    %{
    for p = 1:size(mi_at_peaks_sig_pairs,1)
        figure(2)
        set(gcf,'Position',[496 767 1048 489])
        clf
        ylims = [0 0];
        for type = 1:2
            subplot(1,2,type)
            hold on
            mean_vals = squeeze(nanmean(squeeze(mi_at_peaks_sig_pairs(p,:,type,:))));
            std_err_vals = squeeze(std_err(squeeze(mi_at_peaks_sig_pairs(p,:,type,:))));
            shadedErrorBar(tvals,mean_vals,std_err_vals,'b')
            
            thisylim = get(gca,'ylim');
            ylims(1) = min(ylims(1),thisylim(1));
            ylims(2) = max(ylims(2),thisylim(2));
        end
        for type = 1:2
            subplot(1,2,type)
            set(gca,'ylim',ylims)
        end
        %{
        ylims = [0 0];
        for type = 1:2
            subplot(3,2,2+type)
            hold on
            
            [~,t_ind] = max(squeeze(nanmean(squeeze(mi_at_peaks_sig_pairs(p,:,type,:)))));
            
            mean_vals = squeeze(nanmean(squeeze(mi_vs_tau_sig_pairs(p,:,type,t_ind,:))));
            std_err_vals = squeeze(std_err(squeeze(mi_vs_tau_sig_pairs(p,:,type,t_ind,:))));
     
            shadedErrorBar(lags,mean_vals,std_err_vals,'b');

            
            thisylim = get(gca,'ylim');
            ylims(1) = min(ylims(1),thisylim(1));
            ylims(2) = max(ylims(2),thisylim(2));
        end
        for type = 1:2
            subplot(3,2,2+type)
            set(gca,'ylim',ylims)
        end
        
        for type = 1:2
            subplot(3,2,4+type)
            hold on
            
            [~,t_ind] = min(squeeze(nanmean(squeeze(mi_at_peaks_sig_pairs(p,:,type,:)))));
            
            mean_vals = squeeze(nanmean(squeeze(mi_vs_tau_sig_pairs(p,:,type,t_ind,:))));
            std_err_vals = squeeze(std_err(squeeze(mi_vs_tau_sig_pairs(p,:,type,t_ind,:))));
     
            shadedErrorBar(lags,mean_vals,std_err_vals,'b');
            
        end
        for type = 1:2
            subplot(3,2,4+type)
            set(gca,'ylim',ylims)
        end
        %}
        elnums = pair_to_inds{pair_inds(p)};
        el1 = els(elnums(1));
        el2 = els(elnums(2));
        
        suptitle([el1.Loc1(1) ' ' el1.Loc3 '  ' el2.Loc1(1) ' ' el2.Loc3]);
        el1toplot = el1;
        el2toplot = el2;
        print(['~/Desktop/figures/' pat_s '/corr_task_effects_4/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) '.eps'],'-depsc','-loose');
        
        pause(.1)
        %movefile(['~/Desktop/figures/' pat_s '/connected/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) '.eps'],['~/Desktop/figures/' pat_s '/well_connected/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) '.eps']);
        %movefile(['~/Desktop/figures/' pat_s '/connected/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) '_brain.eps'],['~/Desktop/figures/' pat_s '/well_connected/' pat_s '_' num2str(elnums(1)) '_' num2str(elnums(2)) '_brain.eps']);
        
        
    end
    %}
    
    toc;
end
