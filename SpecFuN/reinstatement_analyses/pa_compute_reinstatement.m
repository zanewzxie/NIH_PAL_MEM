%This script computes the power for the paired associates task

close all
clear all

%Comment this out to run analysis on all patients
%Otherwise, analysis is run just on patients listed in new_pat_ids
%new_pat_ids = {'NIH016','NIH017','NIH018','TJ060','TJ064','TJ065'};

%new_pat_ids = {'NIH022','NIH023','NIH024','NIH025','NIH026','NIH028'};

%new_pat_ids = {'NIH032'};

if exist('new_pat_ids','var')
    add_pats_to_all_pat_ids(new_pat_ids);
    add_pats_to_all_pat_ids_behav(new_pat_ids);
    set_pat_ids(new_pat_ids);
    %createDirectories;
else
    set_pat_ids;
end

load pat_ids.mat

%pat_ids = {'NIH034'};
%pat_ids = {'NIH016','NIH018','NIH023','NIH025','NIH026','NIH029','NIH030','NIH032','NIH034'};
npats = numel(pat_ids);

pars = set_params9();

froot = pars.froot;

for n = 1:npats
%for n = 1:6
    tic;
    pat_s = pat_ids{n};
    display(pat_s)
    %set path to subject directory
    pat_dir_out = [froot pat_s];
    
    load([pat_dir_out '/' pat_s '_info.mat'],'electrodes','bp_electrodes','eventsByCondition','all_events','num_each')
    
    %set powDir
    powDir = [pat_dir_out '/powDir'];    
    load([powDir '/avgPowData500monopolar.mat'],'avgPowData','pars')
    %load([powDir '/avgPowData500.mat'],'avgPowData','pars')
    
    corr_ind = 1:num_each(1);
    incorr_ind = num_each(1)+1:sum(num_each);
    nmatch = min(numel(corr_ind),numel(incorr_ind));
    
    output = cell(2,1);
    output_match = cell(2,1);
    output_match_cue = cell(2,1);
    output_theta = cell(2,1);
    output_hg = cell(2,1);
    output_theta_cue = cell(2,1);
    output_hg_cue = cell(2,1);
    r_each = cell(2,1);
    r_each_theta = cell(2,1);
    r_each_hg = cell(2,1);
    r_each_theta_cue = cell(2,1);
    r_each_hg_cue = cell(2,1);
    r_each_match = cell(2,1);
    r_each_match_cue = cell(2,1);
    for zz = 1:2
        if zz==1
            inds = corr_ind;
        else
            inds = incorr_ind;
        end
        inds_match = randsample(inds,nmatch);
        encMatCorr = squeeze(avgPowData(inds,1,:,:,:));
        recMatCorr = squeeze(avgPowData(inds,3,:,:,:));
        %encMatCorr is trial X electrode X freq X time

        encMatCorr = permute(encMatCorr,[1 3 4 2]);%trial X freq X time X electrode
        recMatCorr = permute(recMatCorr,[1 3 4 2]);
        
        encMatCorrTheta = select_features_by_freq(encMatCorr,1);
        recMatCorrTheta = select_features_by_freq(recMatCorr,1);
        
        [output_theta{zz}, r_each_theta{zz}] = compute_reinstatement(encMatCorrTheta,recMatCorrTheta);
        
        encMatCorrHG = select_features_by_freq(encMatCorr,5);
        recMatCorrHG = select_features_by_freq(recMatCorr,5);
        
        [output_hg{zz}, r_each_hg{zz}] = compute_reinstatement(encMatCorrHG,recMatCorrHG);
        
        encMatCorr = select_all_features(encMatCorr);
        recMatCorr = select_all_features(recMatCorr);

        [output{zz}, r_each{zz}] = compute_reinstatement(encMatCorr,recMatCorr);
        
        encMatCorrMatch = squeeze(avgPowData(inds_match,1,:,:,:));
        recMatCorrMatch = squeeze(avgPowData(inds_match,3,:,:,:));
        %encMatCorr is trial X electrode X freq X time

        encMatCorrMatch = permute(encMatCorrMatch,[1 3 4 2]);
        recMatCorrMatch = permute(recMatCorrMatch,[1 3 4 2]);

        encMatCorrMatch = select_all_features(encMatCorrMatch);
        recMatCorrMatch = select_all_features(recMatCorrMatch);

        [output_match{zz}, r_each_match{zz}] = compute_reinstatement(encMatCorrMatch,recMatCorrMatch);
        
        encMatCorrMatch = squeeze(avgPowData(inds_match,1,:,:,:));
        recMatCorrMatch = squeeze(avgPowData(inds_match,2,:,:,:));
        %encMatCorr is trial X electrode X freq X time

        encMatCorrMatch = permute(encMatCorrMatch,[1 3 4 2]);
        recMatCorrMatch = permute(recMatCorrMatch,[1 3 4 2]);

        encMatCorrMatch = select_all_features(encMatCorrMatch);
        recMatCorrMatch = select_all_features(recMatCorrMatch);

        [output_match_cue{zz}, r_each_match_cue{zz}] = compute_reinstatement(encMatCorrMatch,recMatCorrMatch);
        
        encMatCorrMatch = squeeze(avgPowData(inds_match,1,:,:,:));
        recMatCorrMatch = squeeze(avgPowData(inds_match,2,:,:,:));
        %encMatCorr is trial X electrode X freq X time

        encMatCorrMatch = permute(encMatCorrMatch,[1 3 4 2]);
        recMatCorrMatch = permute(recMatCorrMatch,[1 3 4 2]);

        encMatCorrMatch = select_features_by_freq(encMatCorrMatch,5);
        recMatCorrMatch = select_features_by_freq(recMatCorrMatch,5);

        [output_hg_cue{zz}, r_each_hg_cue{zz}] = compute_reinstatement(encMatCorrMatch,recMatCorrMatch);
        
        encMatCorrMatch = squeeze(avgPowData(inds_match,1,:,:,:));
        recMatCorrMatch = squeeze(avgPowData(inds_match,2,:,:,:));
        %encMatCorr is trial X electrode X freq X time

        encMatCorrMatch = permute(encMatCorrMatch,[1 3 4 2]);
        recMatCorrMatch = permute(recMatCorrMatch,[1 3 4 2]);

        encMatCorrMatch = select_features_by_freq(encMatCorrMatch,1);
        recMatCorrMatch = select_features_by_freq(recMatCorrMatch,1);

        [output_theta_cue{zz}, r_each_theta_cue{zz}] = compute_reinstatement(encMatCorrMatch,recMatCorrMatch);
        
        
    end
    save([pat_dir_out '/' pat_s '_reinstatement.mat'],'output*','r_each*')
end
