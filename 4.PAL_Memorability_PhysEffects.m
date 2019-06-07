clear all

cd('/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/NIH_PAL_MEM')


%add toolboxes
addpath(genpath('/Volumes/Zane/Matlab/eeg_toolbox/trunk'));
addpath(genpath('/Volumes/Zane/Matlab/dungeon_toolbox_17a'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/EEG_Preprocessing'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/IEM_ester'));
addpath(genpath('/Volumes/Zane/NIH_NINDS/Data_InProgress/SpecFuN'));

curpath='/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/NIH_PAL_MEM'
addpath([curpath '/SpecFuN']);



%% extract PAL behaviroal data
allsub=[26:64 66];
rootEEGdir = '/Volumes/Zane/NIH_FRNU_ROOT';                      %office-local

load('SubjTable_palRAMword.mat')
load('PAL_Memo_PALRAM.mat')
load('WordtoAnalyze_ATL_2ROIs.mat');

ToAnalyzeSub=allsub;
medianMemor = median(PAL_Memo.Responsememorability);
uniquewordsID_10=PAL_Memo.uniquewordsID_10;

load('WordtoAnalyzeTReduced.mat');

%%
%%- FILTERING OPTIONS
WHICH_DECOMP    = 'wavelet';  %'wavelet' 'multitaper' 'multitaper200'  % MT assumes a 500ms window, use mt200 for 200ms window; 'hilbert' not supported yet
WHICH_BANDS     = 'standard';  %WHICH_BANDS = 'standard' 'NoiseAnalysis' 'LowHigh' 'HighMidLow' 'HighMidLow2' 'JustHigh'
%- select bands for decomposition
[waveletFreqs waveletWidth waveletFreqLabels freqBandAr freqBandYticks freqBandYtickLabels hilbertFreqs multitaperObj decompStruct] = prepFreqList(WHICH_BANDS,WHICH_DECOMP);

freqBandYticksSplit=[exp(linspace(log(4),log(16),10));exp(linspace(log(70),log(150),10))]; % only bin the data into 3 bins.

waveletFreqs3 = exp(linspace(log(4),log(150),30));

freqBandYticks = [4   8   16   32   70   150]; % overwrite the old ticks.
WordtoAnalyzeTReduced
mkdir('data_ATL_3Event_100')
curpath=pwd;
outpath=fullfile(curpath,'data_ATL_3Event_100');


% loop through subjects.
%% Preprocessing Power spectrum 

for isubtem=1:length(WordtoAnalyzeTReduced.subID)
    
    
    clear SubjResults
    isub = WordtoAnalyzeTReduced.subID(isubtem);
    % duration variable across subjects
    % nanmean(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>50)) + 500 + 1000
    duration = 10000;
    offset   = -5000;
    buffer   = 500; % the butter is built-in in the duration
    resamp   = 250;
    time  =  downsample(linspace(offset,duration+offset,duration),1000/resamp);
    
    chan=WordtoAnalyzeTReduced.chan{isubtem};
    
    % only analyze the middle and inferior temporal gryus. The superior sites may be influenced by sounds
    electrodeToAna=ismember(WordtoAnalyzeTReduced.chan{isubtem},WordtoAnalyzeTReduced.ATLchan{isubtem});
    %     electrodeToAna=ismember(WordtoAnalyzeTReduced.chan{isubtem},WordtoAnalyzeTReduced.MTLchan{isubtem});
    chan=chan(electrodeToAna);
    %%
    if length(chan)>=2 % Ony analyze the data from available channel.
        
        temChans=[];
        allcleanTrials=[];
        
        Study_Event=WordtoAnalyzeTReduced.SubjTable{isubtem}.alignedEvents.STUDY_PAIR(WordtoAnalyzeTReduced.InDX{isubtem});
        TEST_PROBE_Event=WordtoAnalyzeTReduced.SubjTable{isubtem}.alignedEvents.TEST_PROBE(WordtoAnalyzeTReduced.InDX{isubtem});
        Recall_Event=WordtoAnalyzeTReduced.SubjTable{isubtem}.alignedEvents.RecallEvent(WordtoAnalyzeTReduced.InDX{isubtem});
        
        allssession=[TEST_PROBE_Event.session];
        uniquesessions=unique(allssession);
        
        for iss=1:length(uniquesessions)
            iss
            sessionGA_Study=gete_ms('global_avg_good',Study_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[200 ],'low',2,resamp);
            sessionGA_Probe=gete_ms('global_avg_good',TEST_PROBE_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[200 ],'low',2,resamp);
            sessionGA_Recall=gete_ms('global_avg_good',Recall_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[200 ],'low',2,resamp);
            
            clear ProbeRawEEG ProbePowerlog WavePow_log RecallRawEEG StudyRawEEG
            parfor ic=1:length(chan)
                % get EEG
                ic
                StudyRawEEG(ic,:,:)=gete_ms(chan{ic},Study_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[200 ],'low',2,resamp);
                StudyRawEEG(ic,:,:)=squeeze(StudyRawEEG(ic,:,:))-sessionGA_Study; % reference to global average of all electrode within that session.'
                
                
                ProbeRawEEG(ic,:,:)=gete_ms(chan{ic},TEST_PROBE_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[200 ],'low',2,resamp);
                ProbeRawEEG(ic,:,:)=squeeze(ProbeRawEEG(ic,:,:))-sessionGA_Probe; % reference to global average of all electrode within that session.'
                
                RecallRawEEG(ic,:,:)=gete_ms(chan{ic},Recall_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[200 ],'low',2,resamp);
                RecallRawEEG(ic,:,:)=squeeze(RecallRawEEG(ic,:,:))-sessionGA_Recall; % reference to global average of all electrode within that session.'
            end
            
            StudyRawEEGAll{iss} = StudyRawEEG;
            ProbeRawEEGAll{iss} = ProbeRawEEG;
            RecallRawEEGAll{iss} = RecallRawEEG;
            
            clnWeights=[2.3 2.3]; % cleaning weight
            clfig=figure(2000);clf
            FIG_TITLE= ['Subj' num2str(isub) '--Session' num2str(iss) 'Study_cleaning'];
            [iChanClean1,iEvClean1,strClean] = jwCleanEEGevents_v01(StudyRawEEG,clfig,FIG_TITLE,clnWeights);
            disp(strClean);
            print(clfig,fullfile(outpath,[FIG_TITLE '.pdf']),'-dpdf','-bestfit');
            
            clfig=figure(2001);clf
            FIG_TITLE= ['Subj' num2str(isub) '--Session' num2str(iss) 'Probe_cleaning'];
            [iChanClean2,iEvClean2,strClean] = jwCleanEEGevents_v01(ProbeRawEEG,clfig,FIG_TITLE,clnWeights);
            disp(strClean);
            print(clfig,fullfile(outpath,[FIG_TITLE '.pdf']),'-dpdf','-bestfit');
            
            clfig=figure(2002);clf
            FIG_TITLE= ['Subj' num2str(isub) '--Session' num2str(iss) 'Recall_cleaning'];
            [iChanClean3,iEvClean3,strClean] = jwCleanEEGevents_v01(RecallRawEEG,clfig,FIG_TITLE,clnWeights);
            disp(strClean);
            print(clfig,fullfile(outpath,[FIG_TITLE '.pdf']),'-dpdf','-bestfit');
            
            temChans=[temChans iChanClean1 iChanClean2 iChanClean3];
            
            temTrials=1:length([TEST_PROBE_Event(allssession==uniquesessions(iss))]);
            allcleanTrials=[allcleanTrials ismember(temTrials,iEvClean1) & ismember(temTrials,iEvClean2) & ismember(temTrials,iEvClean3)];
            SessionCleanTrials{iss}=temTrials(ismember(temTrials,iEvClean1) & ismember(temTrials,iEvClean2) & ismember(temTrials,iEvClean3));
            
            clear sessionGA ProbeRawEEG
        end
        allcleanChans=unique(temChans);
        
        
        %% only did the time frequency analysis on clean trials and clean
        % electrodes
        StudyPowerlogAllSession=[]; clear StudyPowerlog
        ProbePowerlogAllSession=[]; clear ProbePowerlog
        RecallPowerlogAllSession=[]; clear RecallPowerlog
        
        % Time-frequency analysis 
        for iss=1:length(uniquesessions)
            iss
            StudyCleanEEG{iss}=StudyRawEEGAll{iss}(allcleanChans,SessionCleanTrials{iss},:); % only include clean chans in both sessions for analysis.
            ProbeCleanEEG{iss}=ProbeRawEEGAll{iss}(allcleanChans,SessionCleanTrials{iss},:); % only include clean chans in both sessions for analysis.
            RecallCleanEEG{iss}=RecallRawEEGAll{iss}(allcleanChans,SessionCleanTrials{iss},:); % only include clean chans in both sessions for analysis
            
            clear ProbePowerlog RecallPowerlog StudyPowerlog
            
            parfor ic=1:size(allcleanChans,2)
                % get
                ic
                [~,WavePow_Raw] = multiphasevec3(waveletFreqs3,squeeze(StudyCleanEEG{iss}(ic,:,:)),resamp,waveletWidth);
                WavePow_log = 10*log10(WavePow_Raw);%log transformed raw power
                % zscoring normalize within each frequency band before averaging
                baseline = mean(WavePow_log(:,:,time>-1000 & time<-500),3);
                StudyPowerlog(ic,:,:,:)  = bsxfun(@rdivide,bsxfun(@minus,WavePow_log,mean(baseline)),std(baseline)); % bsxfun(@minus,WavePow_log,mean(baseline));%
                % electrodes, trial, timepoints.
                
                [~,WavePow_Raw] = multiphasevec3(waveletFreqs3,squeeze(ProbeCleanEEG{iss}(ic,:,:)),resamp,waveletWidth);
                WavePow_log = 10*log10(WavePow_Raw);%log transformed raw power
                % zscoring normalize within each frequency band before averaging
                baseline = mean(WavePow_log(:,:,time>-1000 & time<-500),3);
                ProbePowerlog(ic,:,:,:)  = bsxfun(@rdivide,bsxfun(@minus,WavePow_log,mean(baseline)),std(baseline)); % bsxfun(@minus,WavePow_log,mean(baseline));%
                % electrodes, trial, timepoints.
                
                [~,WavePow_Raw] = multiphasevec3(waveletFreqs3,squeeze(RecallCleanEEG{iss}(ic,:,:)),resamp,waveletWidth);
                WavePow_log = 10*log10(WavePow_Raw);%log transformed raw power
                % zscoring normalize within each frequency band before averaging
                RecallPowerlog(ic,:,:,:)  = bsxfun(@rdivide,bsxfun(@minus,WavePow_log,mean(baseline)),std(baseline)); % bsxfun(@minus,WavePow_log,mean(baseline));%
                % electrodes, trial, timepoints.
            end
            StudyPowerlogAllSession=[StudyPowerlogAllSession single(StudyPowerlog)];
            ProbePowerlogAllSession=[ProbePowerlogAllSession single(ProbePowerlog)];
            RecallPowerlogAllSession=[RecallPowerlogAllSession single(RecallPowerlog)];
            
            clear ProbePowerlog RecallPowerlog StudyPowerlog
        end
        
        
        SubjResults.allcleanChans=chan(allcleanChans);
        SubjResults.MeanRT=nanmean(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>0));
        SubjResults.MedianRT=nanmedian(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>0));
        SubjResults.allcleanTrials=allcleanTrials;
        
        cleanmask=SubjResults.allcleanTrials;
        AccuracyMask=[SubjTable(isub).PALTable.correct{1:end}]>0;
        AccuracyMask=AccuracyMask(cleanmask>0);
        
        
        SubjResults.Mean_Cor_StudyPowerlogAllSession=squeeze(mean(StudyPowerlogAllSession(:,AccuracyMask,:,time>=-1000 & time <=4000),2));
        SubjResults.Mean_Cor_ProbePowerlogAllSession=squeeze(mean(ProbePowerlogAllSession(:,AccuracyMask,:,time>=-1000 & time <=4000),2));
        SubjResults.Mean_Cor_RecallPowerlogAllSession=squeeze(mean(RecallPowerlogAllSession(:,AccuracyMask,:,time>=-4000 & time <=1000),2)); % chan, trial, time;
        
        SubjResults.Mean_InCor_StudyPowerlogAllSession=squeeze(mean(StudyPowerlogAllSession(:,~AccuracyMask,:,time>=-1000 & time <=4000),2));
        SubjResults.Mean_InCor_ProbePowerlogAllSession=squeeze(mean(ProbePowerlogAllSession(:,~AccuracyMask,:,time>=-1000 & time <=4000),2));
        SubjResults.Mean_InCor_RecallPowerlogAllSession=squeeze(mean(RecallPowerlogAllSession(:,~AccuracyMask,:,time>=-4000 & time <=1000),2)); % chan, trial, time;        
        
        clear AccuracyMask        
        
        %% moving average approach
        timewindowlength=0.500;overlappercept=0.90;Fs=resamp;
        
        clear StudyPowerlogAllSessionNormtime  ProbePowerlogAllSessionNormtime  RecallPowerlogAllSessionNormtime
        parfor ic=1:size(ProbePowerlogAllSession,1)
            ic
            [BinTimeCENTER, StudyPowerlogAllSessionNormtime(ic,:,:,:)]= windowed_average_with_freq(squeeze(StudyPowerlogAllSession(ic,:,:,:)),timewindowlength,overlappercept,Fs);
            [BinTimeCENTER, ProbePowerlogAllSessionNormtime(ic,:,:,:)]= windowed_average_with_freq(squeeze(ProbePowerlogAllSession(ic,:,:,:)),timewindowlength,overlappercept,Fs);
            [BinTimeCENTER, RecallPowerlogAllSessionNormtime(ic,:,:,:)]= windowed_average_with_freq(squeeze(RecallPowerlogAllSession(ic,:,:,:)),timewindowlength,overlappercept,Fs);
        end
        [BinTimeCENTER, RecallPowerlogAllSessionNormtime(1,:,:,:)]= windowed_average_with_freq(squeeze(RecallPowerlogAllSession(1,:,:,:)),timewindowlength,overlappercept,Fs);
        
        BinTimeCENTER=(BinTimeCENTER*1000)-5000;
        BinTime=BinTimeCENTER;
        
        clear Low_Study_DataToDecode High_Study_DataToDecode  Low_Probe_DataToDecode Low_Recall_DataToDecode High_Probe_DataToDecode High_Recall_DataToDecode Probe_DataToDecode Recall_DataToDecode
        clear StudyPowerlogAllSessionNormtime_b  ProbePowerlogAllSessionNormtime_b  RecallPowerlogAllSessionNormtime_b
        
        %freq binning
        for itick=1:length(freqBandYticks)-1
            StudyPowerlogAllSessionNormtime_b(:,:,itick,:) = mean(StudyPowerlogAllSessionNormtime(:,:,waveletFreqs3>=freqBandYticks(itick) & waveletFreqs3<=freqBandYticks(itick+1),:),3);
            ProbePowerlogAllSessionNormtime_b(:,:,itick,:) = mean(ProbePowerlogAllSessionNormtime(:,:,waveletFreqs3>=freqBandYticks(itick) & waveletFreqs3<=freqBandYticks(itick+1),:),3);
            RecallPowerlogAllSessionNormtime_b(:,:,itick,:) = mean(RecallPowerlogAllSessionNormtime(:,:,waveletFreqs3>=freqBandYticks(itick) & waveletFreqs3<=freqBandYticks(itick+1),:),3);
        end
        
        
        clear Study_DataToDecode  Probe_DataToDecode  Recall_DataToDecode
        for itrl=1:size(ProbePowerlogAllSessionNormtime,2)
            for itime=1:size(ProbePowerlogAllSessionNormtime,4)
                Study_DataToDecode(itrl,:,itime)=reshape(StudyPowerlogAllSessionNormtime_b(:,itrl,:,itime),size(allcleanChans,2)*size(StudyPowerlogAllSessionNormtime_b,3),1);
                Probe_DataToDecode(itrl,:,itime)=reshape(ProbePowerlogAllSessionNormtime_b(:,itrl,:,itime),size(allcleanChans,2)*size(ProbePowerlogAllSessionNormtime_b,3),1);
                Recall_DataToDecode(itrl,:,itime)=reshape(RecallPowerlogAllSessionNormtime_b(:,itrl,:,itime),size(allcleanChans,2)*size(RecallPowerlogAllSessionNormtime_b,3),1);
            end
        end
        
        
        % all low<=16 and high >=32 frequency
        for itrl=1:size(ProbePowerlogAllSessionNormtime,2)
            for itime=1:size(ProbePowerlogAllSessionNormtime,4)
                Low_Study_DataToDecode(itrl,:,itime)=reshape(StudyPowerlogAllSessionNormtime(:,itrl,waveletFreqs3<16.5,itime),size(allcleanChans,2)*sum(waveletFreqs3<16.5),1);
                Low_Probe_DataToDecode(itrl,:,itime)=reshape(ProbePowerlogAllSessionNormtime(:,itrl,waveletFreqs3<16.5,itime),size(allcleanChans,2)*sum(waveletFreqs3<16.5),1);
                Low_Recall_DataToDecode(itrl,:,itime)=reshape(RecallPowerlogAllSessionNormtime(:,itrl,waveletFreqs3<16.5,itime),size(allcleanChans,2)*sum(waveletFreqs3<16.5),1);
                
                High_Study_DataToDecode(itrl,:,itime)=reshape(StudyPowerlogAllSessionNormtime(:,itrl,waveletFreqs3>32,itime),size(allcleanChans,2)*sum(waveletFreqs3>32),1);
                High_Probe_DataToDecode(itrl,:,itime)=reshape(ProbePowerlogAllSessionNormtime(:,itrl,waveletFreqs3>32,itime),size(allcleanChans,2)*sum(waveletFreqs3>32),1);
                High_Recall_DataToDecode(itrl,:,itime)=reshape(RecallPowerlogAllSessionNormtime(:,itrl,waveletFreqs3>32,itime),size(allcleanChans,2)*sum(waveletFreqs3>32),1);
            end
        end
        
        SubjResults.Low_Study_DataToDecode=Low_Study_DataToDecode;
        SubjResults.Low_Probe_DataToDecode=Low_Probe_DataToDecode;
        SubjResults.Low_Recall_DataToDecode=Low_Recall_DataToDecode;
        
        SubjResults.High_Study_DataToDecode=High_Study_DataToDecode;
        SubjResults.High_Probe_DataToDecode=High_Probe_DataToDecode;
        SubjResults.High_Recall_DataToDecode=High_Recall_DataToDecode;
        
        SubjResults.Study_DataToDecode=Study_DataToDecode;
        SubjResults.Probe_DataToDecode=Probe_DataToDecode;
        SubjResults.Recall_DataToDecode=Recall_DataToDecode;
        
        SubjResults.BinTime=BinTimeCENTER;
        SubjResults.ProbeCleanEEG=ProbeCleanEEG;
        SubjResults.RecallCleanEEG=RecallCleanEEG;
        
        % get accuracy mask and memorability mask based on the clean trials and clearn electrodes;
        AccuracyMask=[SubjTable(isub).PALTable.correct{1:end}]>0;
        HighMemMask = WordtoAnalyzeTReduced.ResponseMemorability{isubtem} > medianMemor; %  >= MemPercent(3); %
        LowMemMask = WordtoAnalyzeTReduced.ResponseMemorability{isubtem} < medianMemor; %  = MemPercent(2);  %
        
        MemorabilityScore =  WordtoAnalyzeTReduced.ResponseMemorability_z{isubtem};
        MemorabilityScore = MemorabilityScore(SubjResults.allcleanTrials>0);
        
        HighMemMask=HighMemMask(allcleanTrials>0);
        LowMemMask=LowMemMask(allcleanTrials>0);
        AccuracyMask=AccuracyMask(allcleanTrials>0);
        
        SubjResults.AccuracyMask=AccuracyMask;
        SubjResults.HighMemMask=HighMemMask;
        SubjResults.LowMemMask=LowMemMask;
        SubjResults.MemorabilityScore=MemorabilityScore;
        SubjResults.cleanRTs=SubjTable(isub).RTtemp(allcleanTrials>0);
        
        savefilename=fullfile(outpath,['SubjResults_ATL' num2str(isub) '.mat']);
        save(savefilename,'SubjResults','-v7.3');
        
    end
    
end

