clear all
%add toolboxes
addpath(genpath('/Volumes/Zane/Matlab/eeg_toolbox/trunk'));
addpath(genpath('/Volumes/Zane/Matlab/dungeon_toolbox_17a'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/EEG_Preprocessing'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/IEM_ester'));
addpath(genpath('/Volumes/Zane/NIH_NINDS/Data_InProgress/SpecFuN'));

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

% freqBandYticksSplit=[exp(linspace(log(2),log(16),16));exp(linspace(log(70),log(150),16))]; % only bin the data into 3 bins.

freqBandYticksSplit=[exp(linspace(log(4),log(16),10));exp(linspace(log(70),log(150),10))]; % only bin the data into 3 bins.

waveletFreqs2 = reshape(freqBandYticksSplit',1,20);
   
waveletFreqs3 = exp(linspace(log(4),log(150),25));

freqBandYticks = [ 4   8   16   32   70   150]; % overwrite the old ticks.
WordtoAnalyzeTReduced
mkdir('ATL_3Event_500')
curpath=pwd;
outpath=fullfile(curpath,'ATL_3Event_500');


% loop through subjects.
%%

for isubtem=2:length(WordtoAnalyzeTReduced.subID)
    
    
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
   
    % Time-frequency
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
    
    
    %%
%     StudyPowerlogAllSession_s = smoothdata(StudyPowerlogAllSession,4,'gaussian',250);    
%     ProbePowerlogAllSession_s = smoothdata(ProbePowerlogAllSession,4,'gaussian',250);
%     RecallPowerlogAllSession_s = smoothdata(RecallPowerlogAllSession,4,'gaussian',250);
%        
    SubjResults.allcleanChans=chan(allcleanChans);
    SubjResults.MeanRT=nanmean(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>50));
    SubjResults.MedianRT=nanmedian(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>0));
    SubjResults.allcleanTrials=allcleanTrials;  

    %% moving average approach 
    timewindowlength=0.500;overlappercept=0.80;Fs=resamp;
    
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
    freqBandYticks=[ 4 8 16 32 70 150];
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
    
    
    % all low frequency 
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



%% decode using SubjResults

allsub=[26:64 66];
rootEEGdir = '/Volumes/Zane/NIH_FRNU_ROOT';                      %office-local

load('SubjTable_palRAMword.mat')
load('PAL_Memo_PALRAM.mat')

ToAnalyzeSub=allsub;
% PAL_Memo.Combinedmemorability = zscore(PAL_Memo.Probememorability) + (PAL_Memo.Responsememorability - nanmean(PAL_Memo.Responsememorability))./nanstd(PAL_Memo.Responsememorability);
% PAL_Memo.Combinedmemorability = PAL_Memo.Combinedmemorability/2;

medianMemor = median(PAL_Memo.Responsememorability);

datapath=outpath;

load('WordtoAnalyzeTReduced.mat')

fprintf('loading SubjResults_ATL...\n')

ATLIDs=WordtoAnalyzeTReduced.subID;
accutrials=[WordtoAnalyzeTReduced.WordCount{1:end}].*[WordtoAnalyzeTReduced.meanACC{1:end}];
inacctrial=[WordtoAnalyzeTReduced.WordCount{1:end}].*(1-[WordtoAnalyzeTReduced.meanACC{1:end}]);
ATLIDs=ATLIDs(accutrials >= 30 & inacctrial>=30 );

%%
for isubtem = 2:length(ATLIDs)
    
    isub=ATLIDs(isubtem);
    
    isub
     savefilename=fullfile(datapath,['SubjResults_ATL' num2str(isub) '.mat']);
    
    
    if exist(savefilename)>0
        
        clear SubjResults
        
        load(savefilename);
        

        %list based analysis
        Curlist=[SubjTable(isub).PALTable.list{1:end}];
        sessions=[[SubjTable(isub).alignedEvents.RecallEvent.session]];
        
%         ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.probe)

        for iw=1:length(SubjTable(isub).PALTable.probe)
            ProbememPreict(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.probe(iw)));
            ExpectmemPreict(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.expected(iw)));
        
            Item1mem(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.item1(iw)));
            Item2mem(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.item2(iw)));
        end

        AverageStudyMem = (ProbememPreict+ExpectmemPreict)/2;
        
        uniquelist = unique(Curlist);
        uniquesession=unique(sessions);
        
        Curlistupdate=[];
        for iss=1:length(unique(uniquesession))
                
            if iss==1
            Curlistupdate = [Curlistupdate Curlist(sessions==uniquesession(iss))];
            else
                
            Curlistupdate = [Curlistupdate Curlist(sessions==uniquesession(iss))*100];
            end
        
        end
        
        uniquelist=unique(Curlistupdate);
        
        AccuracyMask=[SubjTable(isub).PALTable.correct{1:end}]>0;
        RTs=[SubjTable(isub).PALTable.RT{1:end}];
        RTs(RTs<0)=nan;
        

         
        clear listmemscore  ListMem  listmemscore  MemlistACC RTnanmeans  MemlistACC
        for ils=1:length(uniquelist)
            ListMem(Curlistupdate==uniquelist(ils)) = mean(AverageStudyMem(Curlistupdate==uniquelist(ils)));
            RTnanmeans(ils)=nanmean(log(RTs(Curlistupdate==uniquelist(ils))));
            listmemscore(ils)=nanmean(AverageStudyMem(Curlistupdate==uniquelist(ils)));
            MemlistACC(ils)=nanmean(AccuracyMask(Curlistupdate==uniquelist(ils)));
        end
        nonnanin=~isnan(RTnanmeans);
        nanmeanRTall(isub)=nanmean(RTs);
        
        correlationRTmemList(isub) = corr((RTnanmeans(nonnanin))',listmemscore(nonnanin)','type','pearson');
        correlationACCmemList(isub) = corr(MemlistACC(nonnanin)',listmemscore(nonnanin)','type','pearson');

 
        cleanmask=SubjResults.allcleanTrials;
        AccuracyMask=[SubjTable(isub).PALTable.correct{1:end}]>0;
        
        AverageStudyMem = AverageStudyMem(cleanmask>0);
        ExpectmemPreict = ExpectmemPreict(cleanmask>0);
        ProbememPreict=ProbememPreict(cleanmask>0);
        ListMem=ListMem(cleanmask>0);

%         HighMemMask = WordtoAnalyzeTReduced.ResponseMemorability{isubtem} > medianMemor; %  >= MemPercent(3); %
%         LowMemMask = WordtoAnalyzeTReduced.ResponseMemorability{isubtem} < medianMemor; %  = MemPercent(2);  %
        
        AccuracyMask=AccuracyMask(cleanmask>0);
%         HighACCIndex= find(HighMemMask==1);
%         LowACCIndex= find(LowMemMask==1);
   
       %        subtriallength(isub,:)= [length(HighACCIndex) length(LowACCIndex)  length(HighWrongIndex) length(LowWrongIndex) ]



       TrialsPloted=[sum(AccuracyMask==1 & ExpectmemPreict<=medianMemor) sum(AccuracyMask==1 & ExpectmemPreict>medianMemor)];
       min([sum(AccuracyMask==1 & ExpectmemPreict<=medianMemor) sum(AccuracyMask==1 & ExpectmemPreict>medianMemor)])

% % RSA
% clear MemoModelTarget  MemoModelstudy
% %model % based on study
% HighMemMask = AverageStudyMem > medianMemor; %  >= MemPercent(3); %
% LowMemMask = AverageStudyMem < medianMemor; %  = MemPercent(2);  %
% MemoModelstudy=zeros(length(HighMemMask),length(HighMemMask));
% for i=1:length(HighMemMask)
%     for j=1:length(HighMemMask)
%         if HighMemMask(i) == 0 & HighMemMask(j) ==0 % low and high memory are highly similiar within each group
%             MemoModelstudy(i,j)=1;
%         elseif HighMemMask(i) == 1 & HighMemMask(j) ==1
%             MemoModelstudy(i,j)=0;
%         else
%             MemoModelstudy(i,j)=0;
%         end
%     end
% end
% 
%             
% %model % based on expected word
% clear MemoModelTarget
% MemPercent2=prctile(ListMem,[0 1/3 2/3 1]*100);
% HighMemMask = ListMem >= MemPercent2(3); % medianMemor; % 
% LowMemMask = ListMem <= MemPercent2(2);  %%  
% Remaintrials= find(HighMemMask | LowMemMask);
% MemoModelTarget=zeros(length(Remaintrials),length(Remaintrials));
% for i=1:length(Remaintrials)
%     for j=1:length(Remaintrials)
%         if HighMemMask(Remaintrials(i)) == 0 & HighMemMask(Remaintrials(j)) ==0 % low and high memory are highly similiar within each group
%             MemoModelTarget(i,j)=1;
%         elseif HighMemMask(Remaintrials(i)) == 1 & HighMemMask(Remaintrials(j)) ==1
%             MemoModelTarget(i,j)=0;
%         else
%             MemoModelTarget(i,j)=0;
%         end
%     end
% end
%             
%             
% %ACC
% clear accMemoModel
% for i=1:length(AccuracyMask)
%     for j=1:length(AccuracyMask)
%         if AccuracyMask(i) == 0 & AccuracyMask(j) ==0 % ACC model, accurate trials are similiar, inaccurate trials are similiar
%             accMemoModel(i,j)=1;
%         elseif AccuracyMask(i) == 1 & AccuracyMask(j) ==1
%             accMemoModel(i,j)=0;
%         else
%             accMemoModel(i,j)=0.5;  % in accurate vs acurrate highest diff model, accurate trials are similiar, inaccurate trials are similiar
%         end
%     end
% end
% 
%             
% clear a1 a2 b1  Ptempdata b2 High_a2 High_a1 low_a2 low_a1 STUDY_Mem  Probe_ACC Recall_Mem Recall_ACC
% 
% reshapeInxT=find(tril(ones(size(MemoModelTarget,1)),-1));
% reshapeInx=find(tril(ones(size(MemoModelTarget,1)),-1));
% 
% for itime=1:size(SubjResults.High_Study_DataToDecode,3)
%     
%     RSA_Study = squareform(pdist(SubjResults.Low_Study_DataToDecode(:,:,itime),'cosine'));
%     RSA_Probe = squareform(pdist(SubjResults.Low_Probe_DataToDecode(Remaintrials,:,itime),'cosine'));
%     RSA_Recall = squareform(pdist(SubjResults.Low_Recall_DataToDecode(:,:,itime),'cosine'));
%     
%     STUDY_Mem(itime)=corr(RSA_Study(reshapeInx),MemoModelstudy(reshapeInx),'type','Spearman');
%     STUDY_ACC(itime)=corr(RSA_Study(reshapeInx),accMemoModel(reshapeInx),'type','Spearman');
%     
%     Probe_Mem(itime)=corr(RSA_Probe(reshapeInxT),MemoModelTarget(reshapeInxT),'type','Spearman');
% %     Probe_ACC(itime)=corr(RSA_Probe(reshapeInx),accMemoModel(reshapeInx),'type','Spearman');
%     
%     Recall_Mem(itime)=corr(RSA_Recall(reshapeInx),MemoModelTarget(reshapeInx),'type','Spearman');
%     Recall_ACC(itime)=corr(RSA_Recall(reshapeInx),accMemoModel(reshapeInx),'type','Spearman');    
% end
%             
%             
% 
% RSA_Probe_Mem_normbased(isubtem,:)=Probe_Mem;
% RSA_Probe_ACC(isubtem,:)=Probe_ACC;
% 
% RSA_Study_Mem_normbased(isubtem,:)=STUDY_Mem;
% RSA_Study_ACC(isubtem,:)=STUDY_ACC;
% 
% RSA_Recall_Mem_normbased(isubtem,:)=Recall_Mem;
% RSA_Recall_ACC(isubtem,:)=Recall_ACC;
% 

%%

    %Change the dimensions to be feature X time X event
    encMat = SubjResults.Study_DataToDecode;
    recMat = SubjResults.Recall_DataToDecode;
    numEvents=size(encMat,1)
    clear output;
    for event = 1:numEvents
        encData = squeeze(encMat(event,:,:));
        recData = squeeze(recMat(event,:,:));
        s = (encData')*recData;
        normEnc = sqrt(sum(encData.^2));
        normRec = sqrt(sum(recData.^2));
        norms = (normEnc')*normRec;
        output(event,:,:) = s./norms;
    end

    STUDY_Probe_Similarity=output;
    
    

BinTime=SubjResults.BinTime;  
% 
% 
% %Obtain across subject average
% AccurateTrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1,:,:)));
% Gl_Accurate_highmem_TrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 & ExpectmemPreict>medianMemor,:,:)));
% Gl_Accurate_lowmem_TrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 & ExpectmemPreict<=medianMemor,:,:)));
% InAccurateTrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==0,:,:)));
% 
% lc_medianMemor(isubtem)= median(ExpectmemPreict(AccuracyMask==1));
% lc_Accurate_highmem_TrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 & ExpectmemPreict>lc_medianMemor(isubtem),:,:)));
% lc_Accurate_lowmem_TrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 & ExpectmemPreict<=lc_medianMemor(isubtem),:,:)));
% 
% MedianRT(isubtem)=prctile(SubjResults.cleanRTs(AccuracyMask==1),50);
% MeanRT(isubtem)=mean(SubjResults.cleanRTs(AccuracyMask==1));
% % 
% 
% 
% figure;
% xyimagesc(BinTime,BinTime,squeeze(((lc_Accurate_lowmem_TrialsSim(isubtem,:,:))))); %
% axis xy
% 
% 
% figure;
% xyimagesc(BinTime,BinTime,squeeze(mean(InAccurateTrialsSim))'); %
% 
% figure;
% xyimagesc(BinTime,BinTime,squeeze(mean(AccurateTrialsSim))'-squeeze(mean(InAccurateTrialsSim))'); %
% 
% 
% 
%%
% STUDY_Probe_Similarity1=STUDY_Probe_Similarity- mean(mean(STUDY_Probe_Similarity(:,BinTime<-500 & BinTime>-1000,BinTime>-4000 & BinTime<-3500),2),3);
BinTime=SubjResults.BinTime;  

STUDY_Probe_Similarity1=STUDY_Probe_Similarity;
AccuracyMask
ExpectmemPreict

minRT=prctile(SubjResults.cleanRTs(AccuracyMask==1),50);

h=figure;
subplot(2,4,1);
xyimagesc(BinTime,BinTime,squeeze(mean(STUDY_Probe_Similarity1(AccuracyMask==1 ,:,:)))); %
axis([-4000 1000 -1000 4000]);
hold on;plot([-minRT -minRT],[-1000 4000],'r-','linewidth',2)
hold on;plot([0 0],[-1000 4000],'r.-','linewidth',1)
ylabel('Study Onset');caxis([-0.1 0.1]) 
xlabel('Retrieval Onset')
title('accurate Trials')

subplot(2,4,2);
xyimagesc(BinTime,BinTime,squeeze(mean(STUDY_Probe_Similarity1(AccuracyMask==0 ,:,:)))); %& ExpectmemPreict<=medianMemor
axis([-4000 1000 -1000 4000]);
hold on;plot([-minRT -minRT],[-1000 4000],'r-','linewidth',2)
hold on;plot([0 0],[-1000 4000],'r.-','linewidth',1)
ylabel('Study Onset');caxis([-0.1 0.1]) 
xlabel('Retrieval Onset')
title('Inaccurate Trials')

subplot(2,4,3);
xyimagesc(BinTime,BinTime,squeeze(mean(STUDY_Probe_Similarity1(AccuracyMask==1 ,:,:)))-squeeze(mean(STUDY_Probe_Similarity1(AccuracyMask==0 ,:,:))));
caxis([-0.1 0.1]) 
hold on
[h,p]=ttest(STUDY_Probe_Similarity1(AccuracyMask==1 ,:,:)-squeeze(mean(mean(STUDY_Probe_Similarity1(AccuracyMask==1,BinTime<-3000 & BinTime>-4000, BinTime<0 & BinTime>-1000),2),3)),0,'tail','right','alpha',0.05); 
% [h,p]=ttest2(STUDY_Probe_Similarity1(AccuracyMask==1 ,:,:),STUDY_Probe_Similarity1(AccuracyMask==0 ,:,:),'tail','right','alpha',0.05); 
xyimagesc(BinTime,BinTime,squeeze(h));
axis([-4000 1000 -1000 4000]);
hold on;plot([-minRT -minRT],[-1000 4000],'r-','linewidth',2)
hold on;plot([0 0],[-1000 4000],'r.-','linewidth',1)

TimeROI=BinTime(BinTime>0);
TimeSig=(mean(squeeze(h(:,BinTime>0,BinTime<1000 & BinTime>=-minRT))>0,2)>0);
TimeSig=TimeROI(TimeSig>0);
hold on;plot([-5000 5000],[min(TimeSig) min(TimeSig)],'r-','linewidth',2);
hold on;plot([-5000 5000],[max(TimeSig) max(TimeSig)],'r-','linewidth',2);

A = squeeze(mean(STUDY_Probe_Similarity1(:,BinTime>=min(10000,min(TimeSig)) & BinTime<=max(2,max(TimeSig)),:),2))

% A = squeeze(mean(STUDY_Probe_Similarity1(:,BinTime>=min(TimeSig) & BinTime<=max(TimeSig),:),2))
xlabel('Study Onset');
ylabel('Retrieval Onset')

subplot(2,4,5);
xyimagesc(BinTime,BinTime,squeeze(mean(STUDY_Probe_Similarity1(AccuracyMask==1 & ExpectmemPreict>medianMemor ,:,:)))); %
axis([-4000 1000 -1000 4000]);
hold on;plot([-minRT -minRT],[-1000 4000],'r-','linewidth',2)
hold on;plot([0 0],[-1000 4000],'r.-','linewidth',1)
ylabel('Study Onset');caxis([-0.1 0.1]) 
xlabel('Retrieval Onset')
title('HighMem ACC Trials')


subplot(2,4,6);
xyimagesc(BinTime,BinTime,squeeze(mean(STUDY_Probe_Similarity1(AccuracyMask==1 & ExpectmemPreict<=medianMemor ,:,:)))); %
axis([-4000 1000 -1000 4000]);
hold on;plot([-minRT -minRT],[-1000 4000],'r-','linewidth',2)
hold on;plot([0 0],[-1000 4000],'r.-','linewidth',1)
ylabel('Study Onset');caxis([-0.1 0.1]) 
xlabel('Retrieval Onset')
title('LowMem ACC Trials')


% Ar=A(:,BinTime<1000 & BinTime >= -minRT);
AZ=A;%mean(A,2);
subplot(2,4,7);
plot(BinTime,mean(AZ(AccuracyMask==1 & ExpectmemPreict>medianMemor,:)));
hold on;
plot(BinTime,mean(AZ(AccuracyMask==1 & ExpectmemPreict<=medianMemor,:)))
plot([0 0],[-0.2 0.2],'k-');axis([min(BinTime) 1000 -0.12 0.12])
plot([-minRT -minRT],[-0.2 0.2],'k.-');axis([min(BinTime) 1000 -0.12 0.12])
legend('HighMem','LowMem','recall onset','location','southeast')
title('Averaged Similarity')
xlabel('Percentage of Recall Time')


% AccuracyMask==1 & ExpectmemPreict>medianMemor
RoundBinTime=round(BinTime/25)*25;
% normalized by trial-by-trial reaction time. 
AllAccRT=SubjResults.cleanRTs;
clear AZpercent
for itr=1:size(AllAccRT,2)
    if AllAccRT(itr)==-999
        currentRT=-3450;
    else
        currentRT=-AllAccRT(itr);
        currentRT=round(currentRT/25)*25;
    end
    y = AZ(itr,RoundBinTime>=currentRT & RoundBinTime<=-currentRT);
    x = linspace(-100, 100, length(y));    
    xi = [-100:1:100]; %linspace(-100, 100, 20);    
    
    AZpercent(itr,:)= interp1(x,y,xi);
end

% pertime=linspace(-100, 100, 20);
pertime=[-100:1:100]

subplot(2,4,8);
plot(pertime,mean(AZpercent(AccuracyMask==1 & ExpectmemPreict>medianMemor,:)));
hold on;
plot(pertime,mean(AZpercent(AccuracyMask==1 & ExpectmemPreict<=medianMemor,:)))
plot([0 0],[-0.2 0.2],'k-');axis([min(pertime) 50 -0.08 0.08])
plot([-minRT -minRT],[-0.2 0.2],'k.-');axis([min(pertime) 50 -0.12 0.12])
legend('HighMem','LowMem','recall onset','location','southeast')
title('Averaged Similarity')
xlabel('Percentage of Recall Time')

    h=gcf;
    set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
    printfilename=fullfile(outpath,['Decode_Sub' num2str(isub) '_Reinstate_Mem.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
    
corr(mean(AZpercent(AccuracyMask==1,pertime<0),2),ExpectmemPreict(AccuracyMask==1))

% Az=zscore(A,0,1);
% bx=find(AccuracyMask==1);
% [x,sortedBX]=sort(ExpectmemPreict(bx));
% bx=bx(sortedBX);
% figure;
% for ib=1:length(bx)
%    scaledA(ib,:) = rescale(Az(bx(ib),BinTime<1000),0,1)
%    subplot(length(bx),1,ib) 
% 
%    plot(BinTime(BinTime<1000),scaledA(ib,:));
% end

% & ExpectmemPreict>medianMemor
% & ExpectmemPreict<=medianMemor
% DecodeMemLabels=AverageStudyMem>medianMemor;
% alltimetasktimeinx=1:size(SubjResults.Low_Study_DataToDecode,3);
% DecodeMemInx=1:size(SubjResults.Low_Study_DataToDecode,1);
% 
% [Memo.Low_Study_DecodeACC(iter,:), Memo.Low_Study_DecodeAUC(iter,:),~,Memo.Low_Study_DecodeAUC_sh(iter,:)] = ...
%     DecodeSVM_PAL(SubjResults.High_Study_DataToDecode, DecodeMemInx, alltimetasktimeinx, DecodeMemLabels);
% 
% LowMemState=mean(squeeze(mean(SubjResults.Low_Study_DataToDecode(AverageStudyMem<medianMemor,:,BinTime>2000 & BinTime<3000))),2);
% HighMemState=mean(squeeze(mean(SubjResults.Low_Study_DataToDecode(AverageStudyMem>medianMemor,:,BinTime>2000 & BinTime<3000))),2);
% 
% TimeBin=[-1000:250:2000];
% 
% CurProbeData=SubjResults.Low_Probe_DataToDecode;
% for itime=1:size(CurProbeData,3)
%     for item=1:size(CurProbeData,1)
%         LowMemState_Corr(item,itime)=1-pdist([squeeze(CurProbeData(item,:,itime));LowMemState'],'spearman');
%         HighMemState_Corr(item,itime)=1-pdist([squeeze(CurProbeData(item,:,itime));HighMemState'],'spearman');
%     end
% end
% HighMem=HighMemState_Corr-LowMemState_Corr;
% for it=1:length(TimeBin)-1
%     PredictTargetMem(it)=corr(mean(HighMem(:,BinTime>TimeBin(it) & BinTime<TimeBin(it+1)),2),ExpectmemPreict','type','spearman')
% end

    end
    
            
            
            %%
           
        
       
       end
%     end

%%

acount = [WordtoAnalyzeTReduced.meanACC{1:end}].*[WordtoAnalyzeTReduced.WordCount{1:end}];
includeATL =  acount>7;



%%
plotid=(find(sum(subTRIALS<9,2)==0));

% plotid=1:length(plotid);

for ip=1:length(plotid)
    
    isubtem=plotid(ip);
  subTRIALS(isubtem,:) =   [size(NeuSim(isubtem).low_a1,1) size(NeuSim(isubtem).High_a1,1)];
    
figure(isubtem);


subplot(1,2,1);
plot(BinTime,mean(NeuSim(isubtem).low_a1))
hold on;plot(BinTime,mean((NeuSim(isubtem).High_a1)))

subplot(1,2,2);
plot(BinTime,mean(NeuSim(isubtem).low_a2))
hold on;plot(BinTime,mean((NeuSim(isubtem).High_a2)))

meanLowMem_low(ip,:)=mean(NeuSim(isubtem).low_a2);
meanHighMem_low(ip,:)=mean(NeuSim(isubtem).High_a2);

end

%%





    clear Memo Corr

    iter=4;
    % Pass subject specific value 
    alltimetasktimeinx=SubjResults.alltimetasktimeinx;
    Low_Probe_DataToDecode=SubjResults.Low_Probe_DataToDecode;
    Low_Recall_DataToDecode=SubjResults.Low_Recall_DataToDecode;
    High_Probe_DataToDecode=SubjResults.High_Probe_DataToDecode;
    High_Recall_DataToDecode=SubjResults.High_Recall_DataToDecode;    
    
    % define memorability mask in the iteration

    MemPercent=prctile(PAL_Memo.Responsememorability,[0 33.33 66.67 100]);
    medianMemor = median(PAL_Memo.Responsememorability);
    HighMemMask = WordtoAnalyzeTReduced.ResponseMemorability{isubtem} > medianMemor ;%   >= MemPercent(3); %   ;
    LowMemMask = WordtoAnalyzeTReduced.ResponseMemorability{isubtem} < medianMemor ;%   <= MemPercent(2);  %    medianMemor; % 
    
    HighMemInx=find(HighMemMask);
    LowMemInx=find(LowMemMask);
    minlengthMem= min(length(HighMemInx),length(LowMemInx))-1;
    DecodeMemInx=[randsample(HighMemInx,minlengthMem,false),randsample(LowMemInx,minlengthMem,false)];
    DecodeMemLabels=[ones(1,minlengthMem),zeros(1,minlengthMem)];
 
    % define accuracy mask in the iteration
    accInx=find(SubjResults.AccuracyMask);
    InaccInx=find(~SubjResults.AccuracyMask);
    minlengthMem= min(length(accInx),length(InaccInx))-1;
    DecodeACCInx=[randsample(accInx,minlengthMem,false),randsample(InaccInx,minlengthMem,false)];
    DecodeACCLabels=[ones(1,minlengthMem),zeros(1,minlengthMem)];

    [Memo.Low_Probe_DecodeACC(iter,:), Memo.Low_Probe_DecodeAUC(iter,:),~,Memo.Low_Probe_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(Low_Probe_DataToDecode, DecodeMemInx, alltimetasktimeinx, DecodeMemLabels);
    [Memo.Low_Recall_DecodeACC(iter,:), Memo.Low_Recall_DecodeAUC(iter,:),~,Memo.Low_Recall_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(Low_Recall_DataToDecode, DecodeMemInx, alltimetasktimeinx, DecodeMemLabels);
    [Memo.High_Probe_DecodeACC(iter,:), Memo.High_Probe_DecodeAUC(iter,:),~,Memo.High_Probe_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(High_Probe_DataToDecode, DecodeMemInx, alltimetasktimeinx, DecodeMemLabels);
    [Memo.High_Recall_DecodeACC(iter,:), Memo.High_Recall_DecodeAUC(iter,:),~,Memo.High_Recall_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(High_Recall_DataToDecode, DecodeMemInx, alltimetasktimeinx, DecodeMemLabels);
     
    [Corr.Low_Probe_DecodeACC(iter,:), Corr.Low_Probe_DecodeAUC(iter,:),~,Corr.Low_Probe_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(Low_Probe_DataToDecode, DecodeACCInx, alltimetasktimeinx, DecodeACCLabels);
    [Corr.Low_Recall_DecodeACC(iter,:), Corr.Low_Recall_DecodeAUC(iter,:),~,Corr.Low_Recall_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(Low_Recall_DataToDecode, DecodeACCInx, alltimetasktimeinx, DecodeACCLabels);
    [Corr.High_Probe_DecodeACC(iter,:), Corr.High_Probe_DecodeAUC(iter,:),~,Corr.High_Probe_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(High_Probe_DataToDecode, DecodeACCInx, alltimetasktimeinx, DecodeACCLabels);
    [Corr.High_Recall_DecodeACC(iter,:), Corr.High_Recall_DecodeAUC(iter,:),~, Corr.High_Recall_DecodeAUC_sh(iter,:)] = DecodeSVM_PAL(High_Recall_DataToDecode, DecodeACCInx, alltimetasktimeinx, DecodeACCLabels);    
   

SubjResults.Memo=Memo;
SubjResults.Corr=Corr;




figure;
subplot(1,2,1);
hold on;plot(BinTime(tasktime),mean(Memo.Low_Probe_DecodeAUC));
hold on;plot(BinTime(tasktime),(Corr.Low_Probe_DecodeAUC));
hold on;plot(BinTime(tasktime),(low_PAUC_Memorability_sh));
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   
legend('memorability','memory accuracy','location','southoutside');

subplot(1,2,2);
plot(BinTime(tasktime),mean(low_R_AUC_Memorability));
hold on;plot(BinTime(tasktime),mean(low_R_AUC_accuracy));
hold on;plot(BinTime(tasktime),mean(low_R_AUC_Memorability_sh));
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   
legend('memorability','memory accuracy','location','southoutside');

h=gcf;

set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
printfilename=fullfile(outpath,['Decode_Sub' num2str(isub) '_SVM_allfreatures.pdf']);
print(h,printfilename,'-dpdf','-bestfit');
        




%%
includeID=WordtoAnalyzeTReduced.subID([WordtoAnalyzeTReduced.meanACC{1:end}]>=0.15);

a = SubResults_PACC_Memorab(includeID,:)
b = SubResults_PAUC_Memorab(find(sum(SubResults_PACC_Memorab'>0)),:)
% figure;
hold on;
prettyplot(BinTime(tasktime),mean(a),std(a)/sqrt(size(a,1)),[0  1 0],1)
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   


addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/waveform_comparisons_matlab'));
[acorr,kmin]=gutautocorr(a);
[minlen]=gsgutsims_paried_twotailed(size(a,1),size(a,2),0.05,acorr,2000,0.05);


for i=1:2000
    meana(i,:)=mean(a(randsample(1:size(a,1),size(a,1),'true'),:),1);
end
sigtest=(1-mean(meana>0.5)<0.05);
xtime=BinTime(tasktime);
[B, N, BI]=RunLength(sigtest');
plotix=(BI(B==1 & N>minlen));
for iplox=1:length(plotix)
    if xtime(plotix(iplox))>0
        hold on;plot(xtime(plotix(iplox):plotix(iplox)+N(BI==plotix(iplox))-1),ones(N(BI==plotix(iplox)),1)*0.5,'bo');
    end
end

axis([-500 1500 0.4 0.7])


%% === NOT USE FOR NOW. 


for i=1:10
%% SVM decoding
    nter=10;
    clear low_PAUC_Memorability low_PAUC_Accuracy low_RAUC_Memorability low_RAUC_Accuracy itime
    clear high_PAUC_Memorability high_PAUC_Accfwuracy high_RAUC_Memorability high_RAUC_Accuracy itime
    
    alltimetasktimeinx=find(tasktime);
    for iter=1:nter
        parfor itime=1:sum(tasktime)
            itime
            %%===low freq
            Ptempdata=squeeze(Low_Probe_DataToDecode(:,:,alltimetasktimeinx(itime)));
            Rtempdata=squeeze(Low_Recall_DataToDecode(:,:,alltimetasktimeinx(itime)));
            
            % probe
            T=array2table([Ptempdata HighMemMask']);
            [~,~,low_PAUC_Memorability(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_svm_guassian(T); %
            
            T=array2table([Ptempdata Shuffle(HighMemMask)']);
            [~,~,low_PAUC_Memorability_sh(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_svm_guassian(T); %
           
            T=array2table([Ptempdata (AccuracyMask)']);
            [~, ~,low_PAUC_Accuracy(iter,itime)] = trainClassifier_svm_guassian(T);

%             T=array2table([Ptempdata Shuffle(AccuracyMask)']);
%             [~, ~,low_PAUC_Accuracy_sh(iter,itime)] = trainClassifier_svm_guassian(T);
%             
            % response
            T=array2table([Rtempdata HighMemMask']);
            [~, ~,low_RAUC_Memorability(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_LDA(T); %
            
            T=array2table([Rtempdata Shuffle(HighMemMask)']);
            [~, ~,low_RAUC_Memorability_sh(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_LDA(T); %
            
            T=array2table([Rtempdata (AccuracyMask)']);
            [~,~,low_RAUC_Accuracy(iter,itime)] = trainClassifier_svm_guassian(T);

%             T=array2table([Rtempdata Shuffle(AccuracyMask)']);
%             [~,~,low_RAUC_Accuracy_sh(iter,itime)] = trainClassifier_svm_guassian(T);
            
            %%===high freq
            Ptempdata=squeeze(High_Probe_DataToDecode(:,:,itime));
            Rtempdata=squeeze(High_Recall_DataToDecode(:,:,itime));
            
            % probe
            T=array2table([Ptempdata HighMemMask']);
            [~, ~,high_PAUC_Memorability(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_LDA(T); %

            T=array2table([Ptempdata Shuffle(HighMemMask)']);
            [~, ~,high_PAUC_Memorability_sh(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_LDA(T); %
            
            T=array2table([Ptempdata (AccuracyMask)']);
            [~, ~,high_PAUC_Accuracy(iter,itime)] = trainClassifier_svm_guassian(T);
            
%             T=array2table([Ptempdata Shuffle(AccuracyMask)']);
%             [~, ~,high_PAUC_Accuracy_sh(iter,itime)] = trainClassifier_svm_guassian(T);
                        
            % response
            T=array2table([Rtempdata HighMemMask']);
            [~, ~,high_RAUC_Memorability(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_LDA(T); %

            T=array2table([Ptempdata Shuffle(HighMemMask)']);
            [~, ~,high_RAUC_Memorability_sh(iter,itime)] = trainClassifier_svm_guassian(T); %trainClassifier_LDA(T); %
            
            T=array2table([Rtempdata (AccuracyMask)']);
            [~, ~,high_RAUC_Accuracy(iter,itime)] = trainClassifier_svm_guassian(T);            
 
            T=array2table([Rtempdata Shuffle(AccuracyMask)']);
            [~, ~,high_RAUC_Accuracy_sh(iter,itime)] = trainClassifier_svm_guassian(T);              
        end
    end
    
    GroupResults.SubjResults_PAUC_Memorability(isub,1,:)=mean(low_PAUC_Memorability);
    GroupResults.SubjResults_PAUC_Accuracy(isub,1,:)=mean(low_PAUC_Accuracy);
    GroupResults.SubjResults_RAUC_Memorability(isub,1,:)=mean(low_RAUC_Memorability);
    GroupResults.SubjResults_RAUC_Accuracy(isub,1,:)=mean(low_RAUC_Accuracy);
    
    GroupResults.SubjResults_PAUC_Memorability(isub,2,:)=mean(high_PAUC_Memorability);
    GroupResults.SubjResults_PAUC_Accuracy(isub,2,:)=mean(high_PAUC_Accuracy);
    GroupResults.SubjResults_RAUC_Memorability(isub,2,:)=mean(high_RAUC_Memorability);
    GroupResults.SubjResults_RAUC_Accuracy(isub,2,:)=mean(high_RAUC_Accuracy);

%     GroupResults.SubjResults_PAUC_Memorability_sh(isub,1,:)=mean(low_PAUC_Memorability_sh);
%     GroupResults.SubjResults_PAUC_Accuracy_sh(isub,1,:)=mean(low_PAUC_Accuracy_sh);
%     GroupResults.SubjResults_RAUC_Memorability_sh(isub,1,:)=mean(low_RAUC_Memorability_sh);
%     GroupResults.SubjResults_RAUC_Accuracy_sh(isub,1,:)=mean(low_RAUC_Accuracy_sh);
%     
%     GroupResults.SubjResults_PAUC_Memorability_sh(isub,2,:)=mean(high_PAUC_Memorability_sh);
%     GroupResults.SubjResults_PAUC_Accuracy_sh(isub,2,:)=mean(high_PAUC_Accuracy_sh);
%     GroupResults.SubjResults_RAUC_Memorability_sh(isub,2,:)=mean(high_RAUC_Memorability_sh);
%     GroupResults.SubjResults_RAUC_Accuracy_sh(isub,2,:)=mean(high_RAUC_Accuracy_sh);   

    figure(1000);clf
    subplot(2,2,1);plot(BinTime(tasktime),mean(low_PAUC_Memorability));
%     hold on; plot(BinTime(tasktime),mean(low_PAUC_Memorability_sh));
    hold on; plot(BinTime(tasktime),mean(low_PAUC_Accuracy)); 
%     hold on; plot(BinTime(tasktime),mean(low_PAUC_Accuracy_sh));     
    hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');
    title('Low frequency Probe Onset');legend('memorability','memSh','memory accuracy','acc_sh','location','southoutside');
    
    subplot(2,2,2);plot(BinTime(tasktime),mean(high_PAUC_Memorability));
    hold on; plot(BinTime(tasktime),mean(high_PAUC_Accuracy)); 
    hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');
    title('High frequency Probe Onset');legend('memorability','memory accuracy','location','southoutside');    
    
    subplot(2,2,3);plot(BinTime(tasktime),mean(low_RAUC_Memorability));
    hold on; plot(BinTime(tasktime),mean(low_RAUC_Accuracy)); 
    hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   
    title('Low frequency Probe Onset');legend('memorability','memory accuracy','location','southoutside');
    
    subplot(2,2,4);plot(BinTime(tasktime),mean(high_RAUC_Memorability));
    hold on; plot(BinTime(tasktime),mean(high_RAUC_Accuracy)); 
    hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   
    title('High frequency Probe Onset');legend('memorability','memory accuracy','location','southoutside');    
    
    h=gcf;
    set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
    printfilename=fullfile(outpath,['Decode_Sub' num2str(isub) '_SVM_allfreatures.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
%     savefilename=fullfile(outpath,['Resulst_Sub' num2str(isub) '.mat']);
%     save (savefilename, 'SubjResults', '-v7.3')
    
end


GroupResults
    savefilename=fullfile(outpath,['GroupResults'  '.mat']);
    save (savefilename, 'GroupResults', '-v7.3')


%%

excludeID= [46]; % low accuracy  38 53; low trial 62; low electrode 46;
allsubID=[WordtoAnalyzeTReduced.subID];
allsubID=allsubID;
IncludesubID=allsub(~ismember(allsub,excludeID));

IncludesubID=allsubID(ismember(allsubID,IncludesubID))

a = squeeze(GroupResults.SubjResults_PAUC_Memorability(IncludesubID,1,:));
a2 = squeeze(GroupResults.SubjResults_PAUC_Accuracy(IncludesubID,1,:));
a3 = squeeze(GroupResults.SubjResults_RAUC_Memorability(IncludesubID,1,:));
a4 = squeeze(GroupResults.SubjResults_RAUC_Accuracy(IncludesubID,1,:));


addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/waveform_comparisons_matlab'));
[acorr,kmin]=gutautocorr(a);
[minlen]=gsgutsims_paried_twotailed(size(a,1),size(a,2),0.10,acorr,2000,0.05);


for i=1:2000
    meana(i,:)=mean(a(randsample(1:18,18,'true'),:),1);
end

sigtest=(1-mean(meana>0.5)<0.025);



% a = squeeze(GroupResults.SubjResults_PRSA_Memorability4(IncludesubID,1,:));
% a2 = squeeze(GroupResults.SubjResults_PRSA_Accuracy(IncludesubID,1,:));
% a3 = squeeze(GroupResults.SubjResults_RRSA_Memorability4(IncludesubID,1,:));
% a4 = squeeze(GroupResults.SubjResults_PRSA_Accuracy(IncludesubID,1,:));


figure;
subplot(2,2,1);
prettyplot(BinTime(tasktime),mean(a),std(a)/sqrt(length(IncludesubID)),[1 0 0],1);
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');
hold on;plot(BinTime(tasktime),ttest(a,0.5,'tail','right')*0.6)
xtime=BinTime(tasktime);
[B, N, BI]=RunLength(sigtest');
plotix=(BI(B==1 & N>minlen));
for iplox=1:length(plotix)
    if xtime(plotix(iplox))>0
        hold on;plot(xtime(plotix(iplox):plotix(iplox)+N(BI==plotix(iplox))-1),ones(N(BI==plotix(iplox)),1)*0.5,'ko');
    end
end
axis([-1500 1500 0.45 0.65])
title('probe onset')

subplot(2,2,2);
prettyplot(BinTime(tasktime),mean(a2),std(a2)/sqrt(length(IncludesubID)),[0 1 0],1);
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   
hold on;plot(BinTime(tasktime),ttest(a2,0.5,'tail','right')*0.6)

for i=1:2000
    meana(i,:)=mean(a2(randsample(1:18,18,'true'),:),1);
end
sigtest=(1-mean(meana>0.5)<0.025);
xtime=BinTime(tasktime);
[B, N, BI]=RunLength(sigtest');
plotix=(BI(B==1 & N>minlen));
for iplox=1:length(plotix)
    if xtime(plotix(iplox))>0
        hold on;plot(xtime(plotix(iplox):plotix(iplox)+N(BI==plotix(iplox))-1),ones(N(BI==plotix(iplox)),1)*0.5,'ko');
    end
end

axis([-1500 1500 0.45 0.63])
title('probe onset')

subplot(2,2,3);
prettyplot(BinTime(tasktime),mean(a3),std(a3)/sqrt(length(IncludesubID)),[1 0 0],1);
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');
hold on;plot(BinTime(tasktime),ttest(a3,0.5,'tail','right')*0.6)
axis([-1500 1500 0.45 0.63])
title('recall onset')

subplot(2,2,4);
prettyplot(BinTime(tasktime),mean(a4),std(a4)/sqrt(length(IncludesubID)),[0 1 0],1);
hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');   
hold on;plot(BinTime(tasktime),ttest(a4,0.5,'tail','right')*0.6)
axis([-1500 1500 0.45 0.63])
% hold on;plot(BinTime(tasktime),0.5*ttest(a,1))
title('recall onset')



%% 
for isubtem=1:size(WordtoAnalyzeTReduced,1)
    
isub =  WordtoAnalyzeTReduced.subID(isubtem)
all_high_PAUC_Memorability(isubtem,:)=  SubjResults(isub).SubjResults_high_PAUC_Memorability(isub,:);
all_low_PAUC_Memorability(isubtem,:)=  SubjResults(isub).SubjResults_low_PAUC_Memorability(isub,:);

all_high_PAUC_Accuracy(isubtem,:)=  SubjResults(isub).SubjResults_high_PAUC_Accuracy(isub,:);
all_low_PAUC_Accuracy(isubtem,:)=  SubjResults(isub).SubjResults_low_PAUC_Accuracy(isub,:);

end

excludeID=[29 34 37 38 53 57 ];
allsubID=[WordtoAnalyzeTReduced.subID];
IncludesubID=(~ismember(allsubID,excludeID));
IncludesubID=[WordtoAnalyzeTReduced.meanACC{1:end}]>0.1;

figure;
subplot(1,2,1);
prettyplot(BinTime(tasktime),mean(all_low_PAUC_Memorability(IncludesubID,:)),std(all_low_PAUC_Memorability(IncludesubID,:))/sqrt(sum(IncludesubID)),[1 0 0],1);
hold on;plot(BinTime(tasktime),0.4*ttest(all_low_PAUC_Memorability(IncludesubID,:),0.5))

prettyplot(BinTime(tasktime),mean(all_low_PAUC_Memorability(IncludesubID,:)));
hold on;plot(BinTime(tasktime),mean(all_low_PAUC_Accuracy(IncludesubID,:)));
    hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');

subplot(1,2,2);plot(BinTime(tasktime),mean(all_high_PAUC_Memorability(IncludesubID,:)));
hold on;plot(BinTime(tasktime),mean(all_high_PAUC_Accuracy(IncludesubID,:)));
    hold on; plot([min(BinTime(tasktime)) max(BinTime(tasktime))] ,[0.5 0.5],'k-');


%% calculate the average experimental effects.

% excludeID=[33 34 37];
%
% excludeID=[33 43 53 55 57 62 ];
excludeID=[];
allsubID=[WordtoAnalyzeTReduced.subID];
IncludesubID=allsubID(~ismember(allsubID,excludeID));

InteractionPower=[];MainAccPower=[];MainMemorabilityPower=[];
for isub=1:length(IncludesubID)
    
    InteractionPower(isub,:,:)=SubjResults(IncludesubID(isub)).Probe_InteractionEff;
    MainAccPower(isub,:,:)=SubjResults(IncludesubID(isub)).Probe_AccuracyMainEffect;
    MainMemorabilityPower(isub,:,:)=SubjResults(IncludesubID(isub)).Probe_MemorabiltyMainEffect;
    
end

timePercent

figure;
subplot(2,3,1);
y=squeeze(nanmean(MainAccPower))
hold on;contourf(timePercent(tasktime),waveletFreqs,y,200,'linecolor','none');
set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
colorbar;caxis([-0.08 0.08])
title('Main Effect of Memory Accuracy')
xlabel('% of retrieval time');ylabel('freq');

subplot(2,3,4);
y=squeeze(nanmean(MainAccPower))
contourf(timePercent(tasktime),waveletFreqs,y,200,'linecolor','none');
set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
colorbar;caxis([-0.08 0.08])
title('Main Effect of Memory Accuracy')
xlabel('% of retrieval time');ylabel('freq');
hold on;contour(timePercent(tasktime),waveletFreqs,squeeze(ttest(MainAccPower)))


subplot(2,3,2);
y=squeeze(nanmean(MainMemorabilityPower))
contourf(timePercent(tasktime),waveletFreqs,y,200,'linecolor','none');
set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
colorbar;caxis([-0.08 0.08])
title('Main Effect of Target Retrievablity')
xlabel('% of retrieval time');ylabel('freq');

subplot(2,3,5);
y=squeeze(nanmean(MainMemorabilityPower))
contourf(timePercent(tasktime),waveletFreqs,y,200,'linecolor','none');
set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
colorbar;caxis([-0.08 0.08])
title('Main Effect of Target Retrievablity')
xlabel('% of retrieval time');ylabel('freq');
hold on;contour(timePercent(tasktime),waveletFreqs,squeeze(ttest(MainMemorabilityPower)))


subplot(2,3,3);
y=squeeze(nanmean(InteractionPower))
contourf(timePercent(tasktime),waveletFreqs,y,200,'linecolor','none');
set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
colorbar;caxis([-0.08 0.08])
title('Interaction Higher-Lower Memorablilty (Accurate-Inaccurate)')
xlabel('% of retrieval time');ylabel('freq');

subplot(2,3,6);
y=squeeze(nanmean(InteractionPower))
contourf(timePercent(tasktime),waveletFreqs,y,200,'linecolor','none');
set(gca, 'yscale', 'log', 'ytick', 2.^[1:8]);
colorbar;caxis([-0.08 0.08])
title('Interaction Higher-Lower Memorablilty (Accurate-Inaccurate)')
xlabel('% of retrieval time');ylabel('freq');
hold on;contour(timePercent(tasktime),waveletFreqs,-squeeze(ttest(InteractionPower)))

h=gcf;
set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
printfilename=fullfile(outpath,['Probe_Group_AveragePower.pdf']);
print(h,printfilename,'-dpdf','-bestfit');


%%
% duration variable across subjects
% nanmean(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>50)) + 500 + 1000
duration = 3000;
offset   = -1500;
buffer   = 1000;
resamp   = 500;
time  =  downsample(linspace(offset,duration+offset,duration),1000/resamp);

% moving window approach, no overlaps, with 20 ms gap.
points_per_win=round(20/(1000/resamp));
points_per_slide=round(20/(1000/resamp)); % data has been downsampled!
BinTime=windowed_average(time,points_per_win,points_per_slide);
%     SubjResults.percentTimeInMS=20;

tasktime = BinTime>=-1250 & BinTime<=1250;

for isubtem=1:size(WordtoAnalyzeTReduced,1)
    isubtem
    
    clear Probe_WaveFreqTime_Z Recall_WaveFreqTime_Z ProbePowerlogAllSessionNormtime HighMemMask AccuracyMask RecallPowerlogAllSessionNormtime
    
    isub =  WordtoAnalyzeTReduced.subID(isubtem)
    
    ProbePowerlogAllSessionNormtime= SubjResults.ProbePowerlogAllSessionNormtime;
    RecallPowerlogAllSessionNormtime= SubjResults.RecallPowerlogAllSessionNormtime;
    
    %     HighMemMask =  SubjResults.HighMemMask;
    
    MemorabilityScore =  WordtoAnalyzeTReduced.ResponseMemorability{isubtem};
    MemorabilityScore = MemorabilityScore(SubjResults.allcleanTrials==1);
    
    HighMemMask = MemorabilityScore >= medianMemor;
    %     MemorabilityScore = MemorabilityScore+0.005;
    
    AccuracyMask =  SubjResults.AccuracyMask;
    
    clear Probe_DataToDecode Recall_DataToDecode
    freqBandYticksSplit=[exp(linspace(log(2),log(8),4));exp(linspace(log(8),log(32),4));exp(linspace(log(32),log(150),4))]; % only bin the data into 3 bins.
    
    freqBandYticksSplit=[2  4  8  12  32  70  150]; % only bin the data into 3 bins.
    
    electrodeToAna=ismember(SubjResults.allcleanChans,WordtoAnalyzeTReduced.ATLchan{isubtem});
    
    for ifreqs=1
        ifreqs
        freqBandYticks=freqBandYticksSplit(ifreqs,:) % only bin the data into 3 bins.
        clear Probe_WaveFreqTime_Z Recall_WaveFreqTime_Z Probe_DataToDecode Recall_DataToDecode
        for ifreq=1:length(freqBandYticks)-1
            Probe_WaveFreqTime_Z(:,:,ifreq,:) = single(squeeze(mean(ProbePowerlogAllSessionNormtime(electrodeToAna,:,waveletFreqs>=freqBandYticks(ifreq) & waveletFreqs<freqBandYticks(ifreq+1),tasktime),3)));
            Recall_WaveFreqTime_Z(:,:,ifreq,:) = single(squeeze(mean(RecallPowerlogAllSessionNormtime(electrodeToAna,:,waveletFreqs>=freqBandYticks(ifreq) & waveletFreqs<freqBandYticks(ifreq+1),tasktime),3)));
        end
        
        for itrl=1:size(Recall_WaveFreqTime_Z,2)
            for itime=1:size(Recall_WaveFreqTime_Z,4)
                Probe_DataToDecode(ifreqs,itrl,:,itime)=reshape(Probe_WaveFreqTime_Z(electrodeToAna,itrl,:,itime),sum(electrodeToAna)*size(Probe_WaveFreqTime_Z,3),1);
                Recall_DataToDecode(ifreqs,itrl,:,itime)=reshape(Recall_WaveFreqTime_Z(electrodeToAna,itrl,:,itime),sum(electrodeToAna)*size(Recall_WaveFreqTime_Z,3),1);
            end
        end
    end
    
    
    % MeanFreqBinPower=squeeze(mean(ProbePowerlogAllSession(:,:,waveletFreqs>=8 & waveletFreqs<=16,:),3));
    
    %%
    nter=10;
    
    for ifreq=1
        clear PAUC_Memorability PAUC_Memorability_sh PAUC_Accuracy PAUC_Accuracy_sh RAUC_Memorability RAUC_Memorability_sh RAUC_Accuracy RAUC_Accuracy_sh
        
        for iter=1:nter
            parfor itime=1:size(Recall_DataToDecode,4)
                itime
                Ptempdata=squeeze(Probe_DataToDecode(1,:,:,itime));
                Rtempdata=squeeze(Recall_DataToDecode(1,:,:,itime));
                variance=95;
                
                %                 Ptempdata=squeeze(Probe_WaveFreqTime_Z(:,:,ifreq,itime))';
                %                 Rtempdata=squeeze(Recall_WaveFreqTime_Z(:,:,ifreq,itime))';
                
                T=array2table([Ptempdata HighMemMask']);
                [~, PMemorability_acc(iter,itime),PAUC_Memorability(iter,itime)] = trainClassifier_LDA_pca(T,variance); %trainClassifier_LDA(T); %
                
                T=array2table([Ptempdata Shuffle(HighMemMask)']);
                [~, PMemorability_sh_acc(iter,itime),PAUC_Memorability_sh(iter,itime)] = trainClassifier_LDA_pca(T,variance);
                
                T=array2table([Ptempdata (AccuracyMask)']);
                [~, PAccuracy_acc(iter,itime),PAUC_Accuracy(iter,itime)] = trainClassifier_LDA_pca(T,variance);
                
                T=array2table([Ptempdata Shuffle(AccuracyMask)']);
                [~, PAccuracy_sh_acc(iter,itime),PAUC_Accuracy_sh(iter,itime)] = trainClassifier_LDA_pca(T,variance);
                
                
                T=array2table([Rtempdata HighMemMask']);
                [~, RMemorability_acc(iter,itime),RAUC_Memorability(iter,itime)] = trainClassifier_LDA_pca(T,variance); %trainClassifier_LDA(T); %
                
                T=array2table([Rtempdata Shuffle(HighMemMask)']);
                [~, RMemorability_sh_acc(iter,itime),RAUC_Memorability_sh(iter,itime)] = trainClassifier_LDA_pca(T,variance);
                
                T=array2table([Rtempdata (AccuracyMask)']);
                [~, RAccuracy_acc(iter,itime),RAUC_Accuracy(iter,itime)] = trainClassifier_LDA_pca(T,variance);
                
                T=array2table([Rtempdata Shuffle(AccuracyMask)']);
                [~, RAccuracy_sh_acc(iter,itime),RAUC_Accuracy_sh(iter,itime)] = trainClassifier_LDA_pca(T,variance);
            end
        end
        %
        SubjDecoding(isub).PAUC_Memorability(ifreq,:,:)=mean(PAUC_Memorability);
        SubjDecoding(isub).PAUC_Memorability_sh(ifreq,:,:)=mean(PAUC_Memorability_sh);
        SubjDecoding(isub).PAUC_Accuracy(ifreq,:,:)=mean(PAUC_Accuracy);
        SubjDecoding(isub).PAUC_Accuracy_sh(ifreq,:,:)=mean(PAUC_Accuracy_sh);
        SubjDecoding(isub).RAUC_Memorability(ifreq,:,:)=mean(RAUC_Memorability);
        SubjDecoding(isub).RAUC_Memorability_sh(ifreq,:,:)=mean(RAUC_Memorability_sh);
        SubjDecoding(isub).RAUC_Accuracy(ifreq,:,:)=mean(RAUC_Accuracy);
        SubjDecoding(isub).RAUC_Accuracy_sh(ifreq,:,:)=mean(RAUC_Accuracy_sh);
        
        h=figure(100+ifreq);clf
        set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
        subplot(2,2,1);
        plot(BinTime(tasktime),mean(PAUC_Memorability)); hold on;
        plot(BinTime(tasktime),mean(PAUC_Memorability_sh)); hold on;
        axis([min(BinTime(tasktime)) max(BinTime(tasktime)) 0.3 0.8])
        legend('memorability','shuffleMem');title('Decoding AUC timelocked to Probe onset')
        
        subplot(2,2,2);
        plot(BinTime(tasktime),mean(PAUC_Accuracy)); hold on;
        plot(BinTime(tasktime),mean(PAUC_Accuracy_sh)); hold on;
        axis([min(BinTime(tasktime)) max(BinTime(tasktime)) 0.3 0.8])
        legend('recall accuracy','shuffleACC');title('Decoding AUC timelocked to Probe onset')
        
        subplot(2,2,3);
        plot(BinTime(tasktime),mean(RAUC_Memorability)); hold on;
        plot(BinTime(tasktime),mean(RAUC_Memorability_sh)); hold on;
        axis([min(BinTime(tasktime)) max(BinTime(tasktime)) 0.3 0.8])
        legend('memorability','shuffleMem');title('Decoding AUC timelocked to recall onset')
        
        subplot(2,2,4);
        plot(BinTime(tasktime),mean(RAUC_Accuracy)); hold on;
        plot(BinTime(tasktime),mean(RAUC_Accuracy_sh)); hold on;
        axis([min(BinTime(tasktime)) max(BinTime(tasktime)) 0.3 0.8])
        legend('recall accuracy','shuffleACC');title('Decoding AUC timelocked to recall onset')
        
        printfilename=fullfile(outpath,['Subj_' num2str(isub) '_' 'FreqBIN' num2str(ifreq) '.pdf']);
        %         print(h,printfilename,'-dpdf','-bestfit');
        
        
    end
end


%%

save SubjDecoding_ATL_2ROIs.mat SubjDecoding



%%
%
for isubtem=1:size(WordtoAnalyzeTReduced,1)
    
    %%
    
    
    clear Probe_WaveFreqTime_Z Recall_WaveFreqTime_Z ProbePowerlogAllSessionNormtime HighMemMask AccuracyMask RecallPowerlogAllSessionNormtime
    
    isub =  WordtoAnalyzeTReduced.subID(isubtem)
    
    ProbePowerlogAllSessionNormtime= SubjResults.ProbePowerlogAllSessionNormtime;
    RecallPowerlogAllSessionNormtime= SubjResults.RecallPowerlogAllSessionNormtime;
    
    %     HighMemMask =  SubjResults.HighMemMask;
    
    MemorabilityScore =  WordtoAnalyzeTReduced.ResponseMemorability{isubtem};
    MemorabilityScore = MemorabilityScore(SubjResults.allcleanTrials==1);
    
    HighMemMask = MemorabilityScore > medianMemor;
    %     MemorabilityScore = MemorabilityScore+0.005;
    
    AccuracyMask =  SubjResults.AccuracyMask;
    
    clear Probe_DataToDecode Recall_DataToDecode
    freqBandYticksSplit=[exp(linspace(log(2),log(8),4));exp(linspace(log(8),log(32),4));exp(linspace(log(32),log(150),4))]; % only bin the data into 3 bins.
    
    freqBandYticksSplit=[2  4  8  12  32  70  150]; % only bin the data into 3 bins.
    
    electrodeToAna=ismember(SubjResults.allcleanChans,WordtoAnalyzeTReduced.ATLchan{isubtem});
    
    for ifreqs=1
        ifreqs
        freqBandYticks=freqBandYticksSplit(ifreqs,:) % only bin the data into 3 bins.
        clear Probe_WaveFreqTime_Z Recall_WaveFreqTime_Z Probe_DataToDecode Recall_DataToDecode
        for ifreq=1:length(freqBandYticks)-1
            Probe_WaveFreqTime_Z(:,:,ifreq,:) = single(squeeze(mean(ProbePowerlogAllSessionNormtime(electrodeToAna,:,waveletFreqs>=freqBandYticks(ifreq) & waveletFreqs<freqBandYticks(ifreq+1),tasktime),3)));
            Recall_WaveFreqTime_Z(:,:,ifreq,:) = single(squeeze(mean(RecallPowerlogAllSessionNormtime(electrodeToAna,:,waveletFreqs>=freqBandYticks(ifreq) & waveletFreqs<freqBandYticks(ifreq+1),tasktime),3)));
        end
        
        for itrl=1:size(Recall_WaveFreqTime_Z,2)
            for itime=1:size(Recall_WaveFreqTime_Z,4)
                Probe_DataToDecode(ifreqs,itrl,:,itime)=reshape(Probe_WaveFreqTime_Z(electrodeToAna,itrl,:,itime),sum(electrodeToAna)*size(Probe_WaveFreqTime_Z,3),1);
                Recall_DataToDecode(ifreqs,itrl,:,itime)=reshape(Recall_WaveFreqTime_Z(electrodeToAna,itrl,:,itime),sum(electrodeToAna)*size(Recall_WaveFreqTime_Z,3),1);
            end
        end
    end
    % RSA based within trial across electrodes
    clear P_Memo_RSAtime  P_ACC_RSAtime  R_Memo_RSAtime  R_ACC_RSAtime  MemoModel
    
    
    
    for ifreq=1:6
        %             for ic=1:size(Probe_WaveFreqTime_Z,1)
        for itime=1:125
            %                 Ptempdata=squeeze(Probe_DataToDecode(:,:,:,itime));
            %                 Rtempdata=squeeze(Recall_DataToDecode(:,:,:,itime));
            Ptempdata=squeeze(Probe_WaveFreqTime_Z(:,:,ifreq,itime))';
            Rtempdata=squeeze(Recall_WaveFreqTime_Z(:,:,ifreq,itime))';
            %
            %                 MemoModel=1-HighMemMask.*(HighMemMask)';
            AccuracyModel=1-AccuracyMask.*AccuracyMask';
            
            %                 for i=1:length(MemorabilityScore)
            %                         MemoModel(i,:)=MemorabilityScore-MemorabilityScore(i);
            %                 end
            
            for i=1:length(HighMemMask)
                for j=1:length(HighMemMask)
                    MemoModel(i,j)=abs(MemorabilityScore(i)-MemorabilityScore(j));
                    %                         if HighMemMask(i)==HighMemMask(j) % low and high memory are highly similiar within each group
                    %                          MemoModel(i,j)=0;
                    %                         else
                    %                          MemoModel(i,j)=1;
                    %                         end
                end
            end
            
            RSA=(squareform(pdist(Ptempdata,'cosine')));
            a1 = RSA(find(tril(ones(size(RSA)),-1)));
            
            clear RSA
            RSA=(squareform(pdist(Rtempdata,'cosine')));
            a2 = RSA(find(tril(ones(size(RSA)),-1)));
            
            b = MemoModel(find(tril(ones(size(RSA)),-1)));
            c = AccuracyModel(find(tril(ones(size(RSA)),-1)));
            
            P_Memo_RSAtime(itime)=corr(a1,b,'type','Spearman');
            P_ACC_RSAtime(itime)=corr(a1,c,'type','Spearman');
            
            R_Memo_RSAtime(itime)=corr(a2,b,'type','Spearman');
            R_ACC_RSAtime(itime)=corr(a2,c,'type','Spearman');
            
        end
        
        
        all_P_Memo_RSAtime(isub,ifreq,:)=P_Memo_RSAtime;
        all_P_ACC_RSAtime(isub,ifreq,:)=P_ACC_RSAtime;
        all_R_Memo_RSAtime(isub,ifreq,:)=R_Memo_RSAtime;
        all_R_ACC_RSAtime(isub,ifreq,:)=R_ACC_RSAtime;
        
        figure(isub);
        subplot(2,6,ifreq);plot(BinTime(tasktime),[P_Memo_RSAtime;P_ACC_RSAtime]');title('probe');legend('memorability','accuracy');
        subplot(2,6,ifreq+6);plot(BinTime(tasktime),[R_Memo_RSAtime;R_ACC_RSAtime]');title('recall');legend('memorability','accuracy');
    end
end
%
%
%
%             all_P_Memo_RSAtime(isub,ifreq,:)=P_Memo_RSAtime;
%             all_P_ACC_RSAtime(isub,ifreq,:)=P_ACC_RSAtime;
%             all_R_Memo_RSAtime(isub,ifreq,:)=R_Memo_RSAtime;
%             all_R_ACC_RSAtime(isub,ifreq,:)=R_ACC_RSAtime;


%%

excludeID=[35 53]; % 28 37 40 33 43 29 38 62 57
allsubID=[WordtoAnalyzeTReduced.subID];
IncludesubID=allsubID(~ismember(allsubID,excludeID));


a=0.5.*log((1+all_P_Memo_RSAtime)./(1-all_P_Memo_RSAtime));

a = squeeze(a(IncludesubID,1,:))
figure;
prettyplot(BinTime(tasktime),mean(a),std(a)/sqrt(size(a,1)),[0 0 1],0.3)



%%









%% average decoding results
excludeID=[29 33 34 37 38 39 46 53 57];
% excludeID=[33 43 53 55 57 62 ];
excludeID=[53]; % 33 43 46 29 38 53
allsubID=[WordtoAnalyzeTReduced.subID];
IncludesubID=allsubID(~ismember(allsubID,excludeID));


Decoding_ACC=[];Decoding_Memorability=[];Decoding_ShuffleMemo=[];
RDecoding_ACC=[];RDecoding_Memorability=[];RDecoding_ShuffleMemo=[];
PDecoding_ACC=[];
PDecoding_Memorability=[];
PDecoding_ShuffleMemo=[];

for ifreq=1:6
    
    for isub=1:length(IncludesubID)
        
        
        PDecoding_ACC(isub,ifreq,:)=SubjDecoding(IncludesubID(isub)).PAUC_Accuracy(ifreq,:);
        PDecoding_Memorability(isub,ifreq,:)=SubjDecoding(IncludesubID(isub)).PAUC_Memorability(ifreq,:);
        PDecoding_ShuffleMemo(isub,ifreq,:)=SubjDecoding(IncludesubID(isub)).PAUC_Memorability_sh(ifreq,:);
        
        RDecoding_ACC(isub,ifreq,:)=SubjDecoding(IncludesubID(isub)).RAUC_Accuracy(ifreq,:);
        RDecoding_Memorability(isub,ifreq,:)=SubjDecoding(IncludesubID(isub)).RAUC_Memorability(ifreq,:);
        RDecoding_ShuffleMemo(isub,ifreq,:)=SubjDecoding(IncludesubID(isub)).RAUC_Memorability_sh(ifreq,:);
    end
    
end


figure;
for ifreq=1:6
    a=  squeeze(PDecoding_ACC(:,ifreq,:));
    b=  squeeze(PDecoding_Memorability(:,ifreq,:));
    c=  squeeze(PDecoding_ShuffleMemo(:,ifreq,:));
    
    
    subplot(2,3,ifreq)
    prettyplot(BinTime(tasktime),mean(a),std(a)/sqrt(size(a,1)),[1 0 0],1);
    hold on;prettyplot(BinTime(tasktime),mean(b),std(b)/sqrt(size(a,1)),[0 1 0],1);
    hold on;prettyplot(BinTime(tasktime),mean(c),std(c)/sqrt(size(a,1)),[0 0 1],1);
    legend('accuracy', '', 'memorability','', 'shuffle', '','location','northeastoutside')
    axis([-1000 1000 0.4 0.6])
    
end

figure;
for ifreq=1:6
    a=  squeeze(RDecoding_ACC(:,ifreq,:));
    b=  squeeze(RDecoding_Memorability(:,ifreq,:));
    c=  squeeze(RDecoding_ShuffleMemo(:,ifreq,:));
    
    
    subplot(2,3,ifreq)
    prettyplot(BinTime(tasktime),mean(a),std(a)/sqrt(size(a,1)),[1 0 0],1);
    hold on;prettyplot(BinTime(tasktime),mean(b),std(b)/sqrt(size(a,1)),[0 1 0],1);
    hold on;prettyplot(BinTime(tasktime),mean(c),std(c)/sqrt(size(a,1)),[0 0 1],1);
    legend('accuracy', '', 'memorability','', 'shuffle', '','location','northeastoutside')
    axis([-1000 1000 0.4 0.6])
    
end

% print(gcf,'scaled_decodingacc','-dpdf','-bestfit');
%% PLOT ALL Electrodes
close all

bp = brainplotter()
bd = braindata2()
bd.loadAverage()
bp =bd.ezplot(bp)
bp.view('ventral')

datapath='/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/Scripts/3_Mono/data/ATL_05_19';

for isub = allsub
    isub
     savefilename=fullfile(datapath,['SubjResults_ATL' num2str(isub) '.mat']);
    
    
    if exist(savefilename)>0
        
        load(savefilename);

        % identify and limit the electrodes for analysis, only look at the surface ALT electrodes.
        subj     = ['NIH0' num2str(isub)];
        
        % load subject braindata
        bd = braindata2(subj,rootEEGdir)
        
        % load each electrodes mesh info on the standardized brain
        t  = bd.roi.lead_mesh;
        
        clear elecs
        % get example electrodes
        elecs=SubjResults.allcleanChans;
        
        
        % Electrodes to exclude for misalignment
        % NIH026 G27
        % NIH029 G17
        % NIH032 ROF3
        % NIH036 OF3, OF4
        % NIH062 OF4
        % NIH066 TG84, TG123
        
        %         if strcmp(subj,'NIH026')
        %             elecs=elecs(~ismember(elecs,{'G27'}));
        %         elseif strcmp(subj,'NIH029')
        %             elecs=elecs(~ismember(elecs,{'G17'}));
        %         elseif strcmp(subj,'NIH032')
        %             elecs=elecs(~ismember(elecs,{'ROF3'}));
        %         elseif strcmp(subj,'NIH036')
        %             elecs=elecs(~ismember(elecs,{'OF3','OF4'}));
        %         elseif strcmp(subj,'NIH062')
        %             elecs=elecs(~ismember(elecs,{'OF4'}));
        %         elseif strcmp(subj,'NIH066')
        %             elecs=elecs(~ismember(elecs,{'TG84', 'TG123'}));
        %         end
        
        
        for j = 1:length(elecs)
            % get channels index within lead mesh table
            chan = strcmp(t.chanName,elecs{j});
            
            % get surface to plot on (left or right hemisphere)
            s = ['pial_' t.whichHemi{chan}];
            bp.plotPoint(t(chan,:),...
                'radius',1,...
                'mesh_clust_d',0,...
                'alpha',0.5,...
                'color',[1 0 0 ],...
                'surf',s) %'label', elecs{j}
            
            %   mesh_clust_d - If you are passing mesh index lists, by default a point is placed at the first mesh index in the list. However, for sulcus-
            %             spanning electrodes, we may want to plot 2+ points at each group of mesh indices. Pass a distance here if you want to clusterize
            %             such that mesh_indices >=~ mesh_clust_d apart will be plotted as separate points. Default is 3. 0 to just use the first index.
            % ,...
            
            
        end
        
        %keyboard
    end
    
end



% Electrodes to exclude for misalignment
% NIH026 G27
% NIH029 G17
% NIH032 ROF3
% NIH036 OF3, OF4
% NIH062 OF4
% NIH066 TG84, TG123


