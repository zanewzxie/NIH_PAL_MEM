clear all
%add toolboxes
addpath(genpath('/Volumes/Zane/Matlab/eeg_toolbox/trunk'));
addpath(genpath('/Volumes/Zane/Matlab/dungeon_toolbox_17a'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/EEG_Preprocessing'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/IEM_ester'));
addpath(genpath('/Volumes/Zane/NIH_NINDS/Data_InProgress/SpecFuN'));

%% extract PAL behaviroal data
allsub=[26:63 66];
rootEEGdir = '/Volumes/Zane/NIH_FRNU_ROOT';                      %office-local

load('SubjTable_palRAMword.mat')
load('PAL_Memo_PALRAM.mat')

ToAnalyzeSub=allsub;
% PAL_Memo.Combinedmemorability = zscore(PAL_Memo.Probememorability) + (PAL_Memo.Responsememorability - nanmean(PAL_Memo.Responsememorability))./nanstd(PAL_Memo.Responsememorability);
% PAL_Memo.Combinedmemorability = PAL_Memo.Combinedmemorability/2;

medianMemor = median(PAL_Memo.Responsememorability);

for isubtem=1:length(ToAnalyzeSub)
    isubtem
    isub=ToAnalyzeSub(isubtem);
    % for i=1:length(SubjTable(isub).PALTable.response)
    %     if strcmp(num2str(SubjTable(isub).PALTable.response{i}),'-999')| strcmp(SubjTable(isub).PALTable.response{i},'<>') | strcmp(SubjTable(isub).PALTable.response{i},'PASS')
    %         SubjTable(isub).PALTable.response{i} = 'NORESPONSEXX';
    %     end
    % end
     
    if ~isempty(SubjTable(isub).PALTable)
        WordtoAnalyze(isubtem).Word=SubjTable(isub).PALTable.expected(ismember(SubjTable(isub).PALTable.expected,PAL_Memo.uniquewordsID_10));
        WordtoAnalyze(isubtem).InDX=find(ismember(SubjTable(isub).PALTable.expected,PAL_Memo.uniquewordsID_10));
        WordtoAnalyze(isubtem).WordCount=length(WordtoAnalyze(isubtem).InDX);
        
        for iev=1:length(WordtoAnalyze(isubtem).InDX)
            
            %first find the memorability score for when the word being a probe and
            %being a response word
            WordtoAnalyze(isubtem).ProbeMemorability(iev)=PAL_Memo.Probememorability(ismember(PAL_Memo.uniquewordsID_10,WordtoAnalyze(isubtem).Word(iev)));
            WordtoAnalyze(isubtem).ResponseMemorability(iev)=PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,WordtoAnalyze(isubtem).Word(iev)));
            WordtoAnalyze(isubtem).ResponseMemorability_z(iev)=PAL_Memo.Responsememorability_z(ismember(PAL_Memo.uniquewordsID_10,WordtoAnalyze(isubtem).Word(iev)));
            
            %     WordtoAnalyze(isubtem).CombinedMemorability(iev)=PAL_Memo.Combinedmemorability(ismember(PAL_Memo.uniquewordsID_10,WordtoAnalyze(isubtem).Word(iev)));
        end
        
        memratio(isubtem)=sum(WordtoAnalyze(isubtem).ResponseMemorability>medianMemor)/sum(WordtoAnalyze(isubtem).ResponseMemorability>0);
        WordtoAnalyze(isubtem).memratio=memratio(isubtem);
        
        %% identify and limit the electrodes for analysis, only look at the surface ALT electrodes.
        subj     = ['NIH0' num2str(isub)];
        info     = getElementInfo(subj, rootEEGdir);
        jack     = getJackTable(subj, rootEEGdir);
        
        
        % got to subject's folder and identify ATL.
        atlasfilename=fullfile(rootEEGdir,[subj '/tal/atlas/atlas_monopolar_simple.csv']);
        Atlas    = readtable(atlasfilename);
        
        ATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & (strcmp(Atlas.label_desikan,'Inferior temporal gyrus')...
            | strcmp(Atlas.label_desikan,'Middle temporal gyrus')...
            | strcmp(Atlas.label_desikan,'Superior temporal gyrus'))); %
        
        superiorATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Superior temporal gyrus')); %
        middleATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Middle temporal gyrus')); %
        inferiorATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Inferior temporal gyrus')); %
        
        MTLchan     = Atlas.chanName(strcmp(Atlas.label_desikan,'Entorhinal cortex') | strcmp(Atlas.label_desikan,'Parahippocampal gyrus')); %
        Auditorychan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Posterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Superior temporal gyrus')); %
        
        chan=[ATLchan;MTLchan;Auditorychan];
        
        
        allwords=WordtoAnalyze(isubtem).Word;
        
        for iev=1:length(SubjTable(isub).alignedEvents.TEST_PROBE)
            SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile,'noreref','processed');
            SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile,'/Volumes/Shares/FRNU/data/eeg',rootEEGdir);
            
            SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile,'noreref','processed');
            SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile,'/Volumes/Shares/FRNU/data/eeg',rootEEGdir);
        end
        
        WordtoAnalyze(isubtem).SubjTable=SubjTable(isub);
        % here exclude bad channels based on Julio's script.
        uniquesession=unique({SubjTable(isub).alignedEvents.TEST_PROBE.eegfile});
        allssession=[SubjTable(isub).alignedEvents.TEST_PROBE.session];
        
        allbadchan=[];availablechan=[];
        for iss=1:length(uniquesession)
            badchans{iss}=load(fullfile(uniquesession{iss},'bad_chans.mat'))
            allbadchan=[allbadchan;badchans{iss}.bad_chans];
            
            % to the session folder and identify the channels are only sepecific to
            % the avaialble session.
            temptablesfilename =fullfile(uniquesession{iss},'variance.csv');
            temptables    = readtable(temptablesfilename);
            availablechan=[availablechan; temptables.chanName];
        end
        uniqueallchan=unique(availablechan);
        if ~isempty(allbadchan)
        chan=chan(~ismember(chan,allbadchan));
        end
        chan=chan(ismember(chan,uniqueallchan));
        
        % get rid of channels that are at the boundary of ROI and difficult
        % to normalized or being shifted
        if strcmp(subj,'NIH026')
            chan=chan(~ismember(chan,{'G27'}));
        elseif strcmp(subj,'NIH029')
            chan=chan(~ismember(chan,{'G17'}));
        elseif strcmp(subj,'NIH032')
            chan=chan(~ismember(chan,{'ROF3'}));
        elseif strcmp(subj,'NIH036')
            chan=chan(~ismember(chan,{'OF3','OF4','TT_sh6','TT6','AST_sh4','AST4'}));
        elseif strcmp(subj,'NIH062')
            chan=chan(~ismember(chan,{'OF4'}));
        elseif strcmp(subj,'NIH066')
            chan=chan(~ismember(chan,{'TG84', 'TG123'}));
        end
        
        ATLchan=ATLchan(ismember(ATLchan,chan));
        
        superiorATLchan=superiorATLchan(ismember(superiorATLchan,chan));
        middleATLchan=middleATLchan(ismember(middleATLchan,chan));
        inferiorATLchan=inferiorATLchan(ismember(inferiorATLchan,chan));
        Auditorychan=Auditorychan(ismember(Auditorychan,chan));
        MTLchan=MTLchan(ismember(MTLchan,chan));
        
        clear chanHemisphere
        if length(chan)>0
            for ic=1:length(chan)
                chanHemisphere(ic)= Atlas.hemisphere(ismember(Atlas.chanName,chan(ic)))
            end
        else
            chanHemisphere=[];
        end

        WordtoAnalyze(isubtem).chan=chan;
        WordtoAnalyze(isubtem).chanHemisphere=chanHemisphere;        
        WordtoAnalyze(isubtem).ATLchan=ATLchan;
        WordtoAnalyze(isubtem).superiorATLchan=superiorATLchan;
        WordtoAnalyze(isubtem).middleATLchan=middleATLchan;
        WordtoAnalyze(isubtem).inferiorATLchan=inferiorATLchan;
        WordtoAnalyze(isubtem).MTLchan=MTLchan;
        WordtoAnalyze(isubtem).Auditorychan=Auditorychan;       
        WordtoAnalyze(isubtem).chanNum=[length(chan)];
        
    end
    WordtoAnalyze(isubtem).subID=isub;
    WordtoAnalyze(isubtem).meanACC = SubjTable(isub).MeanACC;
end

save WordtoAnalyze_ATL_2ROIs.mat WordtoAnalyze  memratio

%% select the subjects based on availbale number of trials, and available numbers of electrodes
% only analyze subject with greater 5% accuracy. 
% Memory intrusion analysis

uniquewordsID_10=PAL_Memo.uniquewordsID_10;
MissRetrievedWords=[];
medianMemor = median(PAL_Memo.Responsememorability);

MemorabilityMisReportMean(1:66)=nan;

for isub=allsub %1:length(allsub)
    
    if size(SubjTable(isub).PALTable,1)>1
        
        MissRetrievedWords=[SubjTable(isub).PALTable.response(~[SubjTable(isub).PALTable.correct{1:end}])];
        AccResponseWords= unique([SubjTable(isub).PALTable.response([SubjTable(isub).PALTable.correct{1:end}]>0)]);
        
        for i=1:length(MissRetrievedWords)
            if strcmp(num2str(MissRetrievedWords{i}),'-999')| strcmp(MissRetrievedWords{i},'<>')
                MissRetrievedWords{i} = 'NORESPONSE';
            end
        end
        
        uniquemisrepord=unique(MissRetrievedWords);
        
        for iq=1:length(uniquemisrepord)
            if  any(ismember(uniquewordsID_10,uniquemisrepord(iq)))
                MemorabilityMisReport(iq)= PAL_Memo.Responsememorability(ismember(uniquewordsID_10,uniquemisrepord(iq)));
            else
                MemorabilityMisReport(iq)=NaN;
            end 
        end

        MemorabilityMisReportMean(isub)=nanmean(MemorabilityMisReport);
        
         HigMemCount(isub)=sum(MemorabilityMisReport(~isnan(MemorabilityMisReport))>=medianMemor);
         LowMemCount(isub)=sum(MemorabilityMisReport(~isnan(MemorabilityMisReport))<medianMemor);
         RemainsMemCount(isub)=length(uniquemisrepord)-HigMemCount(isub)- LowMemCount(isub);
         allcountintro(isub)=length(uniquemisrepord);
        
         
        [h,p(isub),ci,stat]=ttest(MemorabilityMisReport(~isnan(MemorabilityMisReport)),medianMemor,'tail','right');
        r(isub)= stat.tstat/abs(stat.tstat)* sqrt(stat.tstat^2/(stat.tstat^2+stat.df));

    end
    MissRetrievedWords=[];
    MemorabilityMisReport =[];   
end


% rs=r(WordtoAnalyzeTReduced.subID);
% rs=rs(~isnan(rs));
% zrs=0.5*log((1+rs)./(1-rs));
% [h,z,ci]=signrank(zrs)

[h,z,ci]=signrank(MemorabilityMisReportMean(~isnan(MemorabilityMisReportMean)),medianMemor)
toincludeid=find(~isnan(MemorabilityMisReportMean));
nter=2000;
for i=1:nter
   boostrap_median(i)= median(MemorabilityMisReportMean(randsample(toincludeid,length(toincludeid),'true')));
end
p=1-sum((boostrap_median-medianMemor)>0)/nter
figure;hist(boostrap_median,11)
hold on;plot([ medianMemor medianMemor], [0 200],'r-','linewidth',5)


%%

% WordtoAnalyzeTReduced=WordtoAnalyzeTReduced([WordtoAnalyzeTReduced.meanACC{1:end}]>0.05,:);
WordtoAnalyzeT=struct2table(WordtoAnalyze);
WordtoAnalyzeTReduced=WordtoAnalyzeT(memratio>0,:);
WordtoAnalyzeTReduced=WordtoAnalyzeTReduced([WordtoAnalyzeTReduced.chanNum{1:end}]>=3,:);
WordtoAnalyzeTReduced=WordtoAnalyzeTReduced(WordtoAnalyzeTReduced.subID~=35,:);
% WordtoAnalyzeTReduced=WordtoAnalyzeTReduced([WordtoAnalyzeTReduced.WordCount{1:end}]>60,:);

save WordtoAnalyzeTReduced.mat WordtoAnalyzeTReduced

%%
%%- FILTERING OPTIONS
WHICH_DECOMP    = 'wavelet';  %'wavelet' 'multitaper' 'multitaper200'  % MT assumes a 500ms window, use mt200 for 200ms window; 'hilbert' not supported yet
WHICH_BANDS     = 'standard';  %WHICH_BANDS = 'standard' 'NoiseAnalysis' 'LowHigh' 'HighMidLow' 'HighMidLow2' 'JustHigh'
%- select bands for decomposition
[waveletFreqs waveletWidth waveletFreqLabels freqBandAr freqBandYticks freqBandYtickLabels hilbertFreqs multitaperObj decompStruct] = prepFreqList(WHICH_BANDS,WHICH_DECOMP);

% freqBandYticksSplit=[exp(linspace(log(2),log(16),16));exp(linspace(log(70),log(150),16))]; % only bin the data into 3 bins.

freqBandYticksSplit=[exp(linspace(log(4),log(16),10));exp(linspace(log(70),log(150),10))]; % only bin the data into 3 bins.

waveletFreqs2 = reshape(freqBandYticksSplit',1,20);
   
waveletFreqs3 = exp(linspace(log(4),log(150),30));

freqBandYticks = [2   4   8   16   32   70   150]; % overwrite the old ticks.
WordtoAnalyzeTReduced
mkdir('ATL_datafile')
curpath=pwd;
outpath=fullfile(curpath,'ATL_datafile');
% loop through subjects.
%%
MemPercent=prctile(PAL_Memo.Responsememorability,[0 33 67 100]);

for isubtem=1:size(WordtoAnalyzeTReduced,1)
    
    clear SubjResults
  
    isub =  WordtoAnalyzeTReduced.subID(isubtem)
    % duration variable across subjects
    % nanmean(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>50)) + 500 + 1000
    duration = 4000;
    offset   = -2000;
    buffer   = 500; % the butter is built-in in the duration
    resamp   = 500;
    time  =  downsample(linspace(offset,duration+offset,duration),1000/resamp);
   
    chan=WordtoAnalyzeTReduced.chan{isubtem};
%     electrodeToAna=(ismember(WordtoAnalyzeTReduced.chan{isubtem},WordtoAnalyzeTReduced.inferiorATLchan{isubtem})...
%          | ismember(WordtoAnalyzeTReduced.chan{isubtem},WordtoAnalyzeTReduced.middleATLchan{isubtem}));
     % only analyze the middle and inferior temporal gryus. The superior sites may be influenced by sounds
%          ismember(WordtoAnalyzeTReduced.chanHemisphere{isubtem},{'Right'})';
            
% only analyze the middle and inferior temporal gryus. The superior sites may be influenced by sounds
    electrodeToAna=ismember(WordtoAnalyzeTReduced.chan{isubtem},WordtoAnalyzeTReduced.ATLchan{isubtem}); 
%     electrodeToAna=ismember(WordtoAnalyzeTReduced.chan{isubtem},WordtoAnalyzeTReduced.MTLchan{isubtem}); 
    chan=chan(electrodeToAna);
%     
if length(chan)>=3 % Ony analyze the data from available channel. 
    
    temChans=[];
    allcleanTrials=[];
    
    TEST_PROBE_Event=WordtoAnalyzeTReduced.SubjTable{isubtem}.alignedEvents.TEST_PROBE(WordtoAnalyzeTReduced.InDX{isubtem});
    Recall_Event=WordtoAnalyzeTReduced.SubjTable{isubtem}.alignedEvents.RecallEvent(WordtoAnalyzeTReduced.InDX{isubtem});
 
    allssession=[TEST_PROBE_Event.session];
    uniquesessions=unique(allssession);
    
    for iss=1:length(uniquesessions)
        iss
        sessionGA_Probe=gete_ms('global_avg_good',TEST_PROBE_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[180 ],'low',2,resamp);
        sessionGA_Recall=gete_ms('global_avg_good',Recall_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[180 ],'low',2,resamp);
        
        clear ProbeRawEEG ProbePowerlog WavePow_log RecallRawEEG
        parfor ic=1:length(chan)
            % get EEG
            ic
            ProbeRawEEG(ic,:,:)=gete_ms(chan{ic},TEST_PROBE_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[180 ],'low',2,resamp);
            ProbeRawEEG(ic,:,:)=squeeze(ProbeRawEEG(ic,:,:))-sessionGA_Probe; % reference to global average of all electrode within that session.'
            
            RecallRawEEG(ic,:,:)=gete_ms(chan{ic},Recall_Event(allssession==uniquesessions(iss)),duration,offset,buffer,[180 ],'low',2,resamp);
            RecallRawEEG(ic,:,:)=squeeze(RecallRawEEG(ic,:,:))-sessionGA_Recall; % reference to global average of all electrode within that session.'
        end
        
        ProbeRawEEGAll{iss} = ProbeRawEEG;
        RecallRawEEGAll{iss} = RecallRawEEG;
        
        clfig=figure(2000);clf
        clnWeights=[2.3 2.3];
        FIG_TITLE= ['Subj' num2str(isub) '--Session' num2str(iss) 'Probe_cleaning'];
        [iChanClean1,iEvClean1,strClean] = jwCleanEEGevents_v01(ProbeRawEEG,clfig,FIG_TITLE,clnWeights);
        disp(strClean);
        print(clfig,fullfile(outpath,[FIG_TITLE '.pdf']),'-dpdf','-bestfit');
        
        clfig=figure(2001);clf
        FIG_TITLE= ['Subj' num2str(isub) '--Session' num2str(iss) 'Recall_cleaning'];
        [iChanClean2,iEvClean2,strClean] = jwCleanEEGevents_v01(RecallRawEEG,clfig,FIG_TITLE,clnWeights);
        disp(strClean);
        print(clfig,fullfile(outpath,[FIG_TITLE '.pdf']),'-dpdf','-bestfit');
        
        temChans=[temChans iChanClean1 iChanClean2];
        
        temTrials=1:length([TEST_PROBE_Event(allssession==uniquesessions(iss))]);
        allcleanTrials=[allcleanTrials ismember(temTrials,iEvClean1) & ismember(temTrials,iEvClean2)];
        SessionCleanTrials{iss}=temTrials(ismember(temTrials,iEvClean1) & ismember(temTrials,iEvClean2));
        
        clear sessionGA ProbeRawEEG
    end
    allcleanChans=unique(temChans);
    
   
    %% only did the time frequency analysis on clean trials and clean
    % electrodes
    ProbePowerlogAllSession=[]; clear ProbePowerlog
    RecallPowerlogAllSession=[]; clear RecallPowerlog
   
    % Time-frequency
    for iss=1:length(uniquesessions)
        iss
        ProbeCleanEEG{iss}=ProbeRawEEGAll{iss}(allcleanChans,SessionCleanTrials{iss},:); % only include clean chans in both sessions for analysis.
        RecallCleanEEG{iss}=RecallRawEEGAll{iss}(allcleanChans,SessionCleanTrials{iss},:); % only include clean chans in both sessions for analysis
               
        clear ProbePowerlog RecallPowerlog
        parfor ic=1:size(allcleanChans,2)
            % get
            ic
            [~,WavePow_Raw] = multiphasevec3(waveletFreqs2,squeeze(ProbeCleanEEG{iss}(ic,:,:)),resamp,waveletWidth);
            WavePow_log = 10*log10(WavePow_Raw);%log transformed raw power
            % zscoring normalize within each frequency band before averaging
            baseline = mean(WavePow_log(:,:,time>-500 & time<-300),3);
            ProbePowerlog(ic,:,:,:)  = bsxfun(@rdivide,bsxfun(@minus,WavePow_log,mean(baseline)),std(baseline)); % bsxfun(@minus,WavePow_log,mean(baseline));%
            % electrodes, trial, timepoints.
            
            [~,WavePow_Raw] = multiphasevec3(waveletFreqs2,squeeze(RecallCleanEEG{iss}(ic,:,:)),resamp,waveletWidth);
            WavePow_log = 10*log10(WavePow_Raw);%log transformed raw power
            % zscoring normalize within each frequency band before averaging
            RecallPowerlog(ic,:,:,:)  = bsxfun(@rdivide,bsxfun(@minus,WavePow_log,mean(baseline)),std(baseline)); % bsxfun(@minus,WavePow_log,mean(baseline));%
            % electrodes, trial, timepoints.
        end
        ProbePowerlogAllSession=[ProbePowerlogAllSession single(ProbePowerlog)];
        RecallPowerlogAllSession=[RecallPowerlogAllSession single(RecallPowerlog)];
        
        clear ProbePowerlog RecallPowerlog
    end
    
    
    ProbePowerlogAllSession_s = smoothdata(ProbePowerlogAllSession,4,'gaussian',250);
    RecallPowerlogAllSession_s = smoothdata(RecallPowerlogAllSession,4,'gaussian',250);
       
    SubjResults.allcleanChans=chan(allcleanChans);
    SubjResults.MeanRT=nanmean(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>50));
    SubjResults.MedianRT=nanmedian(SubjTable(isub).RTtemp(SubjTable(isub).RTtemp>0));
    SubjResults.allcleanTrials=allcleanTrials;  
   
    %%
    timewithoutbuffer = linspace(-2000,2000,2000);
    % moving window approach, no overlaps, with 20 ms gap.
    points_per_win=round(20/(1000/resamp));
    points_per_slide=round(20/(1000/resamp)); % data has been downsampled!
    BinTime=windowed_average(timewithoutbuffer,points_per_win,points_per_slide);
    
    clear ProbePowerlogAllSessionNormtime MeanATLPower_Probe RecallPowerlogAllSessionNormtime ProbePowerlogAllSessionNormtime
    parfor ic=1:size(ProbePowerlogAllSession,1)
        ic
        ProbePowerlogAllSessionNormtime(ic,:,:,:)=single(windowed_average_with_freq(squeeze(ProbePowerlogAllSession_s(ic,:,:,:)),points_per_win,points_per_slide));
        RecallPowerlogAllSessionNormtime(ic,:,:,:)=single(windowed_average_with_freq(squeeze(RecallPowerlogAllSession_s(ic,:,:,:)),points_per_win,points_per_slide))
    end
    
    tasktime = BinTime >= -1500 & BinTime <= 1500;
     
    clear Low_Probe_DataToDecode Low_Recall_DataToDecode High_Probe_DataToDecode High_Recall_DataToDecode Probe_DataToDecode Recall_DataToDecode
    % all low frequency 
    for itrl=1:size(ProbePowerlogAllSessionNormtime,2)
        for itime=1:size(ProbePowerlogAllSessionNormtime,4)
            Low_Probe_DataToDecode(itrl,:,itime)=reshape(ProbePowerlogAllSessionNormtime(:,itrl,waveletFreqs2<16.5,itime),size(allcleanChans,2)*size(ProbePowerlogAllSessionNormtime,3)/2,1);
            Low_Recall_DataToDecode(itrl,:,itime)=reshape(RecallPowerlogAllSessionNormtime(:,itrl,waveletFreqs2<16.5,itime),size(allcleanChans,2)*size(RecallPowerlogAllSessionNormtime,3)/2,1);
            
            High_Probe_DataToDecode(itrl,:,itime)=reshape(ProbePowerlogAllSessionNormtime(:,itrl,waveletFreqs2>16.5,itime),size(allcleanChans,2)*size(ProbePowerlogAllSessionNormtime,3)/2,1);
            High_Recall_DataToDecode(itrl,:,itime)=reshape(RecallPowerlogAllSessionNormtime(:,itrl,waveletFreqs2>16.5,itime),size(allcleanChans,2)*size(RecallPowerlogAllSessionNormtime,3)/2,1);
        
            Probe_DataToDecode(itrl,:,itime)=reshape(ProbePowerlogAllSessionNormtime(:,itrl,:,itime),size(allcleanChans,2)*size(ProbePowerlogAllSessionNormtime,3),1);
            Recall_DataToDecode(itrl,:,itime)=reshape(RecallPowerlogAllSessionNormtime(:,itrl,:,itime),size(allcleanChans,2)*size(RecallPowerlogAllSessionNormtime,3),1);
        end
    end
    
   SubjResults.Low_Probe_DataToDecode=Low_Probe_DataToDecode;
   SubjResults.Low_Recall_DataToDecode=Low_Recall_DataToDecode;
   SubjResults.High_Probe_DataToDecode=High_Probe_DataToDecode;
   SubjResults.High_Recall_DataToDecode=High_Recall_DataToDecode;   
   
   SubjResults.BinTime=BinTime;
   SubjResults.tasktime=tasktime;
   alltimetasktimeinx=find(tasktime);   
   SubjResults.alltimetasktimeinx=alltimetasktimeinx;   
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
   
   savefilename=fullfile(outpath,['SubjResults_ATL' num2str(isub) '.mat']);
   save(savefilename,'SubjResults','-v7.3');
        
end 

end



%% decode using SubjResults


    clear Memo Corr

    iter=1;
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
plot(BinTime(tasktime),(Memo.Low_Probe_DecodeACC));
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

for isub = 1:size(SubjResults,2)
    
    
    if ~isempty(SubjResults.allcleanChans)
        
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
                'alpha',0.25,...
                'color',[0.10 0.10 0.10],...
                'surf',s,'label', elecs{j})
            
            %   mesh_clust_d - If you are passing mesh index lists, by default a point is placed at the first mesh index in the list. However, for sulcus-
            %             spanning electrodes, we may want to plot 2+ points at each group of mesh indices. Pass a distance here if you want to clusterize
            %             such that mesh_indices >=~ mesh_clust_d apart will be plotted as separate points. Default is 3. 0 to just use the first index.
            % ,...
            
            
        end
        
        keyboard
    end
    
end



% Electrodes to exclude for misalignment
% NIH026 G27
% NIH029 G17
% NIH032 ROF3
% NIH036 OF3, OF4
% NIH062 OF4
% NIH066 TG84, TG123


