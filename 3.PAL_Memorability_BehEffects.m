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

ToAnalyzeSub=allsub;

medianMemor = median(PAL_Memo.Responsememorability);

for isubtem=1:length(ToAnalyzeSub)
    isubtem
    isub=ToAnalyzeSub(isubtem);
     
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
        
        superiorATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Superior temporal gyrus')  & strcmp(Atlas.label_destrieux,'Superior temporal gyrus')); %
        middleATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Middle temporal gyrus')  & strcmp(Atlas.label_destrieux,'Middle temporal gyrus') ); %
        inferiorATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe') & strcmp(Atlas.label_desikan,'Inferior temporal gyrus') & strcmp(Atlas.label_destrieux,'Inferior temporal gyrus') ); %
        
        MTLchan     = Atlas.chanName(strcmp(Atlas.label_desikan,'Entorhinal cortex') | strcmp(Atlas.label_desikan,'Parahippocampal gyrus')); %        
        ATLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Anterior Temporal Lobe')); %       
        PTLchan     = Atlas.chanName(strcmp(Atlas.label_wittig,'Posterior Temporal Lobe'));
        
        chan=[ATLchan;MTLchan;PTLchan];
 
        %% Assign ECoG file location.       
        allwords=WordtoAnalyze(isubtem).Word;
        
        for iev=1:length(SubjTable(isub).alignedEvents.TEST_PROBE)

            SubjTable(isub).alignedEvents.STUDY_PAIR(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.STUDY_PAIR(iev).eegfile,'noreref','processed');
            SubjTable(isub).alignedEvents.STUDY_PAIR(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.STUDY_PAIR(iev).eegfile,'/Volumes/Shares/FRNU/data/eeg',rootEEGdir);
            
            SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile,'noreref','processed');
            SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.TEST_PROBE(iev).eegfile,'/Volumes/Shares/FRNU/data/eeg',rootEEGdir);
            
            SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile,'noreref','processed');
            SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile =  strrep(SubjTable(isub).alignedEvents.RecallEvent(iev).eegfile,'/Volumes/Shares/FRNU/data/eeg',rootEEGdir);
        end
        
        WordtoAnalyze(isubtem).SubjTable=SubjTable(isub);
        
        %% here exclude bad channels based on Julio's script.
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
        % to normalized
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
        PTLchan=PTLchan(ismember(PTLchan,chan));
        MTLchan=MTLchan(ismember(MTLchan,chan));
        
        superiorATLchan=superiorATLchan(ismember(superiorATLchan,chan));
        middleATLchan=middleATLchan(ismember(middleATLchan,chan));
        inferiorATLchan=inferiorATLchan(ismember(inferiorATLchan,chan));
%         Auditorychan=Auditorychan(ismember(Auditorychan,chan));
        
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
%         WordtoAnalyze(isubtem).Auditorychan=Auditorychan;   
        WordtoAnalyze(isubtem).PTLchan=PTLchan;               
        WordtoAnalyze(isubtem).chanNum=[length(chan)];
        
    end
    WordtoAnalyze(isubtem).subID=isub;
    WordtoAnalyze(isubtem).meanACC = SubjTable(isub).MeanACC;
end

save WordtoAnalyze_ATL_2ROIs.mat WordtoAnalyze  memratio

%% Memory intrusion analysis and RT analysis

uniquewordsID_10=PAL_Memo.uniquewordsID_10;
MissRetrievedWords=[];
medianMemor = median(PAL_Memo.Responsememorability);
MemorabilityMisReportMean(1:66)=nan;

for isubtem=1:length(ToAnalyzeSub) %1:length(allsub)
    
   isub=ToAnalyzeSub(isubtem);

    if size(SubjTable(isub).PALTable,1)>1
        
        %Intrusion. 
        MissRetrievedWords=[SubjTable(isub).PALTable.response(~[SubjTable(isub).PALTable.correct{1:end}])];
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
        
        HigMemCount(isub)=sum(MemorabilityMisReport(~isnan(MemorabilityMisReport))>medianMemor);
        LowMemCount(isub)=sum(MemorabilityMisReport(~isnan(MemorabilityMisReport))<medianMemor);
        RemainsMemCount(isub)=length(uniquemisrepord)-HigMemCount(isub)- LowMemCount(isub);
        allcountintro(isub)=length(uniquemisrepord);
         
        [h,p,ci,stat]=ttest(MemorabilityMisReport(~isnan(MemorabilityMisReport)),medianMemor,'tail','right');
        MemIntrusionR(isub)= stat.tstat/abs(stat.tstat)* sqrt(stat.tstat^2/(stat.tstat^2+stat.df));
        MemInstrusionN(isub)=stat.df+1;

        
% for plotting subject level content
%         x =MemorabilityMisReport(~isnan(MemorabilityMisReport));
%         figure;
%         b=bar([x mean(x)]);
%         b.FaceColor = 'flat';
%         b.CData(40,:)=[0 0 1];
%         b.CData(1:39,:)=repmat([0.5 0.5 0.5],39,1);
%         hold on;errorbar(40,mean(x),std(x)/sqrt(length(x)-1),'linewidth',2);
%         hold on;plot([ 0 40],[medianMemor medianMemor],'r-','linewidth',2);
%         axis([0 41 0.2 0.7]);
%         xlabel('Intruded words on different trials');
%         ylabel('Memorability value of intruded words');
%         title('NIH40 with mean (SE) plotted in the last bin');
%         set(gca,'fontsize',18);
        

        % Retrieval RTs. 
        % extrace the memorability score for each retrieval word. 
        clear Acc RT highmask lowmask memorability SDRT MemorabilityAccResponseWords
        Acc=[SubjTable(isub).PALTable.correct{1:end}]>0;
        RT=[SubjTable(isub).PALTable.RT{1:end}];
        AccResponseWords= ([SubjTable(isub).PALTable.response]); 
        AccResponseWords=AccResponseWords(RT>0);
        
        RemainedRT = RT(RT>0);
        
        RTSD = std(RT(RT >0));
        for i=1:length(AccResponseWords)
            if strcmp(num2str(AccResponseWords{i}),'-999')| strcmp(AccResponseWords{i},'<>')
                AccResponseWords{i} = 'NORESPONSE';
            end
        end
        
        clear MemorabilityAccResponseWords
        
        for iq=1:length(AccResponseWords)
            if  any(ismember(uniquewordsID_10,AccResponseWords(iq)))
                MemorabilityAccResponseWords(iq)= PAL_Memo.Responsememorability(ismember(uniquewordsID_10,AccResponseWords(iq)));
            else
                MemorabilityAccResponseWords(iq)=NaN;
            end
        end

        cormask = ~isnan(MemorabilityAccResponseWords);
        Correlation_MRT_s(isub)=corr(MemorabilityAccResponseWords(cormask)',log(RemainedRT(cormask))','type','pearson');
        countN(isub)=sum(cormask);
        Acc_ofsub(isub)=mean(Acc);
        
        MemRTs(isub) = median(RT(cormask));
        x=MemorabilityAccResponseWords(cormask);
        y=log(RemainedRT(cormask));
%         
%         figure(101);
%         Fit = polyfit(x,y,1)
%         hold on;scatter(x,y,'fill');
%         plot([min(x) max(x)],polyval(Fit,[min(x) max(x)]),'linewidth',10*countN(isub)/265)
%          xlabel('Memorability Values of Responded Words');
%          ylabel('Log Response Time (originally in ms)');
%          title('NIH40');
%          set(gca,'fontsize',18);
         
    end
    
    MissRetrievedWords=[];
    MemorabilityMisReport =[];
    
end

zrs = .5.*[log(1+Correlation_MRT_s)-log(1-Correlation_MRT_s)]
zMemIntrusionR=.5.*[log(1+MemIntrusionR)-log(1-MemIntrusionR)]
%% Group level analysis. 

% RT and memorability correlation
IncludeIDRT=find(countN>0);
IncludedMemRTs = MemRTs(IncludeIDRT);
[Correlation_MRT_s(IncludeIDRT)' countN(IncludeIDRT)']

[h,z,ci]=signrank(MemorabilityMisReportMean(~isnan(MemorabilityMisReportMean)),medianMemor)
toincludeid=find(~isnan(MemorabilityMisReportMean));


nter=2000;
for i=1:nter
   boostrap_median(i)= median(MemorabilityMisReportMean(randsample(toincludeid,length(toincludeid),'true')));
   boostrap_medianRT(i,:)= mean(MemRTs(randsample(IncludeIDRT,length(IncludeIDRT),'true')));
end
p=1-sum((boostrap_median-medianMemor)>0)/nter

figure;hist(boostrap_median,11)
hold on;plot([ medianMemor medianMemor], [0 200],'r-','linewidth',5)


%%

% WordtoAnalyzeTReduced=WordtoAnalyzeTReduced([WordtoAnalyzeTReduced.meanACC{1:end}]>0.05,:);
WordtoAnalyzeT=struct2table(WordtoAnalyze);
WordtoAnalyzeTReduced=WordtoAnalyzeT(memratio>0,:);
WordtoAnalyzeTReduced=WordtoAnalyzeTReduced([WordtoAnalyzeTReduced.chanNum{1:end}]>=2,:);
WordtoAnalyzeTReduced=WordtoAnalyzeTReduced(WordtoAnalyzeTReduced.subID~=35,:);

delete WordtoAnalyzeTReduced
save WordtoAnalyzeTReduced.mat WordtoAnalyzeTReduced
% WordtoAnalyzeTReduced=WordtoAnalyzeTReduced([WordtoAnalyzeTReduced.WordCount{1:end}]>60,:);

