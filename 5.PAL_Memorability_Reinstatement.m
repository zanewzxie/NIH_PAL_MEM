
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

curpath=pwd;
outpath=fullfile(curpath,'data_ATL_3Event_100');
datapath=outpath;


%%
%Analysis using SubjResults

allsub=[26:64 66];
rootEEGdir = '/Volumes/Zane/NIH_FRNU_ROOT';                      %office-local

load('SubjTable_palRAMword.mat')
load('PAL_Memo_PALRAM.mat')
load('WordtoAnalyzeTReduced.mat')

ToAnalyzeSub=allsub;
medianMemor = median(PAL_Memo.Responsememorability);

fprintf('loading SubjResults_ATL...\n')

% select subject
ATLIDs=WordtoAnalyzeTReduced.subID;
accutrials=[WordtoAnalyzeTReduced.WordCount{1:end}].*[WordtoAnalyzeTReduced.meanACC{1:end}];
inacctrial=[WordtoAnalyzeTReduced.WordCount{1:end}].*(1-[WordtoAnalyzeTReduced.meanACC{1:end}]);
ATLIDs=ATLIDs(accutrials >= 10 & inacctrial>=10 & ATLIDs'~=64 ); % 37 does not have sufficent high vs. low memorable trials, %64 can only be used for beh. 

%% Calculate reinstatement group effects, to isolate the encoding time to average. 


clear Study_Probe_Sim_MemFail_ListAverage Study_Probe_Sim_MemFail_ListAverage_Sh STUDY_Probe_Similarity_List_Value_sub
for isubtem = 1:length(ATLIDs)
    
    isub=ATLIDs(isubtem);
    
    isub
    savefilename=fullfile(datapath,['SubjResults_ATL' num2str(isub) '.mat']);
   
    %%
    if exist(savefilename)>0
        
        clear SubjResults
        
        load(savefilename);        
        
        % Get the condition codes. 
        cleanmask=SubjResults.allcleanTrials;
        HighMemMask=SubjResults.HighMemMask;
        LowMemMask=SubjResults.LowMemMask;
        AccuracyMask=SubjResults.AccuracyMask;
        
        % obtain list informaiton
        Curlist=cell2mat(SubjTable(isub).PALTable.list);
        sessions=[[SubjTable(isub).alignedEvents.RecallEvent.session]];
        %         Mem = [WordtoAnalyzeT2.ResponseMemorability{isubtem}+WordtoAnalyzeT2.ProbeMemorability{isubtem}]/2;
        uniquelist = unique(Curlist);
        uniquesession=unique(sessions);
        Curlistupdate=[];
        for iss=1:length(unique(uniquesession))
            if iss==1
                Curlistupdate = [Curlistupdate;Curlist(sessions==uniquesession(iss))];
            else
                Curlistupdate = [Curlistupdate;Curlist(sessions==uniquesession(iss))*100];
            end
        end
        uniquelist=unique(Curlistupdate);
        
        
        % Memorability 
        clear ProbememPreict ExpectmemPreict Item1mem Item2mem
        for iw=1:length(SubjTable(isub).PALTable.probe)
            ProbememPreict(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.probe(iw)));
            ExpectmemPreict(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.expected(iw)));
            Item1mem(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.item1(iw)));
            Item2mem(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.item2(iw)));
        end
        AverageStudyMem = max(ProbememPreict,ExpectmemPreict);
        medianMemor=median(PAL_Memo.Responsememorability);
        
        ExpectmemPreict=ExpectmemPreict(cleanmask>0);
        AverageStudyMem=AverageStudyMem(cleanmask>0);
        
                
        % Overall group level memorability testing
        encMat = SubjResults.Study_DataToDecode;
        recMat = SubjResults.Recall_DataToDecode;
        STUDY_Probe_Similarity = ReinstatementScore(encMat,recMat);
        BinTime=SubjResults.BinTime;
        
        ReducedTimeX=BinTime(BinTime>=-4000 & BinTime<=1000);
        ReducedTimeY=BinTime(BinTime>=-1000 & BinTime<=4000);

        %Obtain across subject average
        AccurateTrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1,:,:)));
        Gl_Accurate_highmem_TrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 & ExpectmemPreict>medianMemor,:,:)));
        Gl_Accurate_lowmem_TrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 & ExpectmemPreict<medianMemor,:,:)));
        InAccurateTrialsSim(isubtem,:,:)=squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==0,:,:)));
        
        MedianRT(isubtem)=prctile(SubjResults.cleanRTs(AccuracyMask==1),50);
        MeanRT(isubtem)=mean(SubjResults.cleanRTs(AccuracyMask==1));

        Curlistupdate=Curlistupdate(cleanmask>0);
        
        t_Enctmask=BinTime>=348 & BinTime<=3248;
        t_Rectmask=BinTime>=-MedianRT(isubtem) & BinTime<=0;        
        
        % List level memorability analaysis
        clear STUDY_Probe_Similarity_List_Corr  ListMemorability STUDY_Probe_Similarity_List_Strength 

        STUDY_Probe_Similarity_List_Value=nan(length(uniquelist),5);
        for ilist=1:length(uniquelist)
            Curlist=find(Curlistupdate==uniquelist(ilist));
            FailedCurlist=find(Curlistupdate==uniquelist(ilist) & AccuracyMask'==0);
            
            TempItemRetrieval_Corr=nan;TempItemRetrieval=nan;
            for iitem=1:length(FailedCurlist)
                clear encMat recMat encMat_low encMat_high HighMem_VS_LowMem TempItemRetrieval_Corr_Sh TempItemRetrieval_Corr
                RecMat = SubjResults.Recall_DataToDecode(FailedCurlist(iitem),:,:); 

                RemainingItems=Curlist(Curlist~=FailedCurlist(iitem));
                RemainingItemMem = AverageStudyMem(RemainingItems);
                [sorted_RemainingItemMem, sortinx]=sort(RemainingItemMem);
                
                clear tempListMemSim tempListMemSimMean
                for irm=1:length(RemainingItems)
                 tempListMemSim(irm,:,:)=squeeze(ReinstatementScore(RecMat,SubjResults.Study_DataToDecode(RemainingItems(sortinx(irm)),:,:)));
                 tempListMemSimMean(irm)= mean(mean(tempListMemSim(irm,t_Enctmask,t_Rectmask),2),3);
                end
                
                [~,Rank_MemSim]= sort(tempListMemSimMean);
                
                if length(tempListMemSimMean)<5
                    TemForgItemMem_SimSORTED(iitem,:)=[tempListMemSimMean nan(1,5-length(tempListMemSimMean))];
                    TempItemRetreieval_Rank(iitem,:)=[Rank_MemSim nan(1,5-length(Rank_MemSim))];
                else
                    TemForgItemMem_SimSORTED(iitem,:)=tempListMemSimMean;
                    TempItemRetreieval_Rank(iitem,:)=Rank_MemSim;
                end
                
                TempItemRetrieval(iitem)=mean(tempListMemSimMean);
                TempItemRetrieval_Corr(iitem)=fisherz(corr(sorted_RemainingItemMem',tempListMemSimMean','type','Spearman'));
                TempItemRetrieval_Corr_Sh(iitem)=fisherz(corr(Shuffle(RemainingItemMem)',tempListMemSimMean','type','Spearman'));
                
            end

            STUDY_Probe_Similarity_List_Rank(ilist,:)=nanmean(TempItemRetreieval_Rank);
            
            STUDY_Probe_Similarity_List_Value(ilist,:)=nanmean(TemForgItemMem_SimSORTED);
            STUDY_Probe_Similarity_List_Strength(ilist)=nanmean(TempItemRetrieval);            
            STUDY_Probe_Similarity_List_Corr(ilist)=nanmean(TempItemRetrieval_Corr);     
            STUDY_Probe_Similarity_List_Corr_Sh(ilist)=nanmean(TempItemRetrieval_Corr_Sh);              
            
        end
        
        STUDY_Probe_Similarity_List_Rank_ListAverage(isubtem,:)=nanmean(STUDY_Probe_Similarity_List_Rank);
        Study_Probe_Sim_MemFail_ListAverage(isubtem)=nanmean(STUDY_Probe_Similarity_List_Corr);
        Study_Probe_Sim_MemFail_ListAverage_Sh(isubtem)=nanmean(STUDY_Probe_Similarity_List_Corr_Sh);
        STUDY_Probe_Similarity_List_Value_sub(isubtem,:)=nanmean(STUDY_Probe_Similarity_List_Value); %STUDY_Probe_Similarity_List_Rank
    end

    
end



%% calculate the cluster-based significant time window.
ReducedTime_Acc_Trial_Sim=AccurateTrialsSim(:,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000 );
ReducedTime_InAcc_Trial_Sim=InAccurateTrialsSim(:,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000 );

ReducedTime_Gl_Accurate_highmem_TrialsSim=Gl_Accurate_highmem_TrialsSim(:,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000 );
ReducedTime_Gl_Accurate_lowmem_TrialsSim=Gl_Accurate_lowmem_TrialsSim(:,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000 );

ReducedTimeX=BinTime(BinTime>=-4000 & BinTime<=1000);
ReducedTimeY=BinTime(BinTime>=-1000 & BinTime<=4000);
% 
% ReducedTime_Xneural=Target_NonT_Highmem_sim(:,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000 );
% ReducedTime_Yneural=Target_NonT_Lowmem_sim(:,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000 );


% for overall accurate vs. inaccurate trials
X=ReducedTime_Acc_Trial_Sim;
Y=ReducedTime_InAcc_Trial_Sim;
xaxislabel=ReducedTimeX;
yaxislabel=ReducedTimeY;
medianRTs=mean(MedianRT);
cluster_based_correction

h=gcf;
set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
printfilename=fullfile(outpath,['Reinstate_Allsub.pdf']);
print(h,printfilename,'-dpdf','-bestfit');

% Two critical variable come out from this analysis
%[min(signYtime) max(signYtime)]

% signYtime=[350 2950];
% signYtime=[1000 3000];
%%
clear EearlyNeuMem LateNeuMem  NeurMemTrace NeurMemTrace NeuMemCorr


for isubtem = 1:length(ATLIDs)
    
    isub=ATLIDs(isubtem);
    
     isub
     savefilename=fullfile(datapath,['SubjResults_ATL' num2str(isub) '.mat']);

    if exist(savefilename)>0
        
        clear SubjResults
        load(savefilename);
        
        % extract all condition code
        Curlist=[SubjTable(isub).PALTable.list{1:end}];
        sessions=[[SubjTable(isub).alignedEvents.RecallEvent.session]];

        % Memorability 
        clear ProbememPreict ExpectmemPreict Item1mem Item2mem
        for iw=1:length(SubjTable(isub).PALTable.probe)
            ProbememPreict(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.probe(iw)));
            ExpectmemPreict(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.expected(iw)));
            Item1mem(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.item1(iw)));
            Item2mem(iw) =  PAL_Memo.Responsememorability(ismember(PAL_Memo.uniquewordsID_10,SubjTable(isub).PALTable.item2(iw)));
        end
        AverageStudyMem = (ProbememPreict+ExpectmemPreict)/2;
        medianMemor=median(PAL_Memo.Responsememorability);
        
        %list based analysis
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

        % 
        cleanmask=SubjResults.allcleanTrials;
        AverageStudyMem = AverageStudyMem(cleanmask>0);
        ExpectmemPreict = ExpectmemPreict(cleanmask>0);
        ProbememPreict=ProbememPreict(cleanmask>0);
        ListMem=ListMem(cleanmask>0);

        AccuracyMask=AccuracyMask(cleanmask>0);
        TrialsPloted=[sum(AccuracyMask==1 & ExpectmemPreict<=medianMemor) sum(AccuracyMask==1 & ExpectmemPreict>medianMemor)];
        min([sum(AccuracyMask==1 & ExpectmemPreict<=medianMemor) sum(AccuracyMask==1 & ExpectmemPreict>medianMemor)])

%%

% subject level 
%Change the dimensions to be feature X time X event
encMat = SubjResults.Study_DataToDecode;
recMat = SubjResults.Recall_DataToDecode;
numEvents=size(encMat,1)
clear output STUDY_Probe_Similarity
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


% X=STUDY_Probe_Similarity(AccuracyMask==1,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000);
% Y=STUDY_Probe_Similarity(AccuracyMask==0,BinTime>=-1000 & BinTime<=4000,BinTime>=-4000 & BinTime<=1000);
% xaxislabel=ReducedTimeX;
% yaxislabel=ReducedTimeY;
% zmapthresh=ClusterCorrect_CompBL_WithinSub(X,Y,xaxislabel,yaxislabel,MedianRT(isubtem));
% signYtime=yaxislabel(mean(zmapthresh,2)>0);
% signXtime=xaxislabel(mean(zmapthresh,1)>0);
% subjectTOI(isubtem).signYtime=signYtime;
% subjectTOI(isubtem).signXtime=signXtime;
% hold on;
% plot([min(xaxislabel) max(xaxislabel)],[min(signYtime) min(signYtime)],'r-')
% plot([min(xaxislabel) max(xaxislabel)],[max(signYtime) max(signYtime)],'r-')
% plot([min(signXtime) min(signXtime)],[min(yaxislabel) max(yaxislabel)],'r-')
% plot([min(signXtime) min(signXtime)],[min(yaxislabel) max(yaxislabel)],'r-')
% keyboard


% sanity check!
h=figure;
subplot(2,2,1);
xyimagesc(BinTime,BinTime,squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 ,:,:)))); %
axis([-4000 1000 -1000 4000]);
axis square

subplot(2,2,2);
tmask=BinTime>=min(signYtime) & BinTime<=max(signYtime);
xyimagesc(BinTime,BinTime(tmask),squeeze(mean(STUDY_Probe_Similarity(AccuracyMask==1 ,tmask,:)))); %
axis([-4000 1000 -1000 4000]);
axis square

medRT=prctile(SubjResults.cleanRTs(AccuracyMask==1),50);
hold on;
plot([-medRT -medRT],[min(BinTime) max(BinTime)],'r-')
plot([-0 -0],[min(BinTime) max(BinTime)],'k-')


A = squeeze(mean(STUDY_Probe_Similarity(:,tmask,:),2))

% Normalized in time 
RoundBinTime=round(BinTime/25)*25;
% normalized by trial-by-trial reaction time. 
AllAccRT=SubjResults.cleanRTs;
clear AZpercent AZ
for itr=1:size(AllAccRT,2)
    if AllAccRT(itr)==-999 | AllAccRT(itr)==0
        currentRT=-4000;
    else
        currentRT=-AllAccRT(itr);
        currentRT=round(currentRT/25)*25;
    end
    y = A(itr,RoundBinTime>=currentRT & RoundBinTime<=-currentRT);
    x = linspace(-100, 100, length(y));    
    xi = [-100:1:100]; %linspace(-100, 100, 20);    
    AZpercent(itr,:)= interp1(x,y,xi);
end


% pertime=linspace(-100, 100, 20);
pertime=[-100:1:100]
subplot(2,2,3);
plot(pertime,mean(AZpercent(AccuracyMask==1 & ExpectmemPreict>medianMemor,:)));
hold on;
plot(pertime,mean(AZpercent(AccuracyMask==1 & ExpectmemPreict<medianMemor,:)))
plot([0 0],[-0.2 0.2],'k-');axis([min(pertime) 100 -0.2 0.2])
plot([-medRT -medRT],[-0.2 0.2],'k.-');axis([min(pertime) 100 -0.2 0.2])
legend('HighMem','LowMem','recall onset','location','southeast')
title('Averaged Similarity')
xlabel('Percentage of Recall Time')
axis square

NeurMemTrace.Normalized(isubtem,1,:)=mean(AZpercent(AccuracyMask==1 & ExpectmemPreict>medianMemor,:));
NeurMemTrace.Normalized(isubtem,2,:)=mean(AZpercent(AccuracyMask==1 & ExpectmemPreict<medianMemor,:));
NeurMemTrace.Raw(isubtem,1,:)=mean(A(AccuracyMask==1 & ExpectmemPreict>medianMemor,:));
NeurMemTrace.Raw(isubtem,2,:)=mean(A(AccuracyMask==1 & ExpectmemPreict<medianMemor,:));

ax=subplot(2,2,4)
x=mean(AZpercent(AccuracyMask==1,pertime>-75 & pertime<-50),2);
y=ExpectmemPreict(AccuracyMask==1)';
p = polyfit(x,y,1);
f = polyval(p,x);
plot(x,y,'bo',x,f,'b-')
NeuMemCorr.Early_r(isubtem) = corr(x,y,'type','spearman');
NeuMemCorr.Early_N(isubtem) = length(x);
axis square

hold on
x=mean(AZpercent(AccuracyMask==1,pertime>-25 & pertime<0),2);
y=ExpectmemPreict(AccuracyMask==1)';
p = polyfit(x,y,1);
f = polyval(p,x);
plot(x,y,'ro',x,f,'r-')
legend('early 50%','','late 50%','')
xlabel('Neural Reinstatement - AccTrial')
ylabel('Observed Memorability - AccTrail')
NeuMemCorr.Late_r(isubtem) = corr(x,y,'type','spearman');
NeuMemCorr.Late_N(isubtem) = length(x);


x=mean(AZpercent(AccuracyMask==1,:),2);
y=ExpectmemPreict(AccuracyMask==1)';
NeuMemCorr.Overall_r(isubtem) = corr(x,y,'type','spearman');
NeuMemCorr.Overall_N(isubtem) = length(x);


    h=gcf;
    set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
    printfilename=fullfile(outpath,['Reinstate_Mem_Sub' num2str(isub) '.pdf']);
    print(h,printfilename,'-dpdf','-bestfit');
    
    end
    
end


%% stat

NSUB=size(NeurMemTrace.Raw,1);
for iter=1:1000
    randsub=randsample(1:NSUB,NSUB,'true');
    LateR(iter)=mean(fisherz(NeuMemCorr.Late_r(randsub)));
    EearlR(iter)=mean(fisherz(NeuMemCorr.Early_r(randsub)));
    DiffR(iter)=EearlR(iter)-LateR(iter);
    
    FailMem_rankcorr(iter)=mean(Study_Probe_Sim_MemFail_ListAverage(randsub));
    FailMem_rankcorr_sh(iter)=mean(Study_Probe_Sim_MemFail_ListAverage_Sh(randsub));

end


%%
% Use a resampling approach
MedianRTs=medianRTs;

HighMemR=squeeze((NeurMemTrace.Raw(:,1,:)));
LowMemR=squeeze((NeurMemTrace.Raw(:,2,:)));Study_Probe_Sim_MemFail_ListAverage
figure;
subplot(1,2,1);
prettyplot(BinTime,mean(HighMemR),std(HighMemR)/sqrt(size(HighMemR,1)),[0 0 1],1)
axis([-medianRTs-500 1000 0 0.1])
hold on;
subplot(1,2,1);prettyplot(BinTime,mean(LowMemR),std(LowMemR)/sqrt(size(LowMemR,1)),[1 0 0],1)
title('Raw Reinstatement Score')

HighMem=squeeze((NeurMemTrace.Normalized(:,1,:)));
LowMem=squeeze((NeurMemTrace.Normalized(:,2,:)));


subplot(1,2,2);
prettyplot(pertime,mean(HighMem),std(HighMem)/sqrt(size(HighMem,1)),[0 0 1],1)
axis([-100 50 0 0.1])
hold on;
subplot(1,2,2);prettyplot(pertime,mean(LowMem),std(LowMem)/sqrt(size(HighMem,1)),[1 0 0],1)
title('Normalized Reinstatement Score')


% find peak, jackknife
ALLSUB=1:NSUB;
clear PeakTime1 PeakTime2

for iter=1:NSUB
Y1=mean(HighMemR(ALLSUB~=iter,:));
Y2=mean(LowMemR(ALLSUB~=iter,:));

K1 = dsearchn(Y1',max(Y1));
PeakTime1(iter)=BinTime(K1);

K2 = dsearchn(Y2',max(Y2));
PeakTime2(iter)=BinTime(K2);
end

[h,p,ci,stat]=ttest(PeakTime2,PeakTime1)

t=stat.tstat/stat.df
p=tpdf(t,stat.df)

%%

% raw
figure;
for iter=1:1000
Y1=mean(HighMemR(randsample(1:NSUB,NSUB,'true'),:));
Y2=mean(LowMemR(randsample(1:NSUB,NSUB,'true'),:));

subplot(1,2,1)
hold on;plot(BinTime,Y1);
subplot(1,2,2)
hold on;plot(BinTime,Y2);

tmask=BinTime <1000 & BinTime>-MedianRTs;
X=BinTime(tmask);
% 
K1 = dsearchn(Y1',max(Y1(tmask)));
PeakTime1(iter)=BinTime(K1);

K2 = dsearchn(Y2',max(Y2(tmask)));
PeakTime2(iter)=BinTime(K2);

CumQ=cumtrapz(X,Y1(tmask));
K1 = dsearchn(CumQ',CumQ(end)/2);
FractTime1(iter)=X(K1);

CumQ=cumtrapz(X,Y2(tmask));
K2 = dsearchn(CumQ',CumQ(end)/2);
FractTime2(iter)=X(K2);
end

p=1-sum((PeakTime1-PeakTime2)<0)/1000
% p=1-sum((FractTime1-FractTime2)<0)/1000


%%


% normalized
figure;
for iter=1:1000
Y1=mean(HighMem(randsample(1:NSUB,NSUB,'true'),:));
Y2=mean(LowMem(randsample(1:NSUB,NSUB,'true'),:));

subplot(1,2,1)
hold on;plot(pertime,Y1);
subplot(1,2,2)
hold on;plot(pertime,Y2);

tmask=pertime <50;
X=pertime(tmask);
% 
K1 = dsearchn(Y1',max(Y1(tmask)));
PeakTime1(iter)=pertime(K1);

K2 = dsearchn(Y2',max(Y2(tmask)));
PeakTime2(iter)=pertime(K2);

CumQ=cumtrapz(X,Y1(tmask));
K1 = dsearchn(CumQ',CumQ(end)/2);
FractTime1(iter)=X(K1);

CumQ=cumtrapz(X,Y2(tmask));
K2 = dsearchn(CumQ',CumQ(end)/2);
FractTime2(iter)=X(K2);
end

p=1-sum((PeakTime1-PeakTime2)<0)/1000
% p=1-sum((FractTime1-FractTime2)<0)/1000


