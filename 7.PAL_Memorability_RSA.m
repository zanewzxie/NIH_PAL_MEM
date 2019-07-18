
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
outpath=fullfile('/Volumes/Zane/NIH_HPC/NIH_PAL_Mem','data_ATL_3Event_100');
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


clear Study_Probe_Sim_MemFail_ListAverage Study_Probe_Sim_MemFail_ListAverage_Sh STUDY_Probe_Similarity_List_Value_sub Probe_MemRSA
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
        cleanRTs=SubjResults.cleanRTs;
        
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
        
                
        
        clear probeMat
        % Overall group level memorability testing
        encMat = SubjResults.Low_Study_DataToDecode; % #trial by #feature by #time matrix. 
        probeMat = SubjResults.Low_Probe_DataToDecode; % #trial by #feature by #time matrix. 
        recMat = SubjResults.Low_Recall_DataToDecode; % #trial by #feature by #time matrix. 
        BinTime=SubjResults.BinTime;
        
%         
        clear ACCmodel a1
        reducedACCMASK=AccuracyMask(cleanRTs>-1000);
        for i=1:length(reducedACCMASK)
            for j=1:length(reducedACCMASK)
                 ACCmodel(i,j)=abs(reducedACCMASK(i)-reducedACCMASK(j));
            end
        end
        a1 = ACCmodel(find(tril(ones(size(ACCmodel)),-1)));
        
        clear ACCmodel a5
        reducedACCMASK=AccuracyMask(cleanRTs>-1000);
        for i=1:length(reducedACCMASK)
            for j=1:length(reducedACCMASK)
                 ACCmodel(i,j)=abs(reducedACCMASK(i)+reducedACCMASK(j))/2;
            end
        end
        a5 = ACCmodel(find(tril(ones(size(ACCmodel)),-1)));
        
        clear MemoModel a2
        MemorabilityScoreReduced=ExpectmemPreict(cleanRTs>-1000);
        for i=1:length(MemorabilityScoreReduced)
            for j=1:length(MemorabilityScoreReduced)
                    MemoModel(i,j)= abs(MemorabilityScoreReduced(i)-MemorabilityScoreReduced(j));
            end
        end
        a2 = MemoModel(find(tril(ones(size(MemoModel)),-1)));
        
        
        clear MemoModel a3
        MemorabilityScoreReduced=ExpectmemPreict(cleanRTs>-1000);
        for i=1:length(MemorabilityScoreReduced)
            for j=1:length(MemorabilityScoreReduced)
                    MemoModel(i,j)= abs(MemorabilityScoreReduced(i)+MemorabilityScoreReduced(j))/2;
            end
        end
        a3 = MemoModel(find(tril(ones(size(MemoModel)),-1)));
        

        clear MemoModel a4
        MemorabilityScoreReduced=ExpectmemPreict(AccuracyMask==0);
        for i=1:length(MemorabilityScoreReduced)
            for j=1:length(MemorabilityScoreReduced)
                    MemoModel(i,j)= abs(MemorabilityScoreReduced(i)+MemorabilityScoreReduced(j))/2;
            end
        end
        a4 = MemoModel(find(tril(ones(size(MemoModel)),-1)));
        
        for itime=1:length(BinTime)
            Ptempdata=squeeze(probeMat(cleanRTs>-1000,:,itime));
            RSA=(squareform(pdist(Ptempdata,'cosine')));
            NeuraA12 = RSA(find(tril(ones(size(RSA)),-1)));
            Probe_ACCRSA_distance(isubtem,itime)=corr(NeuraA12,a1,'type','Spearman');
            
            Ptempdata=squeeze(probeMat(cleanRTs>-1000,:,itime));
            RSA=(squareform(pdist(Ptempdata,'cosine')));
            NeuraA12 = RSA(find(tril(ones(size(RSA)),-1)));
            Probe_MemRSA_distance(isubtem,itime)=corr(NeuraA12,a2,'type','Spearman');   
            
            Ptempdata=squeeze(probeMat(cleanRTs>-1000,:,itime));
            RSA=(squareform(1-pdist(Ptempdata,'cosine')));
            NeuraA12 = RSA(find(tril(ones(size(RSA)),-1)));
            Probe_ACCRSA_norm(isubtem,itime)=corr(NeuraA12,a5,'type','Spearman');                     

            Ptempdata=squeeze(probeMat(AccuracyMask>-1000 ,:,itime));
            RSA=(squareform(1-pdist(Ptempdata,'cosine')));
            NeuraA3 = RSA(find(tril(ones(size(RSA)),-1)));
            Probe_MemRSA_norm(isubtem,itime)=corr(NeuraA3,a3,'type','Spearman');        
            
            Ptempdata=squeeze(probeMat(AccuracyMask==0 ,:,itime));
            RSA=(squareform(1-pdist(Ptempdata,'cosine')));
            NeuraA4 = RSA(find(tril(ones(size(RSA)),-1)));
            Probe_MemRSA_Incc_norm(isubtem,itime)=corr(NeuraA4,a4,'type','Spearman');  
            
        end
    
    end
    
    figure;plot(BinTime,Probe_ACCRSA_distance(isubtem,:))
    hold on;plot(BinTime,Probe_MemRSA_distance(isubtem,:))
    hold on;plot(BinTime,Probe_ACCRSA_norm(isubtem,:))
    hold on;plot(BinTime,Probe_MemRSA_norm(isubtem,:))
    hold on;plot(BinTime,Probe_MemRSA_Incc_norm(isubtem,:))    
    legend('Distance Model ACC','Distance Model Memo','Norm Model ACC','Norm Model Memo','Norm Model Memo Inacc')
    pause(1)
end


%% accuracy models
NSUB=18;
figure;
subplot(121)
prettyplot(BinTime, mean(Probe_ACCRSA_norm),std(Probe_ACCRSA_norm)/sqrt(18),[1 0 0],1)
% hold on
% prettyplot(BinTime, mean(Probe_ACCRSA_norm),std(Probe_ACCRSA_norm)/sqrt(18),[0 1 0],1)
hold on
prettyplot(BinTime, mean(Probe_ACCRSA_distance),std(Probe_ACCRSA_distance)/sqrt(18),[0 0 1],1)
axis([-500 2000 -.05 0.1])
legend('norm-based model','','Two-cluster model','')
title('Accuracy model')

subplot(122)
diffx=Probe_ACCRSA_distance-Probe_ACCRSA_norm;
for iter=1:1000
    diffiter(iter,:)=mean(diffx(randsample((1:NSUB),NSUB,'true'),:));
end
p=mean(1-(diffiter>0));
h=p<0.05;
prettyplot(BinTime, mean(Probe_ACCRSA_distance-Probe_ACCRSA_norm),std(Probe_ACCRSA_distance-Probe_ACCRSA_norm)/sqrt(18),[0 1 0],1)
axis([-500 2000 -.05 0.1])
hold on;plot(BinTime,h*0.02);
legend('Diff (Distance-Norm)')
title('Accuracy model')

        
        h=gcf;
        set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
        printfilename=fullfile(outpath,['RAS_FOR_ACC.pdf']);
        print(h,printfilename,'-dpdf','-bestfit');
        
%%
figure;
subplot(121)
prettyplot(BinTime, mean(Probe_MemRSA_norm),std(Probe_MemRSA_norm)/sqrt(18),[1 0 0],1)
% hold on
% prettyplot(BinTime, mean(Probe_ACCRSA_norm),std(Probe_ACCRSA_norm)/sqrt(18),[0 1 0],1)
hold on
prettyplot(BinTime, mean(Probe_MemRSA_distance),std(Probe_MemRSA_distance)/sqrt(18),[0 0 1],1)
axis([-500 2000 -.05 0.05])
legend('norm-based model','','Two-cluster model','')
title('Memorablity model')

subplot(122)
diffx=Probe_MemRSA_norm-Probe_MemRSA_distance;
for iter=1:1000
    diffiter(iter,:)=mean(diffx(randsample((1:NSUB),NSUB,'true'),:));
end
p=mean(1-(diffiter>0));
h=p<0.05;
prettyplot(BinTime, mean(Probe_MemRSA_norm-Probe_MemRSA_distance),std(Probe_MemRSA_norm-Probe_MemRSA_distance)/sqrt(18),[0 1 0],1)
axis([-500 2000 -.05 0.05])
hold on;plot(BinTime,h*0.02);
legend('Diff (Norm-Distance)')
title('Memorablity model')

        h=gcf;
        set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
        printfilename=fullfile(outpath,['RAS_FOR_MEMO.pdf']);
        print(h,printfilename,'-dpdf','-bestfit');
        

%%
figure;
subplot(121)
prettyplot(BinTime, mean(Probe_MemRSA_norm),std(Probe_MemRSA_norm)/sqrt(18),[1 0 0],1)
diffx=Probe_MemRSA_norm;
axis([-500 2000 -.05 0.05])
for iter=1:1000
    diffiter(iter,:)=mean(diffx(randsample((1:NSUB),NSUB,'true'),:));
end
p=mean(1-(diffiter>0));
h1=p<0.05;
hold on;plot(BinTime,h1*0.02);
title('Memorablity model')
legend('norm-based all memo','')


subplot(122)
prettyplot(BinTime, mean(Probe_MemRSA_Incc_norm),std(Probe_MemRSA_Incc_norm)/sqrt(18),[0 0 1],1)
axis([-500 2000 -.05 0.05])
legend('norm-based inaccurate trials','')
title('Memorablity model')
diffx=Probe_MemRSA_Incc_norm;
for iter=1:1000
    diffiter(iter,:)=mean(diffx(randsample((1:NSUB),NSUB,'true'),:));
end
p=mean(1-(diffiter>0));
h2=p<0.05;
hold on;plot(BinTime,h2*0.03);

        h=gcf;
        set(h,'PaperOrientation','landscape','Position',[50 50 1200 800]);
        printfilename=fullfile(outpath,['RAS_FOR_MEMO_INACC.pdf']);
        print(h,printfilename,'-dpdf','-bestfit');
        