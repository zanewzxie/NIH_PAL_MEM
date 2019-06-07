clear all
%add toolboxes
addpath(genpath('/Volumes/Zane/Matlab/eeg_toolbox/trunk'));
addpath(genpath('/Volumes/Zane/Matlab/dungeon_toolbox_17a'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/EEG_Preprocessing'));
addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/IEM_ester'));
addpath(genpath('/Volumes/Zane/NIH_NINDS/Data_InProgress/SpecFuN'));

%% extract PAL behaviroal data
allsub=26:66;
% exclude=[33 49];
rootEEGdir = '/Volumes/Zane/NIH_FRNU_ROOT';                      %office-local

delete SubjTable_palRAMword.mat

for isub= allsub %1:length(allsub)
    
SubjTable(isub).PALTable=[];
SubjTable(isub).subID=[];

if isub<10
CURSUB=['NIH00' num2str(isub)];
else
CURSUB=['NIH0' num2str(isub)];
end

isub

loadfile{1}=fullfile(rootEEGdir,[CURSUB '/behavioral/pa3']);
loadfile{2}=fullfile(rootEEGdir,[CURSUB '/behavioral/palRam']);
loadfile{3}=fullfile(rootEEGdir,[CURSUB '/behavioral/PAL']);

LoadFILE=[];
for i=1:3
    if exist(loadfile{i})>0
        LoadFILE=loadfile{i};
    end
end

clear events

if exist(fullfile(LoadFILE,'events.mat'))>0 & isub~=45 & isub~=30 & isub~=49 
    % & isub~=35 has lesion, include only for behavior but not for ECoG. 
    % & isub~= 45 and 30 speak spanish
    % & isub 33 only has 24 trials. 
    % & isub 29, 37, 38, 49, 53, 57, has fewer than 10 trial accuracy. 
    % & isub~=29 & isub~=37  & isub~=38  & isub~=49  & isub~=53 & isub~=57 
    
load(fullfile(LoadFILE,'events.mat'));

clear sortedEvents alignedTestProbe alignedEvents

if isub>=49 & isub<62
    types = {'STUDY_PAIR','PROBE_START'};
    sortedEvents = struct();
    for t = types
        sortedEvents.(t{1}) = filterStruct(events,'list~=0 & strncmp(type,varargin{1},8)',t{1});
    end
    sortedEvents.TEST_PROBE = sortedEvents.PROBE_START;
else
    types = {'STUDY_PAIR','TEST_PROBE'};
    sortedEvents = struct();
    for t = types
        sortedEvents.(t{1}) = filterStruct(events,'list~=0 & strncmp(type,varargin{1},8)',t{1});
    end
end
    
%function alignedEvents = alignEvents(sortedEvents)
%This function aligns the events so that each encoding event is aligned
%with the proper recall event
clear alignedEvents alignedStudyPair alignedTestProbe

recEvents = sortedEvents.TEST_PROBE;
sess1=[recEvents.session]
uniquesession=unique(sess1);
sess=uniquesession(1);

% here, identify unique sessions.
excludesessions=zeros(length(uniquesession),1);
clear temev list word_list12
for i=1:length(uniquesession)
    temev{i}=sortedEvents.TEST_PROBE(sess1==(uniquesession(i)));
    tempword{i}={temev{i}.probe_word};
    list{i}=[temev{i}.list];
    word_list12{i}=unique(tempword{i}(list{i}==1 | list{i}==2));
    if length( word_list12{i})<11
        excludesessions(i)=1;
    end
end

if length(uniquesession)>1
for i=2:length(uniquesession)
    if length(word_list12{i})>11
    excludesessions(i)=excludesessions(i)+all(strcmp(word_list12{i},word_list12{1}));
    end
end
end
    
% %for all subjects only extract the first 2 sesssions!
sess=uniquesession(~excludesessions);

% subjects sessions to edit mannual
if isub==38;
    sess=[1];   % 38 session 0 is practice 
elseif isub==37 
     sess=0;   
elseif isub==42
     sess=[0 1];
elseif isub==43 
     sess=[2 3 4 5]; % 43 has no session 0, and session 1 is practice  
elseif isub==36 
    sess=[0 1 2];      
elseif isub==34 
    sess=[0]; 
end

% sess=sess(1);
clear alignedStudyPair alignedTestProbe alignedRecallEvent 

%only extract first 3 sessions. 
if length(sess)<=2
    numsession=length(sess);
else
    numsession=2;
end
% numsession=length(sess);

% sort the order of the events    
for iss1=1:numsession
    SessionIDX=find(sess1==sess(iss1));
    for e = SessionIDX
        recEvent = recEvents(e);
        list = recEvent.list;
        serialpos = recEvent.serialpos;
        thisEv = filterStruct(sortedEvents.STUDY_PAIR,['session==' num2str(sess(iss1)) '&list==' num2str(list) '&serialpos==' num2str(serialpos)]);
        alignedStudyPair(e) = thisEv(1);
        thisEv = filterStruct(sortedEvents.TEST_PROBE,['session==' num2str(sess(iss1)) '&list==' num2str(list) '&serialpos==' num2str(serialpos)]);
        alignedTestProbe(e) = thisEv(1);
    end
end

tempRT=[alignedTestProbe.RT];


% Time Shift to the recall event.  %% eegoffset + RT or median RT. 
for iss1=1:numsession
    SessionIDX=find(sess1==sess(iss1));
    for e = SessionIDX
        alignedRecallEvent(e)=alignedTestProbe(e);
        if alignedTestProbe(e).RT > 0  %% eegoffset + RT
            alignedRecallEvent(e).eegoffset = alignedTestProbe(e).eegoffset + alignedTestProbe(e).RT;
        else % plus the median RT of the trials. 
            alignedRecallEvent(e).eegoffset = alignedTestProbe(e).eegoffset + round(nanmedian(tempRT(tempRT>0)));
        end
    end
end
        
alignedEvents.STUDY_PAIR = alignedStudyPair(~cellfun(@isempty,{alignedStudyPair.subject}));
alignedEvents.TEST_PROBE = alignedTestProbe(~cellfun(@isempty,{alignedStudyPair.subject}));
alignedEvents.RecallEvent = alignedRecallEvent(~cellfun(@isempty,{alignedStudyPair.subject}));

RTtemp=[alignedEvents.TEST_PROBE.RT];
ACCtemp=[alignedEvents.TEST_PROBE.correct];
SubjTable(isub).RTtemp=RTtemp;
SubjTable(isub).MinimumRT=nanmin(RTtemp(RTtemp>0 & RTtemp<10000));
SubjTable(isub).MeanRT=nanmean(RTtemp(RTtemp>0 & RTtemp<10000));
SubjTable(isub).MeanACC=nanmean(ACCtemp(ACCtemp>=0));

SubjTable(isub).alignedEvents=alignedEvents;
SubjTable(isub).PALTable=array2table([{alignedEvents.STUDY_PAIR.list}' {alignedEvents.STUDY_PAIR.serialpos}' ...
    {alignedEvents.TEST_PROBE.probepos}' {alignedEvents.STUDY_PAIR.study_1}' {alignedEvents.STUDY_PAIR.study_2}' {alignedEvents.TEST_PROBE.probe_word}' ...
    {alignedEvents.TEST_PROBE.expecting_word}' {alignedEvents.TEST_PROBE.resp_word}' {alignedEvents.TEST_PROBE.correct}' {alignedEvents.TEST_PROBE.RT}'],...
    'VariableNames',{'list' 'StudyPosition' 'ProbePosition' 'item1' 'item2' 'probe' 'expected' 'response' 'correct','RT'});
 
idx=all(cellfun(@isempty,SubjTable(isub).PALTable{:,:}),2);
SubjTable(isub).PALTable(idx,:)=[];

SubjTable(isub).subID= alignedEvents.TEST_PROBE(1).subject;
end

end

save SubjTable_palRAMword.mat SubjTable;

%% extract all unique probe words
% create stimuli by subject matrix. 

uniqueprobewords = textread('RAM_wordpool.txt', '%s', 'delimiter', '\n', 'whitespace', '');
uniqueresponsewords=uniqueprobewords;

for iw=1:length(uniqueprobewords)
    iw
    curword=uniqueprobewords(iw);
    
    uniqueprobe_sub(iw)=0;
    uniqueprobe_acc(iw)=0;
    
    for isub= allsub%1:length(allsub)
        uniqueprobe_subID_uniqueword(iw,isub) = 0;
        uniqueresponse_subID_uniqueword(iw,isub) = 0;
     
        if size(SubjTable(isub).PALTable,1)>1
            
            if sum(strcmp(SubjTable(isub).PALTable.expected,curword)) > 0 %if the response word match with the unique word
               uniqueresponse_subID_uniqueword(iw,isub) = 1;
               uniqueresponse_subID_uniqueword_acc(iw,isub) = sum([SubjTable(isub).PALTable.correct{1:end}]==1 & strcmp(SubjTable(isub).PALTable.expected,curword)');                         
            end
            
            if sum(ismember(SubjTable(isub).PALTable.probe,curword)) > 0 %if the probe word match with the unique word
                uniqueprobe_sub(iw)=uniqueprobe_sub(iw)+sum(ismember(SubjTable(isub).PALTable.probe,curword));
                uniqueprobe_acc(iw)=uniqueprobe_acc(iw)+sum([SubjTable(isub).PALTable.correct{1:end}]==1 & ismember(SubjTable(isub).PALTable.probe,curword)');
                
                uniqueprobe_subID_uniqueword(iw,isub) = 1;
                uniqueprobe_subID_uniqueword_acc(iw,isub) = sum([SubjTable(isub).PALTable.correct{1:end}]==1 & ismember(SubjTable(isub).PALTable.probe,curword)');          
            end
            
        end
    end
end


%% exclude words and subjects with too little samples.
clear PAL_Memo

k=1;excludeword=[];
for iw=1:length(uniqueprobewords)
    countallwordsubcombination(iw)=sum(uniqueprobe_subID_uniqueword(iw,:));
    if sum(uniqueprobe_subID_uniqueword(iw,:))<10
    excludeword(k)=iw;
    k=k+1;
    end
end
allwords=1:size(uniqueprobe_subID_uniqueword,1);
remainword=allwords(~ismember(allwords,excludeword));
uniqueprobe_subID_uniqueword_wordreduced=uniqueprobe_subID_uniqueword(remainword,:);
uniqueresponse_subID_uniqueword_wordreduced=uniqueresponse_subID_uniqueword(remainword,:);
uniqueprobewords_reducedwordsID=uniqueprobewords(remainword);

k=1; 
excludesub=[];
for isub=allsub %1:length(allsub)
    if sum(uniqueprobe_subID_uniqueword_wordreduced(:,isub))<1
       excludesub(k)=isub;
       k=k+1;
    end
end
remainsub=allsub(~ismember(allsub,excludesub));
uniqueprobe_subID_uniqueword_SUBwordreduced=uniqueprobe_subID_uniqueword_wordreduced(:,remainsub);
uniqueprobe_subID_uniqueword_SUBwordreduced_acc=uniqueprobe_subID_uniqueword_acc(remainword,remainsub)>0;
uniqueresponse_subID_uniqueword_SUBwordreduced = uniqueresponse_subID_uniqueword_wordreduced(:,remainsub);
uniqueresponse_subID_uniqueword_SUBwordreduced_acc=uniqueresponse_subID_uniqueword_acc(remainword,remainsub)>0;

[nword,nsub]=size(uniqueprobe_subID_uniqueword_SUBwordreduced);

nter=2000;
clear ProbeCorrelationFold  ResponseCorrelationFold  Response_probe_CorrelationFold
for iter=1:nter
  
   trIdx = randsample(1:nsub,round(nsub/2));
   teIdx = find(~ismember(1:nsub,trIdx));
    
%    X1=mean(uniqueprobe_subID_uniqueword_SUBwordreduced(:,trIdx),2) .* (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc(:,trIdx),2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced(:,trIdx),2));
%    X2=mean(uniqueprobe_subID_uniqueword_SUBwordreduced(:,teIdx),2) .* (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc(:,teIdx),2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced(:,teIdx),2));

   X1= (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc(:,trIdx),2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced(:,trIdx),2));
   X2= (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc(:,teIdx),2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced(:,teIdx),2));
   
   exclude=(isnan(X1)+isnan(X2))>0;      
   ProbeCorrelationFold(iter)= corr(X1(~exclude),X2(~exclude),'type','spearman');
   ProbeCorrelationFold(iter)=2*ProbeCorrelationFold(iter)/(ProbeCorrelationFold(iter)+1);
   
%    XX1=nanmean(uniqueresponse_subID_uniqueword_SUBwordreduced(:,trIdx),2).* (sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc(:,trIdx),2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced(:,trIdx),2));
%    XX2=nanmean(uniqueresponse_subID_uniqueword_SUBwordreduced(:,teIdx),2).* (sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc(:,teIdx),2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced(:,teIdx),2));
   XX1=(sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc(:,trIdx),2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced(:,trIdx),2));
   XX2=(sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc(:,teIdx),2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced(:,teIdx),2));
        
   exclude=(isnan(XX1)+isnan(XX2))>0;
   ResponseCorrelationFold(iter)= corr(XX1(~exclude),XX2(~exclude),'type','spearman');
   ResponseCorrelationFold(iter)=2*ResponseCorrelationFold(iter)/(ResponseCorrelationFold(iter)+1);

%    trIdx=randsample(1:nsub,nsub,true);
%    teIdx=trIdx;
%    XY1=mean(uniqueprobe_subID_uniqueword_SUBwordreduced(:,trIdx),2) .* (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc(:,trIdx),2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced(:,trIdx),2));   
%    XY2=mean(uniqueresponse_subID_uniqueword_SUBwordreduced(:,teIdx),2).* (sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc(:,teIdx),2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced(:,teIdx),2));
   XY1= (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc(:,trIdx),2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced(:,trIdx),2));   
   XY2= (sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc(:,teIdx),2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced(:,teIdx),2));

   exclude=(isnan(XY1)+isnan(XY2))>0;   
   Response_probe_CorrelationFold(iter)= corr(XY1(~exclude),XY2(~exclude),'type','spearman');
   Response_probe_CorrelationFold(iter)=2*Response_probe_CorrelationFold(iter)/(Response_probe_CorrelationFold(iter)+1);
   
end
p1 = sum(ProbeCorrelationFold<=0)/nter
p2 = sum(ResponseCorrelationFold<=0)/nter
p3 = sum(Response_probe_CorrelationFold<=0)/nter

% Probememorability = mean(uniqueprobe_subID_uniqueword_SUBwordreduced,2) .* (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc,2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced,2));
% Responsememorability = mean(uniqueresponse_subID_uniqueword_SUBwordreduced,2) .* (sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc,2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced,2));

% Here in this version, did not normalize based on the number of tested
% subjects
Probememorability = (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc,2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced,2));
Responsememorability = (sum(uniqueresponse_subID_uniqueword_SUBwordreduced_acc,2)./sum(uniqueresponse_subID_uniqueword_SUBwordreduced,2));

corr(Probememorability(~isnan(Responsememorability)),Responsememorability(~isnan(Responsememorability)),'type','spearman');

JustACC =  (sum(uniqueprobe_subID_uniqueword_SUBwordreduced_acc,2)./sum(uniqueprobe_subID_uniqueword_SUBwordreduced,2));

PAL_Memo.uniquewordsID_10=uniqueprobewords_reducedwordsID;
PAL_Memo.remainsub=remainsub;
PAL_Memo.remainword=remainword;
PAL_Memo.allsub=allsub;
PAL_Memo.uniqueprobe_SUBwordreduced=uniqueprobe_subID_uniqueword_SUBwordreduced;
PAL_Memo.uniqueprobe_SUBwordreduced_ACC=uniqueprobe_subID_uniqueword_SUBwordreduced_acc;
PAL_Memo.uniqueresponse_SUBwordreduced=uniqueresponse_subID_uniqueword_SUBwordreduced;
PAL_Memo.Probememorability=Probememorability;
PAL_Memo.Responsememorability=Responsememorability;
PAL_Memo.Probememorability_z=zscore(Probememorability);
PAL_Memo.Responsememorability_z=zscore(Responsememorability);

delete PAL_Memo_PALRAM.mat
save PAL_Memo_PALRAM.mat PAL_Memo;
%%
% generate plots
% stimuli by subjects mask

[sortMEM sortMEMInx] = sort(PAL_Memo.Responsememorability);
[a sortsub]=sort(mean(uniqueresponse_subID_uniqueword_acc));
sortsub = sortsub(ismember(sortsub,remainsub));

h=figure;
orient(h,'landscape')
fonsize=8;
subplot(2,3,1)
imagesc(uniqueresponse_subID_uniqueword(sortMEMInx,sortsub));
cm=colormap;
cm(1,:)=[1 1 1];
cm(end,:)=[0 0 0];
colormap(cm);
xlabel('Subjects (sorted by memory accuracy)');
ylabel('Target Words (sorted by memorability)');title('Stimuli by subjects occurrence');
set(gca,'FontSize',fonsize)

subplot(2,3,4)
bar(mean(uniqueresponse_subID_uniqueword(sortMEMInx,sortsub),2))
axis([0 301 0 1]);
xlabel('Target Words (sorted by memorability)');
ylabel('probability');title('Probabiliy of occurrence');
set(gca,'FontSize',fonsize)


subplot(2,3,2)
imagesc(uniqueresponse_subID_uniqueword_acc(sortMEMInx,sortsub)>0);
cm=colormap;
cm(1,:)=[1 1 1];
cm(end,:)=[0 0 0];
colormap(cm);
xlabel('Subjects (sorted by memory accuracy)');
ylabel('Target Words (sorted by memorability)');title('Success Recall');
set(gca,'FontSize',fonsize)


subplot(2,3,5)
bar(Responsememorability(sortMEMInx))
axis([0 301 0 1]);
xlabel('Target Words (sorted by memorability)');
ylabel('probability');title('Probabiliy of successful recall');
set(gca,'FontSize',fonsize)


subplot(2,3,3)
imagesc(uniqueprobe_subID_uniqueword_acc(sortMEMInx,sortsub)>0);
cm=colormap;
cm(1,:)=[1 1 1];
cm(end,:)=[0 0 0];
colormap(cm);
xlabel('Subjects (sorted by memory accuracy)');
ylabel('Probe Word (sorted by target memorability)');title('Success Recall');
set(gca,'FontSize',fonsize)

subplot(2,3,6)
bar(Probememorability(sortMEMInx))
axis([0 301 0 1]);
xlabel('Probe Words (sorted by target memorability)');
ylabel('probability');title('Probabiliy of successful recall');
set(gca,'FontSize',fonsize)


print('item_subj_analysis4.pdf','-dpdf','-bestfit');


%% across subject consistency plot
figure;
subplot(3,2,1);
scatter(X1,X2,'filled','MarkerFaceAlpha',0.2)
axis([0 0.6 0 0.6]);
xlabel('Joint accuracy in one half')
ylabel('Joint accuracy in another half')
title('Scatter plot of split-half across-subject consistency')
subplot(3,2,2);
hist(ProbeCorrelationFold,100);
title('Word as probe word')
xlabel('Spearman correlation')
ylabel('Count')
axis([-0.5 0.5 0 200]);
txt = ['p1 = ' num2str(p1)];
text(0.3,150,txt)
hold on;plot([0 0],[0 200],'-r')
 
subplot(3,2,3);
scatter(XX1,XX2,'filled','MarkerFaceAlpha',0.2)
axis([0 0.6 0 0.6]);
xlabel('Joint accuracy in one half')
ylabel('Joint accuracy in another half')
title('Scatter plot of split-half across-subject consistency')
subplot(3,2,4);
hist(ResponseCorrelationFold,100);
title('Word as to-be-retrieved word')
xlabel('Spearman correlation')
ylabel('Count');
 hold on;plot([0 0],[0 200],'-r')
axis([-0.5 0.5 0 200]);
txt = ['p2 = ' num2str(p2)];
text(0.3,150,txt)

subplot(3,2,5);
scatter(XY1,XY2,'filled','MarkerFaceAlpha',0.2)
axis([0 0.6 0 0.6]);
xlabel('Joint accuracy in one half')
ylabel('Joint accuracy in another half')
title('Scatter plot of split-half across-subject consistency')
subplot(3,2,6);
hist(Response_probe_CorrelationFold,100);
title('Response-Probe')
xlabel('Spearman correlation')
ylabel('Count')
axis([-0.5 0.5 0 200])
txt = ['p3 = ' num2str(p3)];
text(0.3,150,txt)
 hold on;plot([0 0],[0 200],'-r')

print('item_subj_split-half4.pdf','-dpdf','-bestfit');
