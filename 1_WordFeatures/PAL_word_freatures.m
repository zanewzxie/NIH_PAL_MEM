% Extract Crtical features to predict memorability behaviors


addpath(genpath('/Volumes/Zane/Matlab/Zane_Toolbox_V1/Corr_toolbox_v2'));

% Load calculated PAL mem scores
load('/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/NIH_PAL_MEM/PAL_Memo_PALRAM.mat')


%% =====1. Word Freq
% [Rank Freq_Word Pos Freq Disp] = textread('COCA.txt','%d %s %s %d %.2f');
% PAL_Word = textread('RAM_wordpool.txt','%s');
% PAL_Word =lower(PAL_Word);
% 
% for ipal=1:length(PAL_Word)
%     if ismember(PAL_Word(ipal),Freq_Word)
%             ipal
%         Freq_WordIndx(ipal)=find(ismember(Freq_Word,PAL_Word(ipal)));
%     end
% end
[Freq_Word LocalRank GloRank Freq] = textread('RAM_wordpool_freq_coca.txt','%s %d %d %d');

%% =====2. concreteness

[temRating ConcreteWord  ] = xlsread('Concreteness_ratings.xlsx');
[ConcreteRating ] = temRating(:,2);
PAL_Word = textread('RAM_wordpool.txt','%s');
PAL_Word =lower(PAL_Word);

for ipal=1:length(PAL_Word)
    ConcreteWordIndx(ipal)=nan;
    if ismember(PAL_Word(ipal),ConcreteWord)
        ipal
        ConcreteWordIndx(ipal)=find(ismember(ConcreteWord,PAL_Word(ipal)));
    end
    nanpal(ipal)=nan;
    if isnan(ConcreteWordIndx(ipal))
        nanpal(ipal)=ipal;
    end
end
Concrete_PAL.rating(find(~isnan(ConcreteWordIndx)))=ConcreteRating(ConcreteWordIndx(~isnan(ConcreteWordIndx)),:);
Concrete_PAL.rating(find(isnan(ConcreteWordIndx)))=mean(ConcreteRating(ConcreteWordIndx(~isnan(ConcreteWordIndx)),:));

Concrete_PAL.Word(find(~isnan(ConcreteWordIndx)))=ConcreteWord(ConcreteWordIndx(~isnan(ConcreteWordIndx)));
Concrete_PAL.Word(find(isnan(ConcreteWordIndx)))=PAL_Word(~isnan(nanpal));

%% == 3. GLOVE 

load Glove_PAL.mat
load('/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/NIH_PAL_MEM/PAL_ReducedBigTable.mat')

ReducedBigTable=table2struct(ReducedBigTable);

includedWord = PAL_Memo.Responsememorability ~= median(PAL_Memo.Responsememorability);
PALMemLabels = PAL_Memo.Responsememorability > median(PAL_Memo.Responsememorability);

T=array2table([Glove_PAL.Feature(includedWord,:) PALMemLabels(includedWord)]);

for iter=1:10
[trainedClassifier, validationAccuracy,AUC(iter)] = trainClassifier_svm_guassian(T);

T_sh=array2table([Glove_PAL.Feature(includedWord,:) Shuffle(PALMemLabels(includedWord))]);
[trainedClassifier, validationAccuracy,AUC_sh(iter)] = trainClassifier_svm_guassian(T_sh);
end


%%
% fixed beta 
IncludedW=1:300;
S_j_i=squareform(1-pdist(Glove_PAL.Feature, 'cosine'));

for jj=1:300 %
    P_Q(jj) =  sum((S_j_i(jj,IncludedW~=jj).^1));  % Q,I; a Q cue for all I. 
end

for ii=1:300
    
    cueInx=IncludedW~=ii;
    
    P_R_i(ii)= sum((S_j_i(IncludedW~=ii,ii)'.^1./P_Q(cueInx))* 1/300); % Q,I; I for all Qs
end

% P_R_i=P_R./299;
robust_correlation(P_R_i,PAL_Memo.Responsememorability')

r_test=corr(P_R_i',PAL_Memo.Responsememorability,'type','spearman');

%% Read the actual big data table

allexpect=lower({ReducedBigTable.expected});
allcue=lower({ReducedBigTable.probe});

clear TempC TempS
uniquetarget=unique(allexpect);

for iresp=1:length(uniquetarget)

     CurrentCue=allcue(ismember(allexpect,uniquetarget(iresp)));   
     cueInx=find(ismember(Glove_PAL.Word,CurrentCue));
     targetInx=find(ismember(Glove_PAL.Word,uniquetarget(iresp)));
     for ic=1:length(cueInx)
        TempS{iresp}(ic)=1-pdist(Glove_PAL.Feature([cueInx(ic),targetInx],:), 'cosine');
     end
    
     meanstrength(iresp)= sum((TempS{iresp}./P_Q(cueInx))*1/300);
%      meanstrength(iresp)= sum((TempC{iresp}./P_Q(targetInx))*1/300);
end

robust_correlation(meanstrength,P_R_i)
robust_correlation(meanstrength,PAL_Memo.Responsememorability')
robust_correlation(P_R_i,PAL_Memo.Responsememorability')


% caclculate memorability based on list
Mem_Full=P_R_i;
Mem_Presented=meanstrength;

save ModelPredictedMem.mat Mem_Full Mem_Presented

%% == 4. GLOVE Online study
clear targetperfdata

load('/Volumes/Zane/NIH_HPC/NIH_PAL_Mem/Scripts/1_WordFeatures/Onlinedata/onlinememorydata-042619.mat')

targetperfdata.DART

uniqueprobewords = textread('RAM_wordpool.txt', '%s', 'delimiter', '\n', 'whitespace', '');
uniqueresponsewords=uniqueprobewords;

fnames = fieldnames(targetperfdata);
clear AMT_MemScore

for iw=1:length(uniqueresponsewords)
    
    curInX=find(ismember(fnames,uniqueresponsewords(iw)));
    
    fieldval = getfield(targetperfdata,fnames{curInX});  
    AMT_MemScore(iw)=fieldval.hr;

%     
%     fieldval2 = getfield(intrusioncounts,fnames{curInX});
%     AMT_IntrusionCount(iw)=fieldval2;
end

robust_correlation(AMT_MemScore,PAL_Memo.Responsememorability');
robust_correlation(AMT_MemScore,meanstrength);
robust_correlation(AMT_MemScore,P_R_i);
robust_correlation(AMT_MemScore,log(AMT_IntrusionCount));



%% list based analysis
WordtoAnalyzeT=struct2table(WordtoAnalyze);
WordtoAnalyzeT2=WordtoAnalyzeT(memratio>0,:);

for isubtem=1:height(WordtoAnalyzeT2)
    isubtem
Curlist=cell2mat(WordtoAnalyzeT2.SubjTable{isubtem}.PALTable.list);
sessions=[[WordtoAnalyzeT2.SubjTable{isubtem}.alignedEvents.RecallEvent.session]];
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

        
clear ListTargetMem
targetwordsall=[];
for ilist=1:length(uniquelist)
    listmask=Curlistupdate==uniquelist(ilist);
    allexpect=lower(WordtoAnalyzeT2.SubjTable{isubtem}.PALTable.expected(find(listmask)));
    allcue=[allcue lower(WordtoAnalyzeT2.SubjTable{isubtem}.PALTable.probe(find(listmask)))];

    if isubtem==23 & uniquelist(ilist)==400
        uniquetarget = allexpect(2:7);
        currcue=[lower(WordtoAnalyzeT2.SubjTable{isubtem}.PALTable.probe(find(listmask)))];
        allcue= [allcue currcue(2:7)];
        allcue= [uniquetarget;allcue];
        inx=find(Curlistupdate==uniquelist(ilist));
        listmask(inx(1))=0;
    else
        uniquetarget = allexpect;
        allcue= [uniquetarget;allcue];
    end
    
clear PQ TempC
for icue=1:length(allcue)
     CurrentCue=allcue(ismember(allcue,allcue(icue)));   
     cueInx=find(ismember(Glove_PAL.Word,CurrentCue));
     targets=allcue(~ismember(allcue,allcue(icue)));  
     targetInx=find(ismember(Glove_PAL.Word,targets));
     for ic=1:length(targetInx)
        TempC{icue}(ic)=1-pdist(Glove_PAL.Feature([cueInx,targetInx(ic)],:), 'cosine');
     end     
     PQ(icue)=sum(TempC{icue});
end     

clear meanstrength  TempS
for iresp=1:length(uniquetarget)
     CurrentCue=allcue(~ismember(allcue,uniquetarget(iresp)));   
     cueInx=find(ismember(Glove_PAL.Word,CurrentCue));
     targetInx=find(ismember(Glove_PAL.Word,uniquetarget(iresp)));
     for ic=1:length(cueInx)
        TempS{iresp}(ic)=1-pdist(Glove_PAL.Feature([cueInx(ic),targetInx],:), 'cosine');
     end
     localcueidx=find(ismember(allcue,CurrentCue));
     meanstrength(iresp)= sum((TempS{iresp}./P_Q(localcueidx))*1/11);
%      meanstrength(iresp)= sum((TempC{iresp}./P_Q(targetInx))*1/300);
end
ListTargetMem(listmask)=meanstrength;
end
        
allexpect=lower(WordtoAnalyzeT2.SubjTable{isubtem}.PALTable.expected);

% sorted by glove word.
Glove_PAL.Word;

AllWordinx=find(ismember(Glove_PAL.Word,allexpect));

sub_by_word_ListMem(isubtem,:)=nan(300,1);
for iw=1:length(AllWordinx)
    extractfromlist= find(ismember(allexpect,Glove_PAL.Word(AllWordinx(iw))));
    sub_by_word_ListMem(isubtem,AllWordinx(iw)) =  nanmean(ListTargetMem(extractfromlist));
end

end


