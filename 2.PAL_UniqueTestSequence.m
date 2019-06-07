
% find unique sequence of PAL association.

allsub=[26:66];
% first make a big table.
BIGTABLE=[];
for isub=allsub
    if size(SubjTable(isub).PALTable,1)>1
        BIGTABLE=vertcat(BIGTABLE,SubjTable(isub).PALTable);
    end
end

uniqueprobewords = textread('RAM_wordpool.txt', '%s', 'delimiter', '\n', 'whitespace', '');
uniqueresponsewords=uniqueprobewords;

BIGTABLE=BIGTABLE(ismember(BIGTABLE.expected,uniqueresponsewords) & ismember(BIGTABLE.probe,uniqueresponsewords),:);

[uniquetable idx]= unique(BIGTABLE(:,6:7), 'rows');

ReducedBigTable=BIGTABLE(sort(idx),:);
% ReducedBigTable = sortrows(ReducedBigTable,'expected')
save PAL_ReducedBigTable.mat ReducedBigTable


% verify memorability
allacc=[ReducedBigTable.correct{1:end}];allacc=allacc>0;
for iw=1:length(uniqueresponsewords)
    uniquewordMemorabilityverify(iw)=mean(allacc(ismember(ReducedBigTable.expected,uniqueresponsewords(iw))));
end
load('PAL_Memo_PALRAM.mat')
figure;scatter(PAL_Memo.Responsememorability,uniquewordMemorabilityverify)



