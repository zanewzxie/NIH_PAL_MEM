  
function output = ReinstatementScore(encMat,recMat)
%Change the dimensions to be event X feature X time
%         encMat = SubjResults.Study_DataToDecode;
%         recMat = SubjResults.Recall_DataToDecode;

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
end