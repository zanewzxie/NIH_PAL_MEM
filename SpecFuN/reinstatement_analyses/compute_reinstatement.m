function [output, r_each] = compute_reinstatement(encMat,recMat)
    %This function computes the similarity matrices between encoding
    %vectors and recall vectors.
    %INPUTS
        %encMat: an numEvents X numEncWin X numElements matrix of z-scored
        %power values
        %recMat: an numEvents X numRecWin X numElements matrix of z-scored
        %power values
        %encMat and recMat must have the same numEvents and numElements
        %(1st and 3rd dimensions)
    %OUTPUTS
        %output: numEvents X numEncWin X numRecWin matrix of similarity
        %values
        %r_each: numElements X numEncWin X numRecWin  matrix of similarity
        %values

    numEncWin = size(encMat,2);
    numRecWin = size(recMat,2);
    numEvents = size(encMat,1);
    numElements = size(encMat,3);
    output = zeros(numEvents,numEncWin,numRecWin);
    r_each = zeros(numElements,numEncWin,numRecWin);
    
    for feature = 1:numElements
        encData = encMat(:,:,feature);
        recData = recMat(:,:,feature);
        s = (encData')*recData;
        normEnc = sqrt(sum(encData.^2));
        normRec = sqrt(sum(recData.^2));
        norms = (normEnc')*normRec;
        r_each(feature,:,:) = s./norms;
    end
    
    %Change the dimensions to be feature X time X event
    encMat = permute(encMat,[3 2 1]);
    recMat = permute(recMat,[3 2 1]);
    
    for event = 1:numEvents
        encData = encMat(:,:,event);
        recData = recMat(:,:,event);
        s = (encData')*recData;
        normEnc = sqrt(sum(encData.^2));
        normRec = sqrt(sum(recData.^2));
        norms = (normEnc')*normRec;
        output(event,:,:) = s./norms;
    end
    
end