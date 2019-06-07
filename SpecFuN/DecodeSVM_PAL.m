function [DecodeACC, DecodeAUC,DecodeACC_sh, DecodeAUC_sh] = DecodeSVM_PAL(DataToDecode, DecodeInx, alltimetasktimeinx, Labels)



    
    parfor itime=1:length(alltimetasktimeinx)
        
        Rtempdata=squeeze(DataToDecode(:,:,alltimetasktimeinx(itime))); 

        T=array2table([Rtempdata(DecodeInx,:) Labels']);
        [~,DecodeACC(itime),DecodeAUC(itime)] = trainClassifier_svm_guassian(T);
        
        T_sh=array2table([Rtempdata(DecodeInx,:) Shuffle(Labels)']);
        [~,DecodeACC_sh(itime),DecodeAUC_sh(itime)] = trainClassifier_svm_guassian(T_sh);        
    end
    
 
end