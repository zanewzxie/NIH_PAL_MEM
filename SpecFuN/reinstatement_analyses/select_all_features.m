function reinstatement_data = select_all_features(data)
    %Select all electrodes and freqs to build feature vectors
    %INPUT: data is a 4D matrix with dimensions trials X freqs X time X electrodes
    %OUTPUT: reinstatement_data is a 3D matrix with dimensions trials X time X features 
    
    num_els = size(data,4);
    num_freqs = size(data,2);
    reinstatement_data = zeros(size(data,1),size(data,3),num_els*num_freqs);
    i = 1;
    for f = 1:num_freqs
        for e = 1:num_els
            reinstatement_data(:,:,i) = data(:,f,:,e);
            i = i + 1;
        end
    end 
end