function reinstatement_data = select_features_by_freq(data,f)
    %Select features for a particular freq   
    num_els = size(data,4);
    reinstatement_data = zeros(size(data,1),size(data,3),num_els);
    i = 1;
    for e = 1:num_els
        reinstatement_data(:,:,i) = data(:,f,:,e);
        i = i + 1;
    end
end