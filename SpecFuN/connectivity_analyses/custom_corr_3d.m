function sim = custom_corr_3d(a,b)
    %This function computes the pearson correlation along the first dimension
    %of two 3d matrices of the same size
    ma = mean(a);
    mb = mean(b);
    
    a = bsxfun(@minus,a,ma);
    b = bsxfun(@minus,b,mb);
    
    %{
    ma = repmat(ma,[1 1 size(a,3)]);
    ma = permute(ma,[3 1 2]);
    a = a - ma;
    
    mb = repmat(mb,[1 1 size(b,1)]);
    mb = permute(mb,[3 1 2]);
    b = b - mb;
    %}
    
    s = sum(a.*b);
    norma = sqrt(sum(a.^2));
    normb = sqrt(sum(b.^2));
    sim = s./(norma.*normb);
end