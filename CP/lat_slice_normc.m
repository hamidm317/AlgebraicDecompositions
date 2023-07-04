function Tn = lat_slice_normc(T)

    [~, ~, K] = size(T);
    Tn = zeros(size(T));
    
    for k = 1 : K
        
        Tn(:, :, k) = normc(T(:, :, k));
        
    end

end

