function [lambda, M_n] = mat_norm(M)

    [~, num_col] = size(M);
    
    for i = 1 : num_col
        
        lambda(i) = norm(M(:, i));
        M_n(:, i) = M(:, i) / lambda(i);
        
    end

end

