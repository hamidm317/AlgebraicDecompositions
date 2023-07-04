function [B, C, e] = NNMF_Mul(A, B0, C0, max_iter, CoCr)

    do_it = 1;
    e = [];
    iteration = 1;
    
    B = B0;
    C = C0;
    
    while do_it
        
        C = C .* (B' * A) ./ (B' * B * C + 10 ^ -9);
        B = B .* (A * C') ./ (B * C * C' + 10 ^ -9);
        
        e = [e norm(A - B * C, 'fro')];
        
        if e(iteration) < CoCr || iteration > max_iter
            
            do_it = 0;
            
        end
        
        iteration = iteration + 1;
        
    end

end

