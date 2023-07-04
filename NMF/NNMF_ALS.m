function [B, C, e] = NNMF_ALS(A, B0, max_iter, CoCr)

    % A = BC + E
    % A = n by m, B = n by j, C = j by m
    
    B = B0;
    do_it = 1;
    iterations = 0;
    e = [];
    
    while do_it
        
        C = pinv(B' * B) * B' * A;
        C = anti_negative(C);
        
        B = (pinv(C * C') * C * A')';
        B = anti_negative(B);
        
        e = [e norm(A - B * C, 'fro') / norm(A, 'fro')];
        
        if e(iterations + 1) < CoCr || iterations == max_iter
            
            do_it = 0;
            
        end
        
        iterations = iterations + 1;
        
    end

end