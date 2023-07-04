function [A, B, C, lambda_] = CP_ALS_alg_M3(T, A0, B0, C0, P)

    B = B0;
    C = C0;
    A = A0;
    do_it = 1;
    iteration = 0;
    max_iter = 1000;
    
    L = lambda_generator(ones(P, 1), P);
    
    while do_it && iteration < max_iter
        
        V = (B' * B) .* (C' * C);
        A = tens2mat(T, 1) * kr(C, B) * pinv(V);
        
        [A_lambda, A] = mat_norm(A);
        
        V = (A' * A) .* (C' * C);
        B = tens2mat(T, 2) * kr(C, A) * pinv(V);
        
        [B_lambda, B] = mat_norm(B);
        
        V = (A' * A) .* (B' * B);
        C = tens2mat(T, 3) * kr(B, A) * pinv(V);
        
        [C_lambda, C] = mat_norm(C);
        
        iteration = iteration + 1;
        
        U = {A, B, C};
        T_est = tmprod(L, U, 1 : 3);

        error = norm(tens2mat(T, 1) - tens2mat(T_est, 1), 'fro');
        
        if error < 0.00001 * norm(tens2mat(T, 1), 'fro')
            
            do_it = 0;
            
        end
        
    end
    
    lambda_ = A_lambda .* B_lambda .* C_lambda;
    error = norm(tens2mat(T, 1) - tens2mat(T_est, 1), 'fro')
    

end

