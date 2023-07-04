function [jU2, jS2, jV2] = Jacobi_svd_2sided(A, off_var)
    
    [U, V, A_] = prepare(A);
    
    [m, n] = size(A_);
    
    do_it = 1;

    delta = off_var * off_norm(A_)
    iter = 0;

    while do_it
        max_val = 0;
        iter = iter + 1
        off_norm(A_)

        for i = 1 : m
            for j = 1 : n

                if abs(A_(i, j)) > max_val && i ~= j

                    max_val = abs(A_(i, j));
                    p = max([i, j]);
                    q = min([i, j]);

                end

            end
        end

        [U, V, A_] = JKmultiplier(U, V, A_, p, q);

        if off_norm(A_) < delta

            do_it = 0;

        end

    end

    jU2 = U;
    jV2 = V;
    jS2 = A_;

    for i = 1 : min(size(A_))

        if jS2(i, i) < 0
            
            jS2(i, i) = -jS2(i, i);
            jU2(:, i) = -jU2(:, i);

        end
    end

end