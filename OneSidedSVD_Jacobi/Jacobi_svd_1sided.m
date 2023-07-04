function [jU1, jS1, jV1] = Jacobi_svd_1sided(A, off_var)
    
    [m, n] = size(A);
    V = eye(n);

    do_it = 1;
    iter = 0;
    delta = off_var ^ 2 * off_norm(A);

    while do_it

        B = A' * A;

        iter = iter + 1
        err = off_norm(B)

        max_value = 0;

        for i = 1 : length(B)
            for j = 1 : length(B)

                if max_value < abs(B(i, j)) && i ~= j

                    max_value = abs(B(i, j));
                    p = max([i, j]);
                    q = min([i, j]);

                end
            end
        end

        [V, A] = JKmultiplier_(V, A, B, p, q);
        
        if off_norm(B) < delta

            do_it = 0;

        end

    end

    % 
    
    jS1 = zeros(size(A));

    for i = 1 : min(size(jS1))

        jS1(i, i) = norm(A(:, i));
    
    end
    
    jV1 = V;

    jU1 = A * pinv(jS1);

end