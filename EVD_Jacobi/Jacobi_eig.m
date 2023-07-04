function [jV, jD] = Jacobi_eig(A, off_var)

    n = length(A);
    V = eye(n);
    delta = off_var * off_norm(A);
    do_it = 1;

    while do_it

        err = off_norm(A)
        max_value = 0;

        for i = 1 : n
            for j = 1 : n
                
                if max_value < abs(A(i, j)) && i ~= j

                    max_value = abs(A(i, j));
                    p = max([i, j]);
                    q = min([i, j]);

                end
            
            end
        end

        [V, A] = EVmultiplier(V, A, p, q);

        if off_norm(A) < delta

            do_it = 0;

        end
    
    end

    jV = V;

    [s , Ind] = sort(diag(A), 'ascend');
            
    jD = diag(s);
    jV = jV(:,Ind);

end