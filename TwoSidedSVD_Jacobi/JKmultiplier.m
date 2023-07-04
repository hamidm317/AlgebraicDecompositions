function [U_new, V_new, A_new] = JKmultiplier(U, V, A, p, q)

    [m, n] = size(A);
    [J, K, d] = JKestim(A, p, q);

    U_ = eye(m);
    V_ = eye(n);
    
    U_(q, q) = J(1, 1);
    U_(p, p) = J(2, 2);
    U_(p, q) = J(2, 1);
    U_(q, p) = J(1, 2);
    
    V_(q, q) = K(1, 1);
    V_(p, p) = K(2, 2);
    V_(p, q) = K(2, 1);
    V_(q, p) = K(1, 2);
    
    U_new = U * U_;
    V_new = V * V_;
    
    A_new = U_' * A * V_;
    
end

