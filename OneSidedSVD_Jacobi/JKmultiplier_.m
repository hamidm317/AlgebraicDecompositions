function [V_new, A_new] = JKmultiplier_(V, A, B, p, q)

    [m, n] = size(B);
    [~, K, d] = JKestim(B, p, q);

    V_ = eye(n);
    
    V_(q, q) = K(1, 1);
    V_(p, p) = K(2, 2);
    V_(p, q) = K(2, 1);
    V_(q, p) = K(1, 2);

    V_new = V * V_;
    
    A_new = A * V_;
    
end

