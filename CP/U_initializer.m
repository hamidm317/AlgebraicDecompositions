function U_0 = U_initializer(T, P)

    [u1, ~, ~] = svd(tens2mat(T, 1));
    U1_0 = u1(:, 1 : P);
    
    [u2, ~, ~] = svd(tens2mat(T, 2));
    U2_0 = u2(:, 1 : P);
    
    [u3, ~, ~] = svd(tens2mat(T, 3));
    U3_0 = u3(:, 1 : P);
    
    U_0 = {U1_0, U2_0, U3_0};

end

