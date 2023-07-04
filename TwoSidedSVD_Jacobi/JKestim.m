function [J, K, D] = JKestim(A, p, q)
    
    x = [A(q, q), A(q, p); A(p, q), A(p, p)];
    
    phi = atan((x(1, 2) - x(2, 1)) / (x(1, 1) + x(2, 2)));
    
    x_sym = rot_mat(phi)' * x;
    
%     theta = 0.5 * atan(2 * x_sym(1, 2) / (x_sym(1, 1) - x_sym(2, 2)));

    taw = (x_sym(2, 2) - x_sym(1, 1)) / (2 * x_sym(1, 2));
    tangent_value = min([-taw + sqrt(taw ^ 2 + 1), -taw - sqrt(taw ^ 2 + 1)]);
    
    theta = atan(tangent_value);
    
    D = rot_mat(theta)' * x_sym * rot_mat(theta);
    
    J = rot_mat(phi) * rot_mat(theta);
    K = rot_mat(theta);

end

