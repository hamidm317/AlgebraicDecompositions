function R = rot_mat(theta)

    R = zeros(2, 2);
    
    R(1, 1) = cos(theta);
    R(2, 2) = R(1, 1);
    
    R(1, 2) = sin(theta);
    R(2, 1) = -1 * R(1, 2);

end

