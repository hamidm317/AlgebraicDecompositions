function tr_mat = truncating_(mat, trun_Num)

    [u, s, v] = svd(mat);
    
    tr_mat = u(:, 1 : trun_Num) * s(1 : trun_Num, 1 : trun_Num) * v(:, 1 : trun_Num)';

end

