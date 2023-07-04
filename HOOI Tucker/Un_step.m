function A_n = Un_step(Z, Rn, mode)

    [u, ~, ~] = svd(tens2mat(Z, mode));
    A_n = u(:, 1 : Rn);

end

