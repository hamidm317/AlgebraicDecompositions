function error = Tensor_error(U, U_est)

    error = norm(tens2mat(cpdgen(U) - cpdgen(U_est), 1), 'fro');

end

