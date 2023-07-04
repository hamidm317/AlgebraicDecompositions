input_evd_mat = floor(rnd_Sym_matrix(4) * 100) / 10;

[jV, jD] = Jacobi_eig(input_evd_mat, 0.0001);

[jV_, jD_] = eig(input_evd_mat);

err_evd = norm(input_evd_mat - jV * jD * jV', 'fro');
err_evd_ = norm(input_evd_mat - jV_ * jD_ * jV_', 'fro');