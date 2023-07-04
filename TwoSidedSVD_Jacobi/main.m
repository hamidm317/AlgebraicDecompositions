m2 = 8;
n2 = 5;

A2 = floor(rand(m2, n2, 1) * 100) / 10;

[U2, S2, V2] = Jacobi_svd_2sided(A2, 0.01);

[Um, Sm, Vm] = svd(A2);

err_2 = norm(A2 - U2 * S2 * V2', 'fro');