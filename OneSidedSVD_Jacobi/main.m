m1 = 9;
n1 = 7;

A1 = floor(rand(m1, n1, 1) * 100) / 10;


[U1, S1, V1] = Jacobi_svd_1sided(A1, 0.01);


[Um, Sm, Vm] = svd(A1);

err_1 = norm(A1 - U1 * S1 * V1', 'fro');