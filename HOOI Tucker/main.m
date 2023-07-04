%% Question 1 - HOOI implementation
% Initialize T!

R = [4, 5, 6];
T_R = [7, 8, 9];

G = randn(R);
U = {normc(randn(T_R(1), R(1))), normc(randn(T_R(2), R(2))), normc(randn(T_R(3), R(3)))};

% T = tmprod(G, U, 1 : 3);
T = ttm(G, U).data;


[G_est, U_est] = HOOI_Tucker(T, R);

% T_est = tmprod(G_est, U_est, 1 : 3);
T_est = ttm(G_est, U_est).data;

%%

TMSE_HOOI_Error = TMSE(U, U_est);
Tensor_Diff_Error = norm(tens2mat(T - T_est, 1), 'fro');