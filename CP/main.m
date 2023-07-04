%% QUESTION 1
%% Only testing the Function of CP_ALS, using random initial values!
clear;

P = 5;
I = P;
J = 5;
K = 7;

A = normc(randn(I, P));
B = normc(randn(J, P));
C = normc(randn(K, P));

lambda = ones(P, 1);
L = lambda_generator(lambda, P);

U = {A, B, C};
T = tmprod(L, U, 1 : 3);

U1_0 = normc(randn(I, P));
U2_0 = normc(randn(J, P));
U3_0 = normc(randn(K, P));

[U1, U2, U3, lambda_] = CP_ALS_alg_M3(T, U1_0, U2_0, U3_0, P);

U_est = {U1, U2, U3};
L_est = lambda_generator(lambda_, P);

T_est = tmprod(L, U_est, 1 : 3);

%% Now, let's initialize loading matrices using SVD of mode-n unfolding!

U_0 = U_initializer(T, P);

[U1s, U2s, U3s, lambda_s] = CP_ALS_alg_M3(T, cell2mat(U_0(1)), cell2mat(U_0(2)), cell2mat(U_0(3)), P);
Us_est = {U1s, U2s, U3s};
%% TMSE calc!

random_init_error = TMSE(U, U_est); % / (norm(A, 'fro') ^ 2 + norm(B, 'fro') ^ 2 + norm(C, 'fro') ^ 2);
svd_init_error = TMSE(U, Us_est); % / (norm(A, 'fro') ^ 2 + norm(B, 'fro') ^ 2 + norm(C, 'fro') ^ 2);

%% QUESTION 2 - Part1/Uncorrelated Loading Matrices
%% first, initializings!

clear;

I = 5;
T = 5;
Q = 5;
P = 3;

N = 30;

U_1 = lat_slice_normc(randn(I, P, N));
U_2 = lat_slice_normc(randn(T, P, N));
U_3 = lat_slice_normc(randn(Q, P, N));

lambda = ones(P, 1);
L = lambda_generator(lambda, P);

Tensor = zeros(I, T, Q);

% for n = 1 : N
%     
%     U = {U_1(:, :, n), U_2(:, :, n), U_3(:, :, n)};
%     T_tmp = tmprod(L, U, 1 : 3);
%     
%     Tensor(:, :, :, n) = T_tmp;
%     
% end

Noise = randn(I, T, Q);

SNR = [0, 20, 40, 60];

%% Try it for different SNR values!

error_algorithm_snr_random = zeros(4, length(SNR));
error_algorithm_snr_svd = zeros(4, length(SNR));

tesnor_diff_error_algorithm_snr_random = zeros(4, length(SNR));
tesnor_diff_error_algorithm_snr_svd = zeros(4, length(SNR));

for i = 1 : length(SNR)
    
    snr = SNR(i);
    
    for n = 1 : N
        
        U = {U_1(:, :, n), U_2(:, :, n), U_3(:, :, n)};
        T_tmp = tmprod(L, U, 1 : 3);
        
        Noise_in = Noise * 10 ^ (-snr / 10) * norm(tens2mat(T_tmp, 1), 'fro');
        
        Noisy_T = T_tmp + Noise_in;
        
        % First, using random matrices:
        
        U_0 = {normc(randn(I, P)), normc(randn(T, P)), normc(randn(Q, P))};
        
        % Using my function
        [U1, U2, U3, ~] = CP_ALS_alg_M3(Noisy_T, cell2mat(U_0(1)), cell2mat(U_0(2)), cell2mat(U_0(3)), P);
        U_est_myFunc = {U1, U2, U3};
        
        U_est_als = cpd_als(Noisy_T, U_0);
        U_est_3sd = cpd3_sd(Noisy_T, U_0);
        U_est_minf = cpd_minf(Noisy_T, U_0);
        
        error_algorithm_snr_random(1, i) = error_algorithm_snr_random(1, i) + 1 / N * TMSE(U, U_est_myFunc);
        error_algorithm_snr_random(2, i) = error_algorithm_snr_random(2, i) + 1 / N * TMSE(U, U_est_als);
        error_algorithm_snr_random(3, i) = error_algorithm_snr_random(3, i) + 1 / N * TMSE(U, U_est_3sd);
        error_algorithm_snr_random(4, i) = error_algorithm_snr_random(4, i) + 1 / N * TMSE(U, U_est_minf);
        
        tesnor_diff_error_algorithm_snr_random(1, i) = tesnor_diff_error_algorithm_snr_random(1, i) + 1 / N * Tensor_error(U, U_est_myFunc);
        tesnor_diff_error_algorithm_snr_random(2, i) = tesnor_diff_error_algorithm_snr_random(2, i) + 1 / N * Tensor_error(U, U_est_als);
        tesnor_diff_error_algorithm_snr_random(3, i) = tesnor_diff_error_algorithm_snr_random(3, i) + 1 / N * Tensor_error(U, U_est_3sd);
        tesnor_diff_error_algorithm_snr_random(4, i) = tesnor_diff_error_algorithm_snr_random(4, i) + 1 / N * Tensor_error(U, U_est_minf);
        
        % Then, using HOSVD initialized matrices:
        
        U_0 = U_initializer(Noisy_T, P);
        
        % Using my function
        [U1, U2, U3, ~] = CP_ALS_alg_M3(Noisy_T, cell2mat(U_0(1)), cell2mat(U_0(2)), cell2mat(U_0(3)), P);
        U_est_myFunc = {U1, U2, U3};
        
        U_est_als = cpd_als(Noisy_T, U_0);
        U_est_3sd = cpd3_sd(Noisy_T, U_0);
        U_est_minf = cpd_minf(Noisy_T, U_0);
        
        error_algorithm_snr_svd(1, i) = error_algorithm_snr_svd(1, i) + 1 / N * TMSE(U, U_est_myFunc);
        error_algorithm_snr_svd(2, i) = error_algorithm_snr_svd(2, i) + 1 / N * TMSE(U, U_est_als);
        error_algorithm_snr_svd(3, i) = error_algorithm_snr_svd(3, i) + 1 / N * TMSE(U, U_est_3sd);
        error_algorithm_snr_svd(4, i) = error_algorithm_snr_svd(4, i) + 1 / N * TMSE(U, U_est_minf);
        
        tesnor_diff_error_algorithm_snr_svd(1, i) = tesnor_diff_error_algorithm_snr_svd(1, i) + 1 / N * Tensor_error(U, U_est_myFunc);
        tesnor_diff_error_algorithm_snr_svd(2, i) = tesnor_diff_error_algorithm_snr_svd(2, i) + 1 / N * Tensor_error(U, U_est_als);
        tesnor_diff_error_algorithm_snr_svd(3, i) = tesnor_diff_error_algorithm_snr_svd(3, i) + 1 / N * Tensor_error(U, U_est_3sd);
        tesnor_diff_error_algorithm_snr_svd(4, i) = tesnor_diff_error_algorithm_snr_svd(4, i) + 1 / N * Tensor_error(U, U_est_minf);
        
    end
    
end

%% QUESTION 2 - Part1/Correlated Loading Matrices
%% first, initializings!

clear;

I = 5;
T = 5;
Q = 5;
P = 3;

N = 30;

U_1 = lat_slice_normc(randn(I, P, N));
U_1(:, 2, :) = U_1(:, 1, :) + 0.5 * lat_slice_normc(randn(I, 1, N));
U_2 = lat_slice_normc(randn(T, P, N));
U_2(:, 2, :) = U_2(:, 1, :) + 0.5 * lat_slice_normc(randn(T, 1, N));
U_3 = lat_slice_normc(randn(Q, P, N));

lambda = ones(P, 1);
L = lambda_generator(lambda, P);

Tensor = zeros(I, T, Q);

% for n = 1 : N
%     
%     U = {U_1(:, :, n), U_2(:, :, n), U_3(:, :, n)};
%     T_tmp = tmprod(L, U, 1 : 3);
%     
%     Tensor(:, :, :, n) = T_tmp;
%     
% end

Noise = randn(I, T, Q);

SNR = [0, 20, 40, 60];

%% Try it for different SNR values!

error_algorithm_snr_random = zeros(4, length(SNR));
error_algorithm_snr_svd = zeros(4, length(SNR));

tesnor_diff_error_algorithm_snr_random = zeros(4, length(SNR));
tesnor_diff_error_algorithm_snr_svd = zeros(4, length(SNR));

for i = 1 : length(SNR)
    
    snr = SNR(i);
    
    for n = 1 : N
        
        U = {U_1(:, :, n), U_2(:, :, n), U_3(:, :, n)};
        T_tmp = tmprod(L, U, 1 : 3);
        
        Noise_in = Noise * 10 ^ (-snr / 10) * norm(tens2mat(T_tmp, 1), 'fro');
        
        Noisy_T = T_tmp + Noise_in;
        
        % First, using random matrices:
        
        U_0 = {normc(randn(I, P)), normc(randn(T, P)), normc(randn(Q, P))};
        
        % Using my function
        [U1, U2, U3, ~] = CP_ALS_alg_M3(Noisy_T, cell2mat(U_0(1)), cell2mat(U_0(2)), cell2mat(U_0(3)), P);
        U_est_myFunc = {U1, U2, U3};
        
        U_est_als = cpd_als(Noisy_T, U_0);
        U_est_3sd = cpd3_sd(Noisy_T, U_0);
        U_est_minf = cpd_minf(Noisy_T, U_0);
        
        error_algorithm_snr_random(1, i) = error_algorithm_snr_random(1, i) + 1 / N * TMSE(U, U_est_myFunc);
        error_algorithm_snr_random(2, i) = error_algorithm_snr_random(2, i) + 1 / N * TMSE(U, U_est_als);
        error_algorithm_snr_random(3, i) = error_algorithm_snr_random(3, i) + 1 / N * TMSE(U, U_est_3sd);
        error_algorithm_snr_random(4, i) = error_algorithm_snr_random(4, i) + 1 / N * TMSE(U, U_est_minf);
        
        tesnor_diff_error_algorithm_snr_random(1, i) = tesnor_diff_error_algorithm_snr_random(1, i) + 1 / N * Tensor_error(U, U_est_myFunc);
        tesnor_diff_error_algorithm_snr_random(2, i) = tesnor_diff_error_algorithm_snr_random(2, i) + 1 / N * Tensor_error(U, U_est_als);
        tesnor_diff_error_algorithm_snr_random(3, i) = tesnor_diff_error_algorithm_snr_random(3, i) + 1 / N * Tensor_error(U, U_est_3sd);
        tesnor_diff_error_algorithm_snr_random(4, i) = tesnor_diff_error_algorithm_snr_random(4, i) + 1 / N * Tensor_error(U, U_est_minf);
        
        % Then, using HOSVD initialized matrices:
        
        U_0 = U_initializer(Noisy_T, P);
        
        % Using my function
        [U1, U2, U3, ~] = CP_ALS_alg_M3(Noisy_T, cell2mat(U_0(1)), cell2mat(U_0(2)), cell2mat(U_0(3)), P);
        U_est_myFunc = {U1, U2, U3};
        
        U_est_als = cpd_als(Noisy_T, U_0);
        U_est_3sd = cpd3_sd(Noisy_T, U_0);
        U_est_minf = cpd_minf(Noisy_T, U_0);
        
        error_algorithm_snr_svd(1, i) = error_algorithm_snr_svd(1, i) + 1 / N * TMSE(U, U_est_myFunc);
        error_algorithm_snr_svd(2, i) = error_algorithm_snr_svd(2, i) + 1 / N * TMSE(U, U_est_als);
        error_algorithm_snr_svd(3, i) = error_algorithm_snr_svd(3, i) + 1 / N * TMSE(U, U_est_3sd);
        error_algorithm_snr_svd(4, i) = error_algorithm_snr_svd(4, i) + 1 / N * TMSE(U, U_est_minf);
        
        tesnor_diff_error_algorithm_snr_svd(1, i) = tesnor_diff_error_algorithm_snr_svd(1, i) + 1 / N * Tensor_error(U, U_est_myFunc);
        tesnor_diff_error_algorithm_snr_svd(2, i) = tesnor_diff_error_algorithm_snr_svd(2, i) + 1 / N * Tensor_error(U, U_est_als);
        tesnor_diff_error_algorithm_snr_svd(3, i) = tesnor_diff_error_algorithm_snr_svd(3, i) + 1 / N * Tensor_error(U, U_est_3sd);
        tesnor_diff_error_algorithm_snr_svd(4, i) = tesnor_diff_error_algorithm_snr_svd(4, i) + 1 / N * Tensor_error(U, U_est_minf);
        
    end
    
end