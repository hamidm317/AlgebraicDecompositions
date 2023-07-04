
clear;
%% Part A - ALS algorithm

m = 6;
n = 4;
j = 3;

iB = rand(m, j);
iC = rand(j, n);

A = iB * iC;

B0 = rand(m, j);


[B, C, e] = NNMF_ALS(A, j, B0, 1000, 0.0001);

%% Part B - Multiplicative algorithm

m = 6;
n = 4;
j = 3;

iB = rand(m, j);
iC = rand(j, n);

A = iB * iC;

B0 = rand(m, j);
C0 = rand(j, n);


[B, C, e] = NNMF_Mul(A, B0, C0, 10000, 0.00001);

%% Part C
%% initialization

clear;

m = 6;
n = 4;
j = 3;
SNR = [-10, 0, 10, 30, 50];

iB = rand(m, j);
iC = rand(j, n);
iE = rand(m, n);

alpha = norm(iB * iC, 'fro') ./ (10 .^ (SNR ./ 20) .* norm(iE, 'fro'));

for i = 1 : length(SNR)
    
    A(i, :, :) = iB * iC + alpha(i) * iE;
    
end

B0 = rand(m, j);
C0 = rand(j, n);

%% Matlab nnmf functions!

% [B_ALS, C_ALS] = nnmf(A, j, 'algorithm', 'als', 'w0', B0, 'h0', C0);
% [B_Mul, C_Mul] = nnmf(A, j, 'algorithm', 'mult', 'w0', B0, 'h0', C0);

%% Error of each algorithm for each SNR value

Error_Values = zeros(4, length(SNR));
num_iter = 10;

for i = 1 : length(SNR)
        
    A_ = reshape(A(i, :, :), m, n);

    [B_myALS, C_myALS, ~] = NNMF_ALS(A_, B0, 10000, 0.0000001);
    [B_myMul, C_myMul, ~] = NNMF_Mul(A_, B0, C0, 10000, 0.000001);

    [B_ALS, C_ALS] = nnmf(A_, j, 'algorithm', 'als', 'w0', B0, 'h0', C0);
    [B_Mul, C_Mul] = nnmf(A_, j, 'algorithm', 'mult', 'w0', B0, 'h0', C0);

    Error_Values(1, i) = norm(A_ - B_myALS * C_myALS, 'fro') / norm(A_, 'fro');
    Error_Values(2, i) = norm(A_ - B_myMul * C_myMul, 'fro') / norm(A_, 'fro');

    Error_Values(3, i) = norm(A_ - B_ALS * C_ALS, 'fro') / norm(A_, 'fro');
    Error_Values(4, i) = norm(A_ - B_Mul * C_Mul, 'fro') / norm(A_, 'fro');
    
end

%% plotting error values

figure()
stem(SNR, Error_Values(1, :))
xlim([-20 60])
title('Error probabilities for my ALS NNMF by SNR')
xlabel('SNR (dB)')
ylabel('Error Probability')

figure()
stem(SNR, Error_Values(2, :))
xlim([-20 60])
title('Error probabilities for my Multiplicative NNMF by SNR')
xlabel('SNR (dB)')
ylabel('Error Probability')

figure()
stem(SNR, Error_Values(3, :))
xlim([-20 60])
title('Error probabilities for ALS NNMF by SNR')
xlabel('SNR (dB)')
ylabel('Error Probability')

figure()
stem(SNR, Error_Values(4, :))
xlim([-20 60])
title('Error probabilities for Multiplicative NNMF by SNR')
xlabel('SNR (dB)')
ylabel('Error Probability')
%% Error of each algorithm for each SNR value

Error_Values = zeros(3, 4, length(SNR));
num_iter = 10;

for j = 2 : 4
    
    for iter = 1 : num_iter

        iB = rand(m, j);
        iC = rand(j, n);
        iE = rand(m, n);

        alpha = norm(iB * iC, 'fro') ./ (10 .^ (SNR ./ 20) .* norm(iE, 'fro'));

        for i = 1 : length(SNR)

            A(i, :, :) = iB * iC + alpha(i) * iE;

        end

        B0 = rand(m, j);
        C0 = rand(j, n);

        for i = 1 : length(SNR)

            A_ = reshape(A(i, :, :), m, n);

            [B_myALS, C_myALS, ~] = NNMF_ALS(A_, B0, 10000, 0.0000001);
            [B_myMul, C_myMul, ~] = NNMF_Mul(A_, B0, C0, 10000, 0.000001);

            [B_ALS, C_ALS] = nnmf(A_, j, 'algorithm', 'als', 'w0', B0, 'h0', C0);
            [B_Mul, C_Mul] = nnmf(A_, j, 'algorithm', 'mult', 'w0', B0, 'h0', C0);

            Error_Values(j - 1, 1, i) = Error_Values(j - 1, 1, i) + norm(A_ - B_myALS * C_myALS, 'fro') / norm(A_, 'fro');
            Error_Values(j - 1, 2, i) = Error_Values(j - 1, 2, i) + norm(A_ - B_myMul * C_myMul, 'fro') / norm(A_, 'fro');

            Error_Values(j - 1, 3, i) = Error_Values(j - 1, 3, i) + norm(A_ - B_ALS * C_ALS, 'fro') / norm(A_, 'fro');
            Error_Values(j - 1, 4, i) = Error_Values(j - 1, 4, i) + norm(A_ - B_Mul * C_Mul, 'fro') / norm(A_, 'fro');

        end

    end
    
end

Error_Values = Error_Values / num_iter;

%% plotting again!

figure();
hold on

plot(SNR, reshape(Error_Values(1, 1, :), 1, 5))
plot(SNR, reshape(Error_Values(2, 1, :), 1, 5))
plot(SNR, reshape(Error_Values(3, 1, :), 1, 5))

legend('j = 2', 'j = 3', 'j = 4')
title('Error probabilities for my ALS NNMF by SNR for different number of factors')
xlabel('SNR (dB)')
ylabel('Error Probability')

hold off

figure();
hold on

plot(SNR, reshape(Error_Values(1, 2, :), 1, 5))
plot(SNR, reshape(Error_Values(2, 2, :), 1, 5))
plot(SNR, reshape(Error_Values(3, 2, :), 1, 5))

legend('j = 2', 'j = 3', 'j = 4')
title('Error probabilities for my Mul NNMF by SNR for different number of factors')
xlabel('SNR (dB)')
ylabel('Error Probability')

hold off

figure();
hold on

plot(SNR, reshape(Error_Values(1, 3, :), 1, 5))
plot(SNR, reshape(Error_Values(2, 3, :), 1, 5))
plot(SNR, reshape(Error_Values(3, 3, :), 1, 5))

legend('j = 2', 'j = 3', 'j = 4')
title('Error probabilities for ALS NNMF by SNR for different number of factors')
xlabel('SNR (dB)')
ylabel('Error Probability')

hold off

figure();
hold on

plot(SNR, reshape(Error_Values(1, 4, :), 1, 5))
plot(SNR, reshape(Error_Values(2, 4, :), 1, 5))
plot(SNR, reshape(Error_Values(3, 4, :), 1, 5))

legend('j = 2', 'j = 3', 'j = 4')
title('Error probabilities for Mul NNMF by SNR for different number of factors')
xlabel('SNR (dB)')
ylabel('Error Probability')

hold off