function [err, err_clean] = factorization_test(algHandle, algParameters)
% Returns error for the given algorithm with default
% parameters, and opts.L1 = smoothL11
%
% Uses B * Phi where Phi is a subset of the features of
% Gaussian generated X. i.e. Phi = X and B = eye, B(:, some_ind) = 0

initialize_testing_paths
global sigma_smooth_L11;
sigma_smooth_L11 = 1e-2;

rand('state', 1);
randn('state', 1);

t = 100;   % Number of samples
n = 10;   % Dimension of X
r = 10;   % Dimension of Y
dim = 5;  % Rank of combined matrix
sigma = 0.05; % Noise on data

% Gaussian generated X. i.e. Phi = X and B = eye, B(:, some_ind) = 0
X = randn(n, t);
Y = randn(r, t);
Phi = [X; Y];
B = eye(n+r);
removedind = randperm(n+r);
removedind = removedind(1:(n+r-dim));
B(:, removedind) = 0;

%First sanity check on clean data
sigma = 0.1;
algParameters.L1 = @smooth_L11;
algParameters.L2 = @smooth_L11;
algParameters.reg_loss = @L21_loss;
algParameters.reg_wgt = 0;

Zclean = B*Phi;
Z = Zclean;

[Xrecon, Yrecon] = algHandle(Z(1:n, :), Z((n+1):(n+r), :), algParameters);
err_clean = factorization_err(Zclean, [Xrecon; Yrecon]);

Z = Zclean + sigma*randn(n+r, t);
[Xrecon, Yrecon] = algHandle(Z(1:n, :), Z((n+1):(n+r), :), algParameters);
err = factorization_err(Zclean, [Xrecon; Yrecon]);

end

