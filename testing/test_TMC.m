%TEST
clear;clc;
initialize_testing_paths

t = 100; 
x_dim = 20;
rnk = 4;
y_dim = 10;
per = .2;
lambda = 1;


% generating features
L = randn(x_dim, rnk);
R = randn(t, rnk);
X0 = L*R';
X0 = X0 / std(X0(:));
X = X0 + 0.1*randn(x_dim,t);
W = randn(y_dim, x_dim)*sqrt(10);
b = randn(y_dim, 1)*sqrt(10);
% masks
Ox = rand(x_dim,t) <= per;
Oy = rand(y_dim,t) <= per;
% labels
Y0 = W*X0 + repmat(b,1,t);
Y = rand(y_dim,t);
Y = Y > 1 ./ (1+(1+exp(-Y0))./(1+exp(Y0)));
Y = sign(Y-1/2);

%%%%%%%%% with bias %%%%%%%%%%
opts = [];
opts.Ox = Ox; opts.Oy = Oy; 
opts.lambda = lambda; opts.L2_wgt = 1; opts.verbose = 1;
[Xhat, Yhat] = transductive_matrix_completion(X, Y, opts);

laberr = nnz((1-Oy).*(Y-sign(Yhat))) / ((y_dim*x_dim)-sum(Oy(:))) * 100;

fprintf(1, 'label error with bias: %f %%\n', laberr);
feaerr = norm((1-Ox).*(Xhat-X),'fro')^2 / norm((1-Ox).*X, 'fro')^2;
fprintf(1, 'feature error with bias: %f \n', feaerr);


%%%%%%%%% without bias %%%%%%%%%%
X = [X; ones(1,t)]; Ox = [Ox; ones(1,t)];% adding 1s
opts.Ox = Ox; opts.method = 0;
[Xhat, Yhat] = transductive_matrix_completion(X, Y, opts);

laberr = nnz((1-Oy).*(Y-sign(Yhat))) / ((x_dim*y_dim)-sum(Oy(:))) * 100;

fprintf(1, 'label error without bias: %f %%\n', laberr);
feaerr = norm((1-Ox).*(Xhat-X),'fro')^2 / norm((1-Ox).*X, 'fro')^2;
fprintf(1, 'feature error without bias: %f \n', feaerr);