cd ../
initialize_matlab
cd toy_example
warning off;

rk = 1;
n = 2;
c = 1;
t = 100;

sigma1 = 2; sigma2 = 1;
X = [sigma1*randn(t, 1) sigma2*randn(t,1)];
Y = zeros(t, 1);
sigma = 0.0;
for i = 1:t
  Y(i) = 1 - X(i,2) + sigma*randn(1,1);
end

% Acca and Bcca as directions of correlation given by CCA
X1 = (X-repmat(mean(X),t,1)) ; Y1 = (Y-repmat(mean(Y),t,1));
Z1 = pinv(sqrt(X1'*X1)) * X1';
Z2 = pinv(sqrt(Y1'*Y1)) * Y1';

% CCA
[Acca, Bcca, r, Ucca, Vcca] = canoncorr(X1, Y1);
sum(Ucca-Vcca)
Acca
Bcca

% PCA with rank 1; v is the eigenvector corresponding to the 
% maximum eigenvalue (i.e 1st principal component)
%[V, D] = eig([X Y]'*[X Y]);
%v = V(:,1);
%v
[U,S,V] = svd([X1 Y1],0);
s = diag(S);
umax = U(:,1);
smax = s(1);
vmax = V(:,1);
w = umax;
%phi = smax*vmax';
%[phi, v] = pca([X Y]);
%v = phi/sum(phi);
v = smax*vmax';
v2 = s(2)*V(:,2)';
v
v2

phi = Ucca'; %[0:0.01:1]';

plot_compare_cca(X1, Y1, phi, [Acca;Bcca], v, 'PCA versus CCA');

%figure
%plot3(X1(:,1),X1(:,2),Y1,'b.');
%hold on;
%plot_direction([Acca; Bcca], phi);
%
%plot3(zeros(length(phi),1), phi/Acca(2), phi/Bcca, 'r-');
%hold off;

%%hold on;
%%x1 = phi/v(1); x2 = phi/v(2); x3 = phi/v(3);
%%plot3(x1, x2, x3, 'g+');
%plot_direction(v, phi);
%%plot_direction(v2, phi);
%
%%x1 = phi/v2(1); x2 = phi/v2(2); x3 = phi/v2(3);
%%plot3(x1, x2, x3, 'y+');
%
%hold off;

% Need sigma for smooth_L11 loss
global sigma_smooth_L11;
sigma_smooth_L11 = 1e-3;
reg_wgt = 1;
latent_dim = 1;
reconLoss = @euclidean_loss;

AlgNames = {'Single-View', 'Multi-View', 'Alternator'};
params_all = getParameters(AlgNames, 1);
algs = getAlgs(AlgNames);

opts_recovery = [];
opts_recovery.reg_wgt = reg_wgt;
opts_recovery.latent_dim = latent_dim;

% Run single-view approach
params = params_all{1}{1}; % Only one set of params
params.recover = true;
params.L1 = reconLoss; params.L2 = reconLoss;
params.reg_wgt = reg_wgt;
params.num_basis = 3;
[X_recons Y_recons pobj runtime B_single W_single Phi_single] = algs{1}(Z1, Z2, params);
[X_tr_recons, Y_tr_recons, ~, ~, B_single, W_single, Ph_singlei, pobj] = ...
    adjust_recovery(Z1, Z2,[],[], B_single,W_single,Phi_single,opts_recovery);
B_single
W_single
pobj
plot_compare_cca(X1, Y1, phi, [Acca;Bcca], [B_single;W_single], 'Relaxed Single-View versus CCA');
%plot_direction([B;W], phi);

% Run multi-view approach
params = params_all{2}{1}; % Only one set params
params.recover = true;
params.L1 = reconLoss; params.L2 = reconLoss;
params.reg_wgt = reg_wgt;
params.num_basis = 3;
params.maxiter = 2000;
[X_recons Y_recons pobj runtime B_multi W_multi Phi_multi] = algs{2}(Z1, Z2, params);
[X_tr_recons, Y_tr_recons, ~, ~, B_multi, W_multi, Phi_multi, pobj] = ...
    adjust_recovery(Z1,Z2,[],[], B_multi,W_multi,Phi_multi,opts_recovery);
B_multi
W_multi
pobj
plot_compare_cca(X1, Y1, phi, [Acca;Bcca], [B_multi;W_multi], 'Relaxed Multi-View versus CCA');
%plot_direction([B;W], phi);

% Run alternating approach
params = params_all{3}{1}; % Only one set params
params.recover = true;
params.L1 = reconLoss; params.L2 = reconLoss;
params.reg_wgt = reg_wgt;
params.num_basis = 3;
[X_recons Y_recons pobj runtime B_alt W_alt Phi_alt] = algs{3}(Z1,Z2, params);
[X_tr_recons, Y_tr_recons, ~, ~, B_alt, W_alt, Phi_alt, pobj] = ...
    adjust_recovery(Z1,Z2,[],[], B_alt,W_alt,Phi_alt,opts_recovery);
B_alt
W_alt
pobj
plot_compare_cca(X1, Y1, phi, [Acca;Bcca], [B_alt;W_alt], 'Relaxed Alternator versus CCA');
%plot_direction([B;W], phi);

hold off;