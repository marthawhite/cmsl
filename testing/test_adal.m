%TEST
% Test the adal algorithm

% Will first generate a random matrix A, sized m*n
% The objective function for test is
%   L1(X) + reg_wgt*L2(X)

% As cutting plane can be slow, so set m <= 100 and n <= 100
% If you don't want to compare the result of ADAL with cutting plane
% then feel free to use a larger value of m and n
initialize_testing_paths
global sigma_smooth_L11;
global proj_dir;

m = 100;  
n = 100;
sigma_smooth_L11 = 1e-2;
diffTol = 1e-3;

A = randn(m, n);

opts = [];
opts.reg_wgt = 0.1;
opts.funTol = 1e-5;
opts.maxiter = 2000;
opts.verbose = VerboseConst.DETAILED_SOLVER;

% Choose one out of the two losses for L1
% Smoothed L11 norm of X - A
L1 = @(X)smooth_L11(A, X);   

% Square Euclidean distance 0.5/n*||X-A||_F^2;
% L1 = @(X)Euclidean_loss(A, X);    

L1_quad_solver = @(A, X0, opts)sol_quad_general(L1, A, X0, opts);

L2 = @trace_norm;
L2_quad_solver = @sol_quad_trace;

%X0 = zeros(m, n);
X0 = A;
% Solve it by ADAL
X = adal_solver(L1_quad_solver, L2_quad_solver, X0, opts);
X_adal = reshape(X, size(A));

obj_adal = L1(X) + opts.reg_wgt * L2(X);

% Solve it with lbfgs for comparison
init_ls = A(:); %zeros(numel(A), 1);
[X, obj_lbfgs, iter, msg] = lbfgs(@(X)(pobj_adal(X,L1,L2,A,opts.reg_wgt)),...
                             init_ls, opts);
X_lbfgs = reshape(X, size(A));
msg2str(msg, 'lbfgs')

% Compare solutions
maxDiff = max(max(abs(X_adal-X_lbfgs)));
objDiff = abs(obj_lbfgs - obj_adal);
if objDiff > diffTol
  cprintf('red', 'ADAL failed, f_diff = %f, max_arg_diff = %f\n', objDiff,maxDiff);
else
  cprintf('blue', 'ADAL succeeded, with f_diff = %f, max_arg_diff = %f\n', objDiff, maxDiff);
end

