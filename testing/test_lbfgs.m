%TEST
% Test the lbfgs algorithms

% Will first generate a random matrix A, sized m*n
% The objective function for test is
%   L1(X) + reg_wgt*L2(X)

initialize_testing_paths
global proj_dir;

m = 100;  
n = 100;
diffTol = 1e-3;

X = randn(m, n);

opts = [];
opts.reg_wgt = 0.01;
opts.funTol = 1e-5;
opts.maxiter = 2000;
opts.verbose = VerboseConst.DETAILED_SOLVER;

L1 = @(Xhat)euclidean_loss(X, Xhat);    

L2 = @trace_norm;

Xinit = X(:);
loss = @(Xhat)(pobj_adal(Xhat,L1,L2,X,opts.reg_wgt));

% Solve it with lbfgs for comparison
rmpath([proj_dir '/solver/lbfgsb']);
addpath([proj_dir '/solver/uofa_optimizers/bfgs']);
starttime = cputime;
[X_lbfgs, obj_lbfgs, iter, msg] = lbfgs(loss, Xinit, opts);
fprintf(1, 'Runtime of UofA LBFGS = %g, iters = %d\n', cputime - starttime, iter);
X_lbfgs = reshape(X_lbfgs, size(X));
msg2str(msg, 'lbfgs')

% Solve it with mex implementation of lbfgs
rmpath([proj_dir '/solver/uofa_optimizers/bfgs']);
addpath([proj_dir '/solver/lbfgsb']);
%[X_lbfgs2,obj_lbfgs2,iter,feval,msg] = lbfgsb(Xinit, -Inf(size(Xinit)), Inf(size(Xinit)),
%loss,[],[],struct('maxIters', 50));
starttime = cputime;
[X_lbfgs2, obj_lbfgs2, iter, msg] = lbfgs(loss, Xinit, opts);
fprintf(1, 'Runtime of MEX LBFGS = %g, iters = %d\n', cputime - starttime, iter);
X_lbfgs2 = reshape(X_lbfgs2, size(X));
msg2str(msg, 'lbfgs')

% Compare solutions for two lbfgs implementations
maxDiff = max(max(abs(X_lbfgs-X_lbfgs2)));
objDiff = abs(obj_lbfgs - obj_lbfgs2);
if objDiff > diffTol
  cprintf('red', 'LBFGS results different, f_diff = %f, max_arg_diff = %f\n', objDiff,maxDiff);
else
  cprintf('blue', 'LBFGS results suitably similar, with f_diff = %f, max_arg_diff = %f\n', objDiff, maxDiff);
end

