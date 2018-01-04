% Return the parameters to be used by various solvers
function [pm_bfgs_dtnorm, pm_pbm_dtnorm, pm_bfgs_oracle, pm_pbm_oracle,pm_bfgs_tcup] ...
    = get_opt_param ()

% parameters for line search in dual triple norm

pm_bfgs_dtnorm = [];
pm_bfgs_dtnorm.maxiter = 100;     % max number of iterations
pm_bfgs_dtnorm.curvTol = 1e-5;      % tolerance of constraint satisfaction
pm_bfgs_dtnorm.funTol = 1e-5;   % final objective function accuracy pmeter
pm_bfgs_dtnorm.m = 20;

pm_pbm_dtnorm = [];
pm_pbm_dtnorm.maxIter = 1000;     % max number of iterations
pm_pbm_dtnorm.maxLsIter = 20;   % max number of line search steps in each iteration
pm_pbm_dtnorm.maxBdl = 30;      % max number of bundles to keep
pm_pbm_dtnorm.maxFnCall = 10^6;  % max number of calling the function
pm_pbm_dtnorm.tolCon = 1e-4;      % tolerance of constraint satisfaction
pm_pbm_dtnorm.tolFun = 1e-4;   % final objective function accuracy pmeter
pm_pbm_dtnorm.verbose = 0;      % print intermediate steps (doesn't seem to work)


% parameters for line search in strong oracle

pm_bfgs_oracle = [];
pm_bfgs_oracle.maxiter = 100;    % max number of iterations
pm_bfgs_oracle.funTol = 1e-5;    % final objective function accuracy pmeter
pm_bfgs_oracle.m = 100;

pm_pbm_oracle = [];
pm_pbm_oracle.maxIter = 100;     % max number of iterations
pm_pbm_oracle.maxLsIter = 20;    % max number of line search steps in each iteration
pm_pbm_oracle.maxBdl = 100;      % max number of bundles to keep
pm_pbm_oracle.maxFnCall = 100;   % max number of calling the function
pm_pbm_oracle.tolCon = 1e-5;     % tolerance of constraint satisfaction
pm_pbm_oracle.tolFun = 1e-5;     % final objective function accuracy pmeter
pm_pbm_oracle.verbose = 0;       % print intermediate steps (doesn't seem to work)


% parameters for totally corrective update (QP)
pm_bfgs_tcup = [];
pm_bfgs_tcup.maxiter = 100;             % max number of iterations
pm_bfgs_tcup.backtrack_maxiter = 100;   % max number of backtrack iterations
pm_bfgs_tcup.funTol = 1e-5;             % final objective function accuracy pmeter
pm_bfgs_tcup.m = 20;

end
