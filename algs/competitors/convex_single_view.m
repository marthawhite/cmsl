function [X_recons, Y_recons, pobj, runtime, B, W, Phi] = ...
                     convex_single_view(X, Y, opts)
% Using global convex approach to solve the sparse coding formulation
%
%    minimize_Phi,U   L1([X; Y],U*Phi) + reg_wgt*||Phi||_2,1
%    s.t.   ||columns of U|| <= gamma
%
%   Approach: minimize equation 4 from 
%      Convex Sparse Coding, Subspace Learning, and Semi-Supervised Extensions
%   min_Zhat L(Zhat, Z) + reg_wgt*sqrt(beta^2 + gamma^2)^{-1/2} ||Zhat||_tr
%      where U = A, Phi = Sigma V' where Zhat = A Sigma V' is SVD of Zhat
%
%   
%   Example: L1 is euclidean loss
%   X has each column as an example
%
%   Requires that all losses also return a gradient
%
%   See DEFAULTS for options.
%
% Author: Xinhua Zhang, University of Alberta, 2011

if nargin < 2
    error('convex_single_view requires at least X and Y');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS STARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULTS.gamma = 1;
DEFAULTS.L1 = @smooth_L11;   % reconstruction loss on X
DEFAULTS.recover = true;
DEFAULTS.reg_wgt = 1e-2;
DEFAULTS.ssl_solver = 'boost';  % 'adal', 'lbfgs', 'boost'
DEFAULTS.verbose = 0;
DEFAULTS.timeout = 1000;

DEFAULTS.solver_pbm_params = struct(...
    'maxIter', 50, ...   	   % max number of iterations
    'maxLsIter', 50, ...	   % max number of line search steps in each iteration
    'maxBdl', 50,	...        % max number of bundles to keep
    'maxFnCall', 10^6, ...     % max number of calling the function
    'tolCon', 1e-5, ...	 % tolerance of constraint satisfaction
    'funTol', 1e-5 ...	   % final objective function accuracy parameter
    );

DEFAULTS.solver_lbfgs_params = struct(...
    'maxiter', 2000, ...		 % max number of iterations
    'funTol', 1e-4, ...		 % final objective function accuracy parameter
    'm', 50 ...
    );

DEFAULTS.solver_adal_params = struct(...
    'maxiter', 2000, ...		 % max number of iterations
    'funTol', 1e-5 ...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
  opts = DEFAULTS;
else  
  opts = getOptions(opts, DEFAULTS);
end
opts = addToAllOptions(opts, struct('reg_wgt', opts.reg_wgt,...
          'timeout', opts.timeout, 'verbose', opts.verbose));

% Rescale reg_wgt for the convex sparse coding optimization
if opts.reg_wgt <= 1e-2
  opts.ssl_solver = 'adal';
end
gamma2 = opts.gamma^2;
opts.reg_wgt = opts.reg_wgt/sqrt(1+gamma2);

t = size(X,2);        % # total examples
n = size(X,1);        % # features 1
r = size(Y,1);       % # features 2
XY = [X; Y];

if opts.verbose >= VerboseConst.BASIC_ALG 
  cprintf(VerboseConst.BASIC_ALG_COLOUR,...
          ['\n\nconvex_single_view -> Starting with solver = %s...\n'], opts.ssl_solver); 
end

start_time = cputime;

sizeZhat = [n+r, t];
%Zhat = zeros(sizeZhat);
Zhat = XY;

if strcmp(opts.ssl_solver, 'fmin_LBFGS') || strcmp(opts.ssl_solver, 'lbfgs')
  [Zhat_opt, pobj, iter, flag] = lbfgs( @single_view_loss, Zhat(:), opts.solver_lbfgs_params); 
  Zhat = reshape(Zhat_opt, sizeZhat);
elseif strcmp(opts.ssl_solver, 'adal')

    L1 = @(Zhat)opts.L1(XY, Zhat);

    L1_quad_solver = @(A, X0, opts)sol_quad_general(L1, A, X0, opts);
    L2_quad_solver = @sol_quad_trace;

    % Solve it by ADAL
    Zhat = adal_solver(L1_quad_solver, L2_quad_solver, Zhat, opts.solver_adal_params);

    pobj = L1(Zhat) + opts.reg_wgt * trace_norm(Zhat);
elseif strcmp(opts.ssl_solver, 'pbm')
  ub = inf((n+r)*t,1);
  [Zhat_opt, pobj, iter, feval, msg] = pbm(@single_view_loss, Zhat(:), -ub, ub, opts.solver_pbm_params);
  Zhat = reshape(Zhat_opt, sizeZhat);  
else % Solve it with boosting, with inner solver as lbfgs
  [Zhat, pobj, iter, msg] = trace_reg_pboost(@(Zhat)(opts.L1(XY, Zhat)), sizeZhat, opts);
end

runtime = cputime - start_time;
X_recons = Zhat(1:n,:);
Y_recons = Zhat(n+1:end, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = cputime;
% recover B, W and Phi from Z
if nargout > 4
    if opts.recover
        [Usvd,Sigma,V] = svd(Zhat, 'econ'); 
        U = Usvd*sqrt(1+gamma2);
        Phi = (Sigma * V')/sqrt(1+gamma2);
        B = U(1:n, :);
        W = U((n+1):end, :);
    else
        B = []; W = []; Phi = [];
    end    
end

runtime = [runtime cputime-start_time];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [f,G] = single_view_loss(Zhatvec) 
    Zhat = reshape(Zhatvec, sizeZhat);
    if nargout < 2
      f1 = opts.L1(XY, Zhat);
      f3 = trace_norm(Zhat);
    else  
      [f1, G1] = opts.L1(XY, Zhat);
      [f3,G3] = trace_norm(Zhat);
      G = G1(:)  + opts.reg_wgt*G3(:);    
    end
    f = f1 + opts.reg_wgt*f3;    
end


end  
