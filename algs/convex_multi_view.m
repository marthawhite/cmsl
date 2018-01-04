function [X_learned, Y_learned, pobj, runtime, B, W, Phi] = ...
    convex_multi_view(X, Yl, opts)
% Using global convex approach to solve the sparse coding formulation
%
%	 minimize_Phi,B,W	L1(X,B*Phi) + L2_wgt*L2(Y_L,W*Phi) + reg_wgt*||Phi||_2,1
%	 s.t.	||columns of B|| <= 1, ||columns of W || <= gamma
%
%	Approach: line search in the primal
%	max_eta min_Q L(E_{eta} Q) + ||Q||_tr
%
%	Example: L1 is l2 loss, L2 is hinge loss
%	X has each column as an example
%
%	Requires that all losses also return a gradient
%
%	See DEFAULTS below for all the optional parameters.
%
% Author: Xinhua Zhang, University of Alberta, 2011



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run

if nargin < 2
  error('convex_multi_view requires at least X and Yl');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS STARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFAULTS.reg_wgt = 1e-3;
DEFAULTS.L2_wgt = 1;
DEFAULTS.gamma = 1;
DEFAULTS.L1 = @euclidean_loss;	 % reconstruction loss on X
DEFAULTS.L2 = @logistic_loss;	 % classification loss on Y
DEFAULTS.recover = true;
DEFAULTS.eta_lower_bound = 0.1;

DEFAULTS.verbose = 0;	% 0: nothing
                      % 1: message of the outer solver (line search)
                      % 2: message of the inner solver (eg. neseterov)
DEFAULTS.timeout = -1;  % Stop optimization after specified cputime in seconds
                        % Only works if outer_solver is lbfgs

DEFAULTS.outer_solver = 'lbfgs';	  % 'pbm', 'lbfgs'
DEFAULTS.inner_solver = 'boost';  % 'boost', 'pbm', 'lbfgs', 'adal', or 'nesterov'

DEFAULTS.outer_solver_pbm_params = struct(...
    'maxIter', 20, ...   	   % max number of iterations
    'maxLsIter', 20, ...	   % max number of line search steps in each iteration
    'maxBdl', 30,	...        % max number of bundles to keep
    'maxFnCall', 20, ...     % max number of calling the function
    'tolCon', 1e-2, ...	 % tolerance of constraint satisfaction
    'funTol', 1e-3 ...	   % final objective function accuracy parameter
    );

DEFAULTS.outer_solver_lbfgs_params = struct(...
    'maxiter', 100, ...	    % max number of iterations
    'funTol', 1e-3, ...		% final objective function accuracy parameter
    'backtrack_init_stepsize', 1e-1, ... % eta more sensitive to gradient, make step size small
    'backtrack_acceptfrac', 1, ... % eta more sensitive to gradient, make step size small
    'backtrack_maxiter', 5, ... % eta more sensitive to gradient, make step size small
    'm', 50 ...
    );

DEFAULTS.inner_solver_nesterov_params = struct(...
    'L_amp', 1.5, ...			     % how much L is multiplied in Nesterov's trial
    'tol_nesterov', 1e-4, ...
    'L0_nest', 1e-3, ...	     % intial L for nesterov
    'maxiter', 1000 ...
    );

DEFAULTS.inner_solver_pbm_params = DEFAULTS.outer_solver_pbm_params;
DEFAULTS.inner_solver_pbm_params.tolCon = 1e-3;
DEFAULTS.inner_solver_pbm_params.funTol = 1e-3;
DEFAULTS.inner_solver_pbm_params.maxIter = 1000;

% These inner parameters are also used for boosting
% which implicitly using lbfgs for the inner loop
DEFAULTS.inner_solver_lbfgs_params = struct(...
    'maxiter', 1000, ...		 % max number of iterations
    'funTol', 1e-5, ...		 % final objective function accuracy parameter
    'm', 50 ...
    );

DEFAULTS.inner_solver_adal_params = struct(...
    'maxiter', 1000, ...		 % max number of iterations
    'funTol', 1e-5 ...
    );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
  opts = DEFAULTS;
else
  opts = getOptions(opts, DEFAULTS);
end
opts = addToAllOptions(opts, struct('timeout', opts.timeout,...
                                    'verbose', opts.verbose));

t = size(X,2);		  % # total examples
tl = size(Yl,2);	  % # labeled examples
tu = t - tl;		  % # unlabeled examples
n = size(X,1);		  % # features
r = size(Yl,1);		  % # classes
g2 = opts.gamma^2;

% Adjust maximum number inner iterations and tolerance
% depending on current gradient for eta, so that we can
% concentrate computation when closing in on correct eta
current_eta_g = 10;
prev_eta = 0;

% Choose the solver to solve the inner problem in Q (eta fixed)
% Ie, evaluate f(eta).
% TODO: Adal currently broken, figure out why
if opts.reg_wgt <= 1e-2
  opts.inner_solver = 'adal';
end
switch opts.inner_solver
  case 'pbm'
    pobj_eta = @pobj_eta_pbm;		  % Solve the inner problem (in Q) by pbm
  case 'nesterov'
    pobj_eta = @pobj_eta_nesterov;	% Solve the inner problem (in Q) by Nesterov
  case {'lbfgs', 'fmin_LBFGS'}
    pobj_eta = @pobj_eta_lbfgs;	  % Solve the inner problem (in Q) by lbfgs       
  case 'adal'
    pobj_eta = @pobj_eta_adal;	  % Solve the inner problem (in Q) by adal   
    opts.inner_solver_adal_params.reg_wgt = opts.reg_wgt;
  case 'boost'
    pobj_eta = @pobj_eta_boost;	  % Solve the inner problem (in Q) by boosting
  otherwise
    error('convex_multi_view -> Invalid inner optimizer %s', opts.inner_solver);
end

switch opts.outer_solver
  case 'pbm'    
    outer_solver = @(init_eta, lb, ub)(pbm(pobj_eta, init_eta, lb, ub, opts.outer_solver_pbm_params));
    if opts.outer_solver_pbm_params.timeout > 0
      warning('convex_multi_view -> Timeout not implement for outer_solver = pbm, use lbfgs.');
    end  
  case {'lbfgs', 'fmin_LBFGS'}
    %outer_solver = @(init_eta, lb, ub)(lbfgs(pobj_eta, init_eta, opts.outer_solver_lbfgs_params));
    outer_solver = @(init_eta, lb, ub)(lbfgs(pobj_eta, init_eta, ...
                                       getOptions(opts.outer_solver_lbfgs_params, ...
                                       struct('A', [1; -1], 'b', [lb; -ub])))); 
  otherwise
    error('convex_multi_view -> Invalid outer optimizer %s', opts.outer_solver);
end

if opts.verbose >= VerboseConst.BASIC_ALG, 
  cprintf(VerboseConst.BASIC_ALG_COLOUR,...
         ['\n\nconvex_multi_view -> Starting with outer solver = %s' ...
         ' and inner solver = %s.\n'], opts.outer_solver, opts.inner_solver); 
end

%%%%%%%%%%%%%%%%%%%%%%%%% START OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables
% Ensure initial guess is bounded within ub and lb
init_eta = 1/(1+opts.gamma);
lb = min(opts.eta_lower_bound, init_eta - 0.05);
ub = max(1 - lb, init_eta);
lb = min(lb, 1-ub);
Q = initQ(X, Yl, init_eta);

% START the outer optimization
if opts.verbose >= VerboseConst.BASIC_ALG, fprintf('Outer line search in eta\n'); end
start_time = cputime;
[eta_opt, pobj, iter, msg] = outer_solver(init_eta, lb, ub);
runtime = cputime - start_time;

% Check if solution is a lower or upper bound, in case bound is too tight
if opts.verbose >= VerboseConst.BASIC_ALG
  if eta_opt <= lb
    cprintf(VerboseConst.WARNING_COLOUR,...
    'convex_multi_view -> optimal eta %g is at lower boundary %g\n', eta_opt, lb);
  elseif eta_opt >= ub  
    cprintf(VerboseConst.WARNING_COLOUR,...
            'convex_multi_view -> optimal eta %g is at upper boundary %g\n', eta_opt, ub);
  end
  fprintf('Outer line-search: iter = %d, obj = %f, eta = %g, runtime = %f sec, msg = %s\n', ...
          iter, pobj, eta_opt, runtime, msg2str(msg,opts.outer_solver));
end

% DONE recover solution Z, X_recons, Yl, Yu_pred from line search result
Q = reshape(Q, [n+r, t]);
[v1, v2] = eta_map(eta_opt);
Z_opt = [v1*Q(1:n,:); v2*Q(n+1:end,:)];
X_recons = Z_opt(1:n,:);
Yl_recons = Z_opt(n+1:end, 1:tl);
Yu_pred = Z_opt(n+1:end, tl+1:end);

pobj = -pobj;	  % negate back since max_{eta}

% Returned reconstructed matrices; in the case of semi-supervisde learning
% Yl is transductively filled.
X_learned = X_recons;
Y_learned = [Yl_recons];
if ~isempty(Yu_pred)
  Y_learned = [Y_learned; Yu_pred];
end  

%%%%%%%%%%%%%%%%%%%%%%  END OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN RECOVERY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = cputime;
if nargout > 4
  if opts.recover
    fprintf('recovering B, W, Phi\n');
    % recover B, Phi and W from X_Pred and Y_pred
    [U Phi] = cone_recover (Z_opt, [n, r, t], opts);
    B = U(1:n, :); W = U(n+1:end, :);
  else
    B = [];	  W = [];	Phi = [];
  end
end
runtime = [runtime cputime-start_time];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END RECOVERY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute a subgradient of the outer function in eta
% Use pbm (guaranteed to find an optimal solution Q)
function [f,g] = pobj_eta_pbm(eta)
  
  % Given eta, optimize over Q by using Nesterov
  ub = Inf(numel(Q), 1);		lb = -ub;
  
  % must call pbm1 (copy the mex file of pbm to pbm1)
  % matlab has problem in using a mex function within itself (remember
  % the outer loop may also use pbm)
  [Q, f, iter, feval, msg] = pbm1(@(W)pobj_Q_tr(W, eta), Q(:), ...
                                  lb, ub, opts.inner_solver_pbm_params);
  
  % After finding the optimal Q, compute the gradient in eta
  f = -f;	  % negate it since outer loop is maximization
  g = -grad_eta(eta, reshape(Q, [n+r,t]));
  
  if opts.verbose >= VerboseConst.BASIC_ALG
    cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
            'out (pbm): eta = %f, f = %f, g= %f, iter = %d, feval = %d\n', ...
            eta, f, g, iter, feval);
    cprintf(VerboseConst.BASIC_ALG_COLOUR, 'msg = %s\n', msg2str(msg, 'pbm'));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a subgradient of the outer function in eta
% Use lbfgs (no guarantee of convergence since the obj is not smooth)
function [f,g] = pobj_eta_lbfgs(eta)
  
  % Given eta, optimize over Q by using lbfgs
  [Q, f, iters, flag] = lbfgs(@(W)pobj_Q_tr(W, eta), Q(:),...
                             opts.inner_solver_lbfgs_params);
  
  % After finding the optimal Q, compute the gradient in eta
  f = -f;	  % negate it since outer loop is maximization
  if nargout > 1
    g = -grad_eta(eta, reshape(Q, [n+r,t]));
  end
   
  if opts.verbose >= VerboseConst.BASIC_ALG
    if nargout > 1
      cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
              'out (lbfgs): eta = %f, f = %f, g= %f, msg = %s\n', ...
              eta, f, g, msg2str(flag, 'lbfgs'));
    else
      cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
              'out (lbfgs): eta = %f, f = %f, msg = %s\n', eta, f, msg2str(flag, 'lbfgs'));
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a subgradient of the outer function in eta 
% Use adal (guaranteed convergence)
function [f,g] = pobj_eta_adal(eta)
  
  % Reduce the amount of time given to adal, given time already used
  if opts.timeout > 0
    opts.inner_solver_adal_params.timeout = max(10, opts.timeout - (cputime - start_time));
  end
  if abs(current_eta_g) > 0.5
    opts.inner_solver_adal_params.maxiter = 100;
    opts.inner_solver_adal_params.funTol = 1e-3;
  else
    opts.inner_solver_lbfgs_params.maxiter = 1000;
    opts.inner_solver_lbfgs_params.funTol = 1e-4;    
  end

  L1 = @(W)pobj_Q(W, eta);

  L1_quad_solver = @(A, X0, opts)sol_quad_general(L1, A, X0, opts);
  L2_quad_solver = @sol_quad_trace;

  % Solve it by ADAL
  %Q = initQ(X, Yl, eta);   
  adal_starttime = cputime;
  [Q Lambda] = adal_solver(L1_quad_solver, L2_quad_solver,...
                           Q, opts.inner_solver_adal_params);
  opts.inner_solver_adal_params.Lambda = Lambda;

  f = -L1(Q) - opts.reg_wgt * trace_norm(Q);
  if nargout > 1
    g = -grad_eta(eta, reshape(Q, [n+r,t]));
    current_eta_g = g;
  else
    g = 0;
  end
  
  if opts.verbose >= VerboseConst.BASIC_ALG
    if nargout > 1
	    cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
            'out (adal): eta = %f, f = %f, g= %f, runtime = %g\n',...
            eta, f, g, cputime-adal_starttime);
    else
	    cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
            'out (adal): eta = %f, f = %f, runtime = %g\n',...
            eta, f, cputime-adal_starttime);
    end
  end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a subgradient of the outer function in eta
% Use boosting
function [f,g] = pobj_eta_boost(eta)
  % Given eta, optimize over Q by using lbfgs
  
  %if abs(current_eta_g) > 0.1
  %  opts.inner_solver_lbfgs_params.maxiter = 100;
  %  opts.inner_solver_lbfgs_params.funTol = 1e-1;
  %else
  %  opts.inner_solver_lbfgs_params.maxiter = 1000;
  %  opts.inner_solver_lbfgs_params.funTol = 1e-4;    
  %end
  %opts.use_local = false;
  sizeX = [n+r; t];
  [Q, f, iter, msg] = trace_reg_pboost(@(Q)pobj_Q(Q, eta), sizeX, opts);
  
  % After finding the optimal Q, compute the gradient in eta
  f = -f;	  % negate it since outer loop is maximization
  if nargout > 1
    g = -grad_eta(eta, reshape(Q, [n+r,t]));
    %current_eta_g = g;
  else
    g = 0;
  end
  
  if opts.verbose >= VerboseConst.BASIC_ALG
    if nargout > 1
      cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
              'out (boost): eta = %f, f = %f, g= %f, msg = %s\n', eta, f, g, msg);
    else
      cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
              'out (boost): eta = %f, f = %f, msg = %s\n', eta, f, msg);
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute a subgradient of the outer function in eta
% Use Nesterov (guaranteed to find an optimal solution Q)
function [f,g] = pobj_eta_nesterov(eta)
  
  % Given eta, optimize over Q by using Nesterov
  Z=Q;	L=opts.inner_solver_nesterov_params.L0_nest/1.5;
  a=1;   i=0;	   f_pre = inf;
  best_obj = inf;	  Q_opt = Q;
  max_iters_worse = 100;
  iters_worse = 0;
  max_inner_iters = 1000;
  
  while i < opts.inner_solver_nesterov_params.maxiter
    [fZ, G] = pobj_Q(Z, eta);
    i = i+1;
    count = 0;
    while count < max_inner_iters
      % Compute p_L(Z)
      [U S V] = svd(Z - (1/L)*G, 'econ');
      diagS = diag(S) - opts.reg_wgt/L;
      rk = sum(diagS > 0);
      diagS = max(0, diagS);
      pLZ = (U.*repmat(diagS', n+r, 1))*V';
      
      % Compute F(pLZ) and Q(pLZ)
      f1 = pobj_Q(pLZ, eta);
      T = pLZ-Z;
      f2 = fZ + sum(sum(G.*T)) + L/2*norm(T,'fro')^2;
      
      if f2 > f1, break; else L = L*opts.inner_solver_nesterov_params.L_amp; end
      count = count +1;
    end

    b=a;	a = 0.5*(1+sqrt(1+4*a^2));
    Z = pLZ + ((b-1)/a)*(pLZ - Q);
    Q = pLZ;
    f = f1 + opts.reg_wgt*sum(diagS);
    if f < best_obj
      best_obj = f;		 Q_opt = Q;
    else
      iters_worse = iters_worse+1;
    end
    if opts.verbose >= VerboseConst.DETAILED_ALG, fprintf('%d: %g, L=%g, rk = %d\n', i, f, L, rk); end
    if (abs(f_pre-f) < opts.inner_solver_nesterov_params.tol_nesterov || ...
        iters_worse > max_iters_worse),  break; end
    f_pre = f;
    L = L / 1.1;
  end
  % After finding the optimal Q, compute the gradient in eta
  f = -best_obj;	 % negate it since outer loop is maximization
  if nargout < 2
    g = 0;
  else 
    g = -grad_eta(eta, Q_opt);
  end  
  Q = Q_opt;
  if opts.verbose >= VerboseConst.BASIC_ALG
    cprintf(VerboseConst.BASIC_ALG_COLOUR, 'out (Nesterov): eta = %f, f = %f, g= %f, L = %g, iter = %d\n', ...
            eta, f, g, L, i);
  end
end

function g = grad_eta(eta, Q)
  [v1 v2] = eta_map(eta);
  Q_X = Q(1:n,:);
  Q_Y = Q(n+1:end,1:tl);
  [~, G1] = opts.L1(X, v1*Q_X);
  [~, G2] = opts.L2(Yl, v2*Q_Y);
  
  [v1 v2] = eta_grad(eta);
  g = v1*sum(sum(Q_X.*G1)) + opts.L2_wgt*v2*sum(sum(Q_Y.*G2));
end

% For a fixed rho, compute the gradient of loss in Q without the trace norm
function [f, G] = pobj_Q(Q, eta)
  [v1 v2] = eta_map(eta);
  if nargout > 1
    [f1, G1] = opts.L1(X, v1*Q(1:n,:));
    [f2, G2] = opts.L2(Yl, v2*Q(n+1:end,1:tl));
    G = [v1*G1; [(opts.L2_wgt*v2)*G2, zeros(r,tu)]];    
  else
    f1 = opts.L1(X, v1*Q(1:n,:));
    f2 = opts.L2(Yl, v2*Q(n+1:end,1:tl));
  end
  f = f1 + opts.L2_wgt*f2;
end

% For a fixed eta, compute the gradient of loss in Q including the trace norm
function [f, g] = pobj_Q_tr(Q, eta)
  Q_mat = reshape(Q, [n+r, t]);
  [v1 v2] = eta_map(eta);
  if nargout < 2
    f1 = opts.L1(X, v1*Q_mat(1:n,:));
    f2 = opts.L2(Yl, v2*Q_mat(n+1:end,1:tl));
    f3 = trace_norm(Q_mat);
  else  
    [f1, G1] = opts.L1(X, v1*Q_mat(1:n,:));
    [f2, G2] = opts.L2(Yl, v2*Q_mat(n+1:end,1:tl));
    [f3, G3] = trace_norm(Q_mat);
    g = [v1*G1; [(opts.L2_wgt*v2)*G2, zeros(r,tu)]] + opts.reg_wgt*G3;
    g = g(:);
  end
  f = f1 + opts.L2_wgt*f2 + opts.reg_wgt*f3;
end

function Q = initQ(X, Yl,eta) 
  [v1 v2] = eta_map(eta);
  Q = [X / v1; [Yl / v2, zeros(r,tu)]];  
  %verbose = opts.verbose; opts.verbose = 0;
  %[Q1,Q2] = convex_single_view(X,Yl,opts);
  %opts.verbose = verbose;
  %Q = [Q1;Q2];
  
  %Q = sol_quad_trace([X/v1; Yl/v2], [], opts);
  
  %Q = zeros(n+r, t); 
  
  %Q = rand(n+r, t); 
end

function [v1 v2] = eta_map(eta)
  v1 = sqrt(1+opts.gamma*(1-eta)/eta);
  v2 = sqrt(g2+opts.gamma*eta/(1-eta));
end

function [v1 v2] = eta_grad(eta)
  v1 = -opts.gamma/(2*eta*sqrt(eta^2+opts.gamma*eta*(1-eta)));
  v2 = opts.gamma/(2*(1-eta)*sqrt(g2*(1-eta)^2+opts.gamma*eta*(1-eta)));
end
end
