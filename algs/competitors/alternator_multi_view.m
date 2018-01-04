function [X_learned, Y_learned, obj, runtime, B, W, Phi] = ... 
        alternator_multi_view(X, Yl, opts)
% Using alternating to solve the sparse coding formulation
%   Note: H = B, W = Phi and U = W in our new notation
%
%    minimize_B,Phi,W   L1(X,B*Phi) + L2_wgt*L2(Y_L, W*Phi_L) + reg_wgt*reg_loss(Phi_i:)
%       s.t.   ||columns of B|| <= 1, ||columns of W || <= gamma
%
%   Approach: 1) Fix B & W, minimize Phi. 2) Fix B and Phi, minimize
%   W. 3) Fix W and Phi and minimize B. Continue three steps until
%   convergence. Notice B and W are independent of each other given Phi.
% 
%   Example: L1 and L2 are the l2 loss and regularizer_loss is l1 loss
%   X has each column as a sample
%
%   Requires that all losses also return a gradient, though it is
%   optional to supply the losses
%   Default: L1 = l2 loss, L2 = hinge loss, reg_loss = 2-1 norm
%
%	See DEFAULTS below for all the optional parameters.
%
% Author: Martha White, University of Alberta, 2012


global col_print;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS STARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sparse coding parameters
DEFAULTS.L2_wgt = 1;  % weight on L2, loss between Yl and W*Phi
DEFAULTS.reg_wgt = 0.1;  % weight on regularizer on Phi
DEFAULTS.num_basis = 50;  % size(Phi,1)
DEFAULTS.gamma = 1;  % constraint on columns of W

% Losses
DEFAULTS.L1 = @euclidean_loss;
DEFAULTS.L2 = @logistic_loss;
DEFAULTS.reg_loss = @L21_loss;

% Iteration parameters
DEFAULTS.num_reps = 3;
DEFAULTS.maxiter = 500;
DEFAULTS.funTol = 1e-5;
DEFAULTS.timeout = -1; % Default is no timeout, i.e. = -1; else, cputime seconds e.g. 10000
DEFAULTS.verbose = 0;

% opts.solver specifies what solver to use for all three variables.
% By default, opts.solver is used for B, Phi and W unless the
% user specifies specific solvers for B, Phi and W.
% Current options include 'lbfgs', 'nesterov' and otherwise
% if not one of these, assume that given a function handle of the form
% Phi = opts.solver_Phi(Phi)
DEFAULTS.solver = 'nesterov';
DEFAULTS.inner_solver_params = [];
DEFAULTS.solver_B = []; 
DEFAULTS.solver_W = []; 
DEFAULTS.solver_Phi = []; 
DEFAULTS.init_with_singleview = 0;  % It's cheating to initialize with
                                    % with single-view, but a useful test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
  opts = DEFAULTS;
else  
  opts = getOptions(opts, DEFAULTS);
end
opts.inner_solver_params.verbose = opts.verbose;

AlgName = 'alternator_multi_view';
if opts.L2_wgt == 0, AlgName = 'staged_alternator_multi_view'; end
if opts.verbose >= VerboseConst.BASIC_ALG, 
  cprintf(VerboseConst.BASIC_ALG_COLOUR,['\n' AlgName '-> Starting...\n\n']); 
end

if opts.init_with_singleview, opts.num_reps = 1; end

% Fixed parameters
t = size(X,2);
tl = size(Yl,2);
tu = t-tl;
n = size(X,1);
r = size(Yl,1);

sizePhi = [opts.num_basis t];
sizeB = [n opts.num_basis];
sizeW = [r, opts.num_basis];

obj_Phi = @obj_Phi_alt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine which solver to use for B, W and Phi
  % Need function to reshape results for lbfgs, with takes vectors not matrices  
  function result = reshape_lbfgs(solver, fun, x, dimensions)
    function [f,g] = fun2(xvec)
      xmat = reshape(xvec, dimensions);
      if nargout > 1
        [f,g] = fun(xmat);
        g = g(:);
      else
        f = fun(xmat);
      end
    end
    result = solver(@fun2,x(:));
    result = reshape(result, dimensions);
  end
if strcmp(opts.solver, 'lbfgs') || strcmp(opts.solver, 'fmin_LBFGS')
  
  opts.solver = @(fun, x)(lbfgs(fun, x, opts.inner_solver_params));		
  if opts.L2_wgt == 0, obj_Phi = @obj_Phi_staged; end
  opts.solver_B = @(x)(reshape_lbfgs(opts.solver, @obj_B, x, sizeB));
  opts.solver_Phi = @(x)(reshape_lbfgs(opts.solver, obj_Phi, x, sizePhi));
  opts.solver_W = @(x)(reshape_lbfgs(opts.solver, @obj_W, x, sizeW));   
  
else % Else does generic nesterov
	phiopts = opts;   phiopts.init_L = 1 / t;   phiopts.reg_wgt = opts.reg_wgt;
  phiopts.verbose = 0;
	Wopts = phiopts;  Wopts.radius_ball = opts.gamma;
	Bopts = phiopts;  Bopts.radius_ball = 1;
  if opts.L2_wgt == 0, 
    obj_Phi = @obj_Phi_noreg_staged;
  else
    obj_Phi = @obj_Phi_noreg;    
  end  

  opts.solver_B = @(B)(solve_Nesterov_generic(@obj_B, [], @prox_op_L2ball, B, Bopts));
  opts.solver_Phi = @(x)(solve_Nesterov_generic(obj_Phi, opts.reg_loss,...
                                                @prox_op_L21, x, phiopts));
  opts.solver_W = @(W)(solve_Nesterov_generic(@obj_W, [], @prox_op_L2ball, W, Wopts));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start looping over random restarts

runtime = zeros(2, 1);  % [1]: sum of time,  [2] argmin
fval = zeros(opts.maxiter, 1);
B_rec = [];
Phi_rec = [];
W_rec = [];

obj_rec = inf;
diff_rec = 0;
diff = 0;
iter_mod = floor((opts.num_reps * opts.maxiter) / 100);

for rep = 1 : opts.num_reps

  start_time = cputime;

  % Mean zero, variance 1 weights
  B = initWeights(sizeB);    
  
  % If opts.L2_wgt = 0, then this step skipped as doing staged learning rather than
  % alternating over all three variables    
  W = [];
  if opts.L2_wgt ~= 0
    W = initWeights(sizeW);
  end

  % Initializing basis according to Lee&Ng's work  
  Phi = 2*(rand(sizePhi) - 0.5);
  
  if opts.init_with_singleview
    old_recover = opts.recover; opts.recover = 1;
    [~,~,~,~,B,W,Phi] = convex_single_view(X,Yl,opts);
    opts.recover = old_recover;
    B = [B zeros(sizeB(1), sizeB(2)-size(B,2))];
    W = [W zeros(sizeW(1), sizeW(2)-size(W,2))];
    Phi = [Phi; zeros(sizePhi(1)-size(Phi,1), sizePhi(2))];
  end
  
  % If not staged alternator, add in the loss associated with Yl
  prevOpt = opts.L1(X,B,Phi,0) + opts.reg_wgt*opts.reg_loss(Phi);
  if opts.L2_wgt ~= 0
    prevOpt = prevOpt +  opts.L2_wgt * opts.L2(Yl,W,Phi(:,1:tl),0);
  end

  for iter = 1 : opts.maxiter

    % Start with Phi; if opts.L2_wgt = 0, does not include L2 loss
    Phi = opts.solver_Phi(Phi);
    
    % Next is B, ignore constant_val = opts.L2(Yl,U,Phi(:,1:tl) + opts.reg_wgt*opts.reg_loss(Phi)    
    B = opts.solver_B(B);

    % Next is W, ignore constant_val = opts.L1(X,B,Phi) + opts.reg_wgt*opts.reg_loss(W)
    % If opts.L2_wgt = 0, then this step skipped as doing staged learning rather than
    % alternating over all three variables
    newOpt = opts.L1(X, B, Phi, 0) + opts.reg_wgt * opts.reg_loss(Phi);
    if opts.L2_wgt ~= 0
      W = opts.solver_W(W);
      newOpt = newOpt + opts.L2_wgt * opts.L2(Yl, W, Phi(:, 1:tl), 0);
    end
    
    fval(iter) = newOpt;
    if opts.verbose >= VerboseConst.DETAILED_ALG && mod(iter-1, iter_mod) == 0
      cprintf(VerboseConst.DETAILED_ALG_COLOUR, '%d: obj = %g, rel change = %g\n',...
              iter, newOpt, abs(newOpt-prevOpt)/min(newOpt, prevOpt));
    end
    
    if newOpt > prevOpt && opts.verbose >= VerboseConst.BASIC_ALG
      cprintf(VerboseConst.WARNING_COLOUR, [AlgName '-> Alternating minimization ' ...
                 'increased function value from %g to %g on iter %u\n'], prevOpt, newOpt, iter);
    end   
    if abs(newOpt-prevOpt)/min(newOpt, prevOpt) < opts.funTol || ...
            (abs(prevOpt) < 1e-7 && abs(newOpt-prevOpt) < opts.funTol)
      diff = abs(newOpt-prevOpt);
      break;
    end  
    
    prevOpt = newOpt;
    
    if opts.timeout > 0 && opts.timeout < (cputime - start_time), break; end
  end
  
  rtime = cputime - start_time;
  runtime(1) = runtime(1) + rtime;
  
  if newOpt < obj_rec
    runtime(2) = rtime;
    obj_rec = newOpt;
    diff_rec = diff;
    B_rec = B;  W_rec = W;  Phi_rec = Phi;
  end

  if opts.verbose >= VerboseConst.BASIC_ALG
    cprintf(VerboseConst.BASIC_ALG_COLOUR, '>>>>>  rep = %d, obj = %g, tol = %g, iter = %d, time = %f\n',...
              rep, newOpt, diff, iter, rtime);
  end
  
  if opts.timeout > 0 && runtime(2) > opts.timeout
    cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
            [AlgName '-> Timeout after %g seconds\n'], opts.timeout);
    break;
  end  
end     % random initialization

% If opts.L2_wgt == 0, then doing staged learning, so compute W now
if opts.L2_wgt == 0
  W = initWeights(sizeW);  
  W_rec = opts.solver_W(W);
end

obj = obj_rec;
B = B_rec;
W = W_rec;
Phi = Phi_rec;
if opts.verbose >= VerboseConst.BASIC_ALG
  cprintf(VerboseConst.BASIC_ALG_COLOUR, ...
          [AlgName '-> Best optimization terminated with tolerance = %g ' ...
              '(lower than min tol %g)\n'], diff_rec, opts.funTol);
end 

% Now predict on unlabeled data
Yu_pred = W * Phi(:,(tl+1):end);

X_learned = B * Phi;
Y_learned = [W * Phi(:, 1:tl)];
if ~isempty(Yu_pred)
  Y_learned = [Y_learned; Yu_pred];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Objective functions

% Assumes opts.L2_wgt ~= 0
function [f, g] = obj_Phi_alt(Phimat)
  if nargout < 2
    f = opts.L1(X, B, Phimat, 0);
    f2 = opts.L2(Yl, W, Phimat(:,1:tl), 0);    
    f3 = opts.reg_loss(Phimat);    
  else  
    [f,g] = opts.L1(X, B, Phimat, 2);
    [f2,g2] = opts.L2(Yl, W, Phimat(:,1:tl), 2);    
    [f3,g3] = opts.reg_loss(Phimat); 
    g = g + [opts.L2_wgt * g2 zeros(opts.num_basis, tu)] ...
          + opts.reg_wgt * g3;
  end  
  [f,g] = opts.L1(X, B, Phimat, 2);
  [f3,g3] = opts.reg_loss(Phimat); 

  f = f + opts.L2_wgt * f2 + opts.reg_wgt * f3;
end

% Assumes opts.L2_wgt == 0
function [f, g] = obj_Phi_staged(Phimat)
  if nargout < 2
    f = opts.L1(X, B, Phimat, 0);
    f3 = opts.reg_loss(Phimat);    
  else  
    [f,g] = opts.L1(X, B, Phimat, 2);
    [f3,g3] = opts.reg_loss(Phimat); 
    g = g + opts.reg_wgt * g3;
  end
  f = f + opts.reg_wgt * f3;
end

% For Phi objective passed to Nesterov, cannot include regularizer as
% it is implicitly in the solver_Nesterov_generic
% Assumes that L2_wgt ~= 0, as else obj_Phi_noreg_staged used
function [f, g] = obj_Phi_noreg(Phimat)
  [f,g] = opts.L1(X, B, Phimat, 2);
  [f2,g2] = opts.L2(Yl, W, Phimat(:,1:tl), 2);
  f = f + opts.L2_wgt * f2;
  g = g + [opts.L2_wgt * g2 zeros(opts.num_basis, tu)];
end

function [f, g] = obj_Phi_noreg_staged(Phimat)
  if nargout > 1
    [f,g] = opts.L1(X, B, Phimat, 2);
  else
    f = opts.L1(X, B, Phimat, 0);
  end
end

function [f,g] = obj_B(B)
  if nargout < 2, 
    f = opts.L1(X, B, Phi, 0);
  else
    [f,g] = opts.L1(X, B, Phi, 1);
  end
end

function [f,g] = obj_W(W)
  if nargout < 2, 
    f = opts.L2(Yl, W, Phi(:, 1:tl), 0);
  else
    [f,g] = opts.L2(Yl, W, Phi(:, 1:tl), 1);
  end
end

end  
  
