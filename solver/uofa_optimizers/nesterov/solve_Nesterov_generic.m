function [X, obj, iter, msg] = solve_Nesterov_generic(loss_fun, reg_fun, prox_op, X0, opts)
%    min_X   fun(X) + reg_wgt*reg(X)
% 
%   fun is the loss function, which returns a function value
%   and the gradient (as a matrix, no vectorized)
% 
% We use Nesterov's method with infinity memory

if nargin < 4
  error('solve_Nesterov_general requires at least 4 arguments');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS STARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sparse coding parameters
DEFAULTS.init_L = 1 / size(X0,2);  % estimate of Lipschitz constant, adjusted during opt
DEFAULTS.reg_wgt = 1;  % weight on regularizer on Phi
DEFAULTS.maxiter = 100; % previously 200, changed for time issues
DEFAULTS.funTol = 1e-6;
DEFAULTS.verbose = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    opts = DEFAULTS;
else
    opts = getOptions(opts, DEFAULTS);
end

if isempty(reg_fun)
  reg_fun = @const_zero;
end

% Initialize all variables for optimization
X = X0;
L = opts.init_L;
R = 0.0;
G = X0;
c = norm(G, 'fro')^2/2;
Z = X;
gamma_d = 1.5;
gamma_u = 1.2;
original_reg_wgt = opts.reg_wgt;

fval = zeros(opts.maxiter, 1);
opts.weight_loss = 1;
f_pre = loss_fun(X0) + opts.reg_wgt*reg_fun(X0);


if opts.verbose > VerboseConst.DETAILED_SOLVER, fprintf('Starting Nesterov for H\n'); end
msg = 'Nesterov_generic quits with max number of iteration';

for iter = 1 : opts.maxiter

  L = L / (gamma_d * gamma_u);
  
  while 1
    
    L = L * gamma_u;
    if L > 100000
      obj = JHkp1; % Some default value, better than crashing or hanging
      return;
    end            
    a = (1 + sqrt(1+4*L*R))/(2*L);
    Rkp1 = R + a;
    U = (a / Rkp1) * Z + (R / Rkp1) * X;
    
    [f, g] = loss_fun(U);
    
    Gkp1 = G - a * g;
    ckp1 = c + a * f - a*sum(sum(g.*U));
    opts.reg_wgt = original_reg_wgt * Rkp1;
    Zkp1 = prox_op(Gkp1, opts);
    
    Xkp1 = (R / Rkp1) * X + (a / Rkp1) * Zkp1;
    JHkp1 = loss_fun(Xkp1) + original_reg_wgt*reg_fun(Xkp1);
    phikp1 = 0.5 * sum(sum(Zkp1 .* (Zkp1 - 2*Gkp1))) + ckp1 ...
             + opts.reg_wgt*reg_fun(Zkp1);
    
    if Rkp1 * JHkp1 <= phikp1
      break; 
    end
  end
  
  G = Gkp1;
  c = ckp1;
  R = Rkp1;
  X = Xkp1;
  Z = Zkp1;
  
  fval(iter) = JHkp1;
  %if opts.verbose >= VerboseConst.DETAILED_SOLVER, fprintf('%d: %g\n', i, JHkp1); end
  
  if abs(f_pre - JHkp1)/min(f_pre, JHkp1) < opts.funTol
    msg = ['Nesterov_generic quits with small relative difference: '... 
           num2str(abs(f_pre - JHkp1)/min(f_pre, JHkp1))];
    break;
  end
  
  f_pre = JHkp1;    
end

%if opts.verbose >= VerboseConst.DETAILED_SOLVER, fprintf('%s\n', msg); end

obj = JHkp1;


function [f] = const_zero(x)
  f = 0;
end

end
