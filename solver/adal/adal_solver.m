function [X, Lambda] = adal_solver(L1_quad_solver, L2_quad_solver, X0, opts)
% minimize over X:
%       L1(X) + reg_wgt * L2(X)
%
% Use alternating direction augmented Lagrangian method
%   1. Relax into
%       L1(X) + reg_wgt * L2(Y) + <Lambda, X - Y> + 1/(2*mu) ||X - Y||^2_F
%   2. Alternating between X, Y, and updating Lambda
%
% Input:
%   L1_quad_solver: a solver which solves (given A and reg_wgt)
%                   0.5*|| X - A ||_F^2 + reg_wgt * L1(X)
%   L2_quad_solver: a solver which solves (given A and reg_wgt)
%                   0.5*|| X - A ||_F^2 + reg_wgt * L2(X)
%   X0:     initial solution
%   opts:   options
%
% Output:
%   X: the optimal solution



if nargin < 3
    error('adal_solver requires at least 3 arguments');
end

DEFAULTS.funTol = 1e-3;  % Termination tolerance of ADAL
DEFAULTS.Lambda = zeros(size(X0));
DEFAULTS.maxiter = 1000;
DEFAULTS.mu = 10;
DEFAULTS.reg_wgt = 0.1;
DEFAULTS.verbose = 0;
DEFAULTS.timeout = -1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
  opts = DEFAULTS;
else
  opts = getOptions(opts, DEFAULTS);
end

X = X0;   % X corresponds to the \Gamma in the note
Y = X0;   % Y corresponds to the \Phi in the note

iter = 0;
reg_wgt = opts.reg_wgt;
Lambda = opts.Lambda;
starttime = cputime;
modno = floor(opts.maxiter / 20.0);

while iter < opts.maxiter
  
  % Initialize with the optimal solution in the last iteration
  opts.reg_wgt = reg_wgt*opts.mu;
  Y = L2_quad_solver(X+opts.mu*Lambda, Y, opts);
  
  opts.reg_wgt = opts.mu;
  X = L1_quad_solver(Y-opts.mu*Lambda, X, opts);
  opts.reg_wgt = reg_wgt;
  
  D = Y - X;
  Lambda = Lambda - D / opts.mu;

  diff = norm(D, 'fro');
  if opts.verbose >= VerboseConst.DETAILED_SOLVER && ...
        (mod(iter, modno) == 0 || diff < opts.funTol)
    cprintf(VerboseConst.DETAILED_SOLVER_COLOUR, 'iter = %d, ADAL diff = %g\n', iter, diff);
  end
  iter = iter + 1;

  % Terminate if the Frobenius norm of X-Y falls below the threshold
  if diff < opts.funTol
    break;
  end
  
  runtime = cputime - starttime;
  if opts.timeout > 0 && runtime > opts.timeout
    if  opts.verbose >= VerboseConst.DETAILED_SOLVER
      cprintf(VerboseConst.DETAILED_SOLVER_COLOUR, ...
            ['adal_solver -> Timeout after %g seconds\n'], runtime);
    end
    break;
  end    
end

end
