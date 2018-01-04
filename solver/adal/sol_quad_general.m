function [X] = sol_quad_general(L, A, X0, opts)
% minimize over X (given A and reg_wgt):
%       0.5*|| X - A ||_F^2 + reg_wgt * L(X)
%
% Input:
%   A: matrix A
%   reg_wgt: weight on regularizer
%   L: loss function, such that by calling
%         [f, g] = L(X)
%       we get the value of L(X) in f, and a subgradient in g
%   X0: initial value of X
%   opts: options
%
% Output:
%   X: the optimal solution

if nargin < 3
    error('sol_quad_general requires at least 3 arguments: L, A and X0');
end

DEFAULTS.reg_wgt = 0.1;  
DEFAULTS.verbose = 0;
DEFAULTS.maxiter = 100;
DEFAULTS.funTol = 1e-5;
DEFAULTS.m = 50;

if nargin < 4
  opts = DEFAULTS;
else
  opts = getOptions(opts, DEFAULTS);
end

[X, f, iter, msg] = lbfgs(@pobj, X0(:), opts); 
X = reshape(X, size(A));

%if opts.verbose >= VerboseConst.DETAILED_SOLVER
%  cprintf(VerboseConst.DETAILED_SOLVER_COLOUR, ...
%          'sol_quad_general (lbfgs): f = %f, iter = %d, msg = %s\n', ...
%            f, iter, msg2str(msg, 'lbfgs'));
%end

  function [f, g] = pobj(X)
    X = reshape(X, size(A));
    G = X - A;
    if nargout < 2
      f = L(X);
    else 
      [f, g] = L(X);
      g = opts.reg_wgt*g + G;
      g = g(:);
    end
    f = 0.5 * norm(G, 'fro')^2 + opts.reg_wgt * f;
  end

end
