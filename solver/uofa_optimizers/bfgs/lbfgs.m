function [x,f,iter,flag] = lbfgs(fun,x0,opts)
% limited memory BFGS
% minimizes fun unconstrained if opts.A, opts.b
% not provided; otherwise, can specify opts.A x >= b
% and/or opts.Aeq x = opts.beq.
% for lb and ub, just set A = [1; -1] and b = [lb; -ub]
%
% Author: Dale Schuurmans, University of Alberta

  if nargin < 3
    [x,f,iter,flag] = fmin_LBFGS(fun, x0);
  else
    [x,f,iter,flag] = fmin_LBFGS(fun, x0, opts);
  end
    
end

