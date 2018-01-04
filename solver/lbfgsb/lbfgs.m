function [x,f,iter,flag] = lbfgs(fun,x0,opts)
% limited memory BFGS wrapper to make the 
% interface same as the lbfgs in uofa_optimizers
%
  
% Add reasonable iteration parameters
lbfgs_opts = [];
lbfgsb_opts.maxIter = 50;      % max number of iterations
lbfgsb_opts.maxFnCall = 1000;  % max number of calling the function
lbfgsb_opts.relCha = 1e-5;     % tolerance of constraint satisfaction
lbfgsb_opts.tolPG = 1e-5;      % final objective function accuracy parameter
lbfgsb_opts.m = 20;
lb = -Inf(size(x0)); ub = Inf(size(x0));

% The opts for fmin_LBFGS are similar with different names
% Here, we override the above opts with the corresponding
% ones that would be passed to fmin_LBFGS
% If opts has parameters matching above, override them
if nargin > 2
  if isfield(opts, 'maxiter')
    lbfgsb_opts.maxIter = opts.maxiter;
    lbfgsb_opts.maxFnCall = opts.maxiter*50;
  elseif isfield(opts, 'funTol')
    lbfgsb_opts.tolPG = opts.funTol;
  elseif isfield(opts, 'curvTol')
    lbfgsb_opts.relCha = opts.curvTol;
  elseif isfield(opts, 'm')
    lbfgsb_opts.m = opts.m;
  end

  % If A = [1; -1], then b = [lb; -ub]
  % Or A = 1, then b = lb, or A = -1 then b = -ub
  if isfield(opts, 'A') && length(opts.A(:)) <= 2
    ind_lb = find(opts.A == 1);
    if ~isempty(ind_lb), lb = opts.b(ind_lb); end
    ind_ub = find(opts.A == -1);
    if ~isempty(ind_ub), ub = -opts.b(ind_ub); end  
  end
end

[x,f,iter,feval,msg] = lbfgsb(x0, lb, ub, fun, [], [], lbfgsb_opts);
flag = msg;

end

