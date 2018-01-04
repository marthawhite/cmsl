function [x,f,iter,flag] = fmin_LBFGS(fun,x0,opts)
% limited memory BFGS
% minimizes fun unconstrained if opts.A, opts.b
% not provided; otherwise, can specify opts.A x >= b
% and/or opts.Aeq x = opts.beq.
% for lb and ub, just set A = [1; -1] and b = [lb; -ub]
%
% Author: Dale Schuurmans, University of Alberta
  
if ~isa(fun,'function_handle'), error('improper function handle'); end

% Can also provide opts.Apd, opts.bpd, opts.funTol, which will be used in backtrack
DEFAULTS.curvTol = 1e-3;  % Tolerance for constraint inaccuracy 
DEFAULTS.funTol = 1e-6;   
DEFAULTS.m = 50;          % number gradients in bundle
DEFAULTS.maxiter = 1000;   
DEFAULTS.timeout = -1;    % By default, no timeout, only maxiters
DEFAULTS.verbose = 0;

% Constraint options; empty if unconstrained
DEFAULTS.A = [];
DEFAULTS.b = [];
DEFAULTS.Aeq = [];
DEFAULTS.beq = [];
DEFAULTS.slope_update_fcn = @inequality_update;  % Set based on constrained or unconstrained
DEFAULTS.mu = 1e-1;
DEFAULTS.maxiter_constrained = 12;
DEFAULTS.mupar = 10;

% Options needed for backtrack
DEFAULTS.backtrack_maxiter = 50;   
DEFAULTS.backtrack_funTol = 1e-4;   
DEFAULTS.backtrack_timeout = -1;   
DEFAULTS.backtrack_init_stepsize = 1;   
DEFAULTS.backtrack_backoff = 0.5;    
DEFAULTS.backtrack_acceptfrac = 0.1;
DEFAULTS.Apd = [];
DEFAULTS.bpd = [];

if nargin < 3
  opts = DEFAULTS;
else
  opts = getOptions(opts, DEFAULTS);
end

% Ensure that x0 is a feasible starting point, depending on the constraints
x = x0;
if nargout(fun) == 3
  [~,~,constraint_opts] = fun(x);
  % constrain_opts contains linear constrains A, b and/or Aeq, beq
  opts = getOptions(opts, constraint_opts);
end    

flag = BFGS_ERRORS.SUCCESS;
if ~isempty(opts.A) && ~isempty(opts.Aeq)
	[x,P,flag] = feastartcon(opts.A,opts.b,opts.Aeq,opts.beq,x0);
  opts.slope_update_fcn = @equality_update;
elseif ~isempty(opts.A)
	[x,flag] = feastartineq(opts.A,opts.b,x0);
elseif ~isempty(opts.Aeq)
	[x,flag] = feastarteq(opts.Aeq,opts.beq,x0);
end
if flag ~= BFGS_ERRORS.SUCCESS
	iter = 0; f = NaN; g = zeros(size(x));
	return;
end

% Ensure timeout enforced in backtrack as well
if opts.timeout > 0 && opts.backtrack_timeout < 0 
  opts.backtrack_timeout = opts.timeout; 
end

t = length(x0);
H0 = speye(t);
Rho = zeros(1,opts.m);
Y = zeros(t,opts.m);
S = zeros(t,opts.m);
inds = [];
slope = Inf;

% damped limited memory BFGS method
starttime = cputime;
iter = 0;
if ~isempty(opts.A)
  % log barrier method
  original_fun = fun;
  original_timeout = opts.timeout;
  if opts.timeout > 0, 
    opts.timeout = ceil(opts.timeout / opts.maxiter_constrained);
  end
  for i = 1:opts.maxiter_constrained
    bar = @(x) barlin(x,opts.A,opts.b,opts.mu);
    fun = @(x)(combine_functions(original_fun,bar,x));
    [f,g] = fun(x);
    for niter = 1:opts.maxiter
      do_break = fmin_update();
      if do_break, break; end
    end   
    opts.mu = opts.mu*opts.mupar;
    iter = iter + niter;
    if opts.timeout > 0 && (cputime - starttime > opts.timeout), break; end    
  end  
else
  [f,g] = fun(x);  
  for iter = 1:opts.maxiter
    do_break = fmin_update();
    if do_break, break; end
  end  
end

if iter >= opts.maxiter
	flag = BFGS_ERRORS.MAXITER;
elseif opts.timeout > 0 && (cputime - starttime > opts.timeout)
	flag = BFGS_ERRORS.TIMEOUT;
end

if opts.verbose >= VerboseConst.DETAILED_SOLVER && flag ~= BFGS_ERRORS.SUCCESS
  cprintf(VerboseConst.DETAILED_SOLVER_COLOUR, ...
          msg2str(flag, 'lbfgs', opts, [', with final slope = ' num2str(full(slope))...
                      ', and runtime = ' num2str(cputime - starttime) '\n']));
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPDATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function do_break = fmin_update()
  do_break = 0;
  % compute search direction
	[slope, dir] = opts.slope_update_fcn();
	%if length(x0) == 1, slope,  f, fun(x+dir) , end
  if -slope < opts.funTol, do_break = 1; return; end
  if opts.backtrack_timeout > 0
    opts.backtrack_timeout = max(1, opts.timeout - (cputime - starttime));
  end  
	[xnew,fnew,flag,gnew] = backtrack(fun,x,f,dir,slope,opts);
	if flag ~= BFGS_ERRORS.SUCCESS, do_break = 1; return; end

	% update memory for estimating inverse Hessian
	s = xnew - x;
	y = gnew - g;
	curvature = y'*s;
	if curvature > opts.curvTol	
		rho = 1/curvature;
		if length(inds) < opts.m
			i = length(inds)+1;
			inds = [inds i];
		else
			i = inds(1);
			inds = [inds(2:end) inds(1)];
  end
  Rho(i) = rho;
  Y(:,i) = y;
  S(:,i) = s;
 end

 x = xnew;
 f = fnew;
 g = gnew;

 if opts.timeout > 0 && (cputime - starttime > opts.timeout), do_break = 1; return; end  
end


function [slope, d] = inequality_update()
% compute search direction
  d = invhessmult(-g,Y,S,Rho,H0,inds,opts.m);
  %if length(x0) == 1, d, end
  slope = d'*g;
end

function [slope, d] = equality_update()
% constrained search direction
  B = invhessmult(opts.Aeq',Y,S,Rho,H0,inds,opts.m);
  hg = invhessmult(g,Y,S,Rho,H0,inds,opts.m);
  y = linsolve(opts.Aeq*B,opts.Aeq*(x-hg) - opts.beq,lin_opts);
  d = -hg - B*y;
  d = d - opts.Aeq'*(P*d);	% re-project onto constraint (numerically helpful)
  
  gproj = g - opts.Aeq'*(P*g);
  slope = d'*gproj;
end

function [f,g] = combine_functions(fun1, fun2, x)
  if nargout < 2
    f1 = fun1(x);
    f2 = fun2(x);
  else
    [f1,g1] = fun1(x);
    [f2,g2] = fun2(x);
    g = g1 + g2;
  end
  f = f1 + f2;
end

end