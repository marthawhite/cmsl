function [x,f,flag,varargout] = backtrack(fun,x0,f0,dir,slope,opts)
% backtrack line search
% varargout contains g, constraint_opts, and iter

if nargin < 6
    error('Must supply 6 arguments to backtrack for fmin_LBFGS.');
end    

flag = BFGS_ERRORS.SUCCESS;
varargout = {};
opts.disp = 0;

if any(imag(x0)), x=x0;f=f0;varargout{1}=-dir;flag=BFGS_ERRORS.BACKTRACK_ERROR;return, end;
if any(imag(f0)), x=x0;f=f0;varargout{1}=-dir;flag=BFGS_ERRORS.BACKTRACK_ERROR;return, end;
if any(imag(dir)), x=x0;f=f0;varargout{1}=-dir;flag=BFGS_ERRORS.BACKTRACK_ERROR;return, end;

% linear inequalities
mode_ineq = 0;
if size(opts.A,1) > 0
	mode_ineq = 1; 
	if min(opts.A*x0 - opts.b) < opts.backtrack_funTol
		%warning('infeasible x0 given to backtrack: violates linear inequalities');
		flag =  BFGS_ERRORS.CONSTRAINT;
		x = x0;
		f = f0;
		varargout = cell(1,nargout-3);
		return;
	end
end

% PSD inequality
mode_psd = 0;
if size(opts.Apd,1) > 0
	mode_psd = 1; 
	s = sqrt(length(opts.bpd));
	if s > 0
		if issparse(opts.Apd)
			x0 = sparse(x0);
		end
		M0 = reshape(opts.Apd*x0-opts.bpd,s,s);
		M0 = (M0+M0')/2;
		[R,p] = chol(M0);
		if p
			%warning('infeasible x0 given to backtrack: Matrix(Apd*x0-bpd) not posdef');
			flag = BFGS_ERRORS.CONSTRAINT;
			x = x0;
			f = f0;
			varargout = cell(1,nargout-3);
			return
		end
	end
end

% backtrack
alpha = opts.backtrack_init_stepsize;
runtime = 0; starttime = cputime;
for iter = 1:opts.backtrack_maxiter

	x = x0 + alpha*dir;
  
  runtime = cputime - starttime;
  if opts.backtrack_timeout > 0 && runtime > opts.backtrack_timeout, break; end
  
	if mode_ineq && min(opts.A*x - opts.b) < opts.backtrack_funTol
		alpha = alpha*opts.backtrack_backoff; 
		continue; 
	end

	if mode_psd
		M = reshape(opts.Apd*sparse(x)-opts.bpd,s,s);
		M = (M+M')/2;
		[R,p] = chol(M);
		if p
			alpha = alpha*opts.backtrack_backoff; 
			continue;
		end
	end

	f = fun(x);
	if imag(f) % cholesky says psd is ok, but its not, so keep backtracking
		alpha = alpha*opts.backtrack_backoff;
		continue
	end

  %if length(x0) == 1, f0, f, alpha, slope, end
	if f < f0 || f < f0 + opts.backtrack_acceptfrac*alpha*slope 
		break; 
	end
	alpha = alpha*opts.backtrack_backoff;
end

% if backtrack_maxiters, check feasibility, restore if necessary
if iter == opts.backtrack_maxiter || (opts.backtrack_timeout > 0 && runtime > opts.backtrack_timeout)
	flag = BFGS_ERRORS.BACKTRACK_ITER;
    
	if mode_ineq && min(opts.A*x - opts.b) < opts.backtrack_funTol
		x = x0;
	end

	if mode_psd
		M = reshape(opts.Apd*x-opts.bpd,s,s);
		[R,p] = chol(M);
		if p, x = x0; end
	end
end

if nargout == 4
	[f,varargout{1}] = fun(x);
elseif nargout >= 5
	[f,varargout{1},varargout{2}] = fun(x);
end
if nargout >= 6
	varargout{3} = iter; 
end

