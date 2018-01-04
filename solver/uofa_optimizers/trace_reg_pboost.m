function [Xk, obj, iter, msg] = trace_reg_pboost(func, sizeX, opts)
% Solve  Min_X  func(X) + reg_wgt * ||X||_*
    
DEFAULTS.inner_solver_lbfgs_params = struct(...
    'maxiter', 30, ...
    'funTol', 1e-6 ...
    );
DEFAULTS.num_basis = 50; % Initialize rank of U and V
DEFAULTS.maxiter_inner = 30;
DEFAULTS.reg_wgt = 0.1;
DEFAULTS.funTol = 1e-4;
DEFAULTS.timeout = -1;
DEFAULTS.use_local = true;
DEFAULTS.verbose = 0;

if nargin < 3
  opts = DEFAULTS;
else
  opts = getOptions(opts, DEFAULTS);
end
  
k = opts.num_basis;  
m = sizeX(1);
n = sizeX(2);

U = randn(m, k);
V = randn(k, n);

starttime = cputime;
for i = 1 : opts.maxiter_inner
  
  k = size(U, 2);
  
  if opts.use_local      
    nel = k*(m+n);
    
    % local search
    [UV, obj, msg, iter] = lbfgs(@obj_UV, [U(:); V(:)], opts.inner_solver_lbfgs_params);

    U = reshape(UV(1:k*m), [m, k]);
    V = reshape(UV(1+k*m:end), [k, n]);
  else
    obj = obj_UV([U(:); V(:)]);
    msg = 'None';
  end
  
  sum_U = sqrt(sum(U.*U))';
  sum_V = sqrt(sum(V.*V, 2));
  sum_UV = sqrt(sum_U .* sum_V);
  idx = (sum_U > 1e-5) & (sum_V > 1e-5);
  nidx = sum(idx);
  if nidx < length(idx)
    U = U(:, idx);
    V = V(idx, :);
    sum_U = sum_U(idx);
    sum_V = sum_V(idx);
    sum_UV = sum_UV(idx);
  end
  
  if opts.verbose >= VerboseConst.DETAILED_SOLVER
    cprintf(VerboseConst.DETAILED_SOLVER_COLOUR,'iter=%d, obj=%g, k=%d, total_time = %g, msg = %s\n',...
            i, obj, k,cputime-starttime, msg);
  end
  
  if nidx > 0
    Xk = U*V;    sk = 0.5*(sum_U'*sum_U + sum_V'*sum_V);
  else
    Xk = zeros(m, n);   sk = 0;
  end

  if i > 1 && abs(pre_obj-obj) / min(pre_obj, obj) < opts.funTol
    msg = 'Stop with small relative change';
    break;
  else
    pre_obj = obj;
  end
  
  [~, G] = func(Xk);
  [u, ~, v] = svds(G, 1);    v = -v;
  
  % line search
  lbfgs_constrained = opts.inner_solver_lbfgs_params;
  lbfgs_constrained.A = eye(2);
  lbfgs_constrained.b = [0; 0];
  [weights, obj, msg, iter] = lbfgs(@obj_ls, [1; 1], lbfgs_constrained);
  weights = sqrt(weights);
  if nidx > 0
    U = [U .* repmat(weights(1)*(sum_UV./sum_U)', [m, 1]), weights(2) * u];
    V = [V .* repmat(weights(1)*sum_UV./sum_V, [1, n]); weights(2) * v']; 
  else
    U = weights(2) * u;
    V = weights(2) * v'; 
  end
  
  if opts.timeout > 0 && opts.timeout < (cputime-starttime)
    msg = ['Hits maximum alloted time ' num2str(cputime-starttime) ' seconds'];
    break;
  end
end

iter = i;
if i == opts.maxiter_inner
  msg = 'Stop with max iteration';
end

function [f, g] = obj_UV(X)
  r = length(X) / (m+n);
  UU = reshape(X(1:r*m), [m, r]);
  VV = reshape(X(1+r*m:end), [r, n]);
  if nargout < 2
    f = func(UU*VV);
  else
    [f, G] = func(UU*VV);
    g1 = G*VV'+opts.reg_wgt*UU;
    g2 = UU'*G+opts.reg_wgt*VV;
    g = [g1(:); g2(:)];
  end
  f = f + 0.5*opts.reg_wgt*norm(X)^2;  
end

function [f, g] = obj_ls(x)
  Y = x(1)*Xk + (x(2) * u) * v';
  if nargout < 2
    f = func(Y);
  else  
    [f, G] = func(Y);
    g = [sum(sum(G.*Xk))+opts.reg_wgt*sk; u'*G*v+opts.reg_wgt];
  end
  f = f + opts.reg_wgt*(sk*x(1) + x(2));
end

end
