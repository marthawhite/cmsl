function [U Phi] = cone_recover (Z, dims, opts)
% Given Z and gamma, recover the U and Phi
% Z is sized (dims(1) + dims(2)) * dims(3)
%
% Author: Yaoliag Yu, University of Alberta, 2012
  
  if nargin < 3, error('cone_recover needs at least 3 arguments.\n'); end

  % the matrix Zhat is sized (n+c)*t
  n = dims(1);
  c = dims(2);
  t = dims(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS STARTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  DEFAULTS.gamma = 1;
  DEFAULTS.verbose = VerboseConst.NONE;
  DEFAULTS.latent_dim = n+c;
  DEFAULTS.diff_thrd = 1e-2;       % exiting threshold of relative error in Frobenius norm
  DEFAULTS.ls_init_dtnorm = 1;     % initial value of line search in dual triple norm
  DEFAULTS.ls_init_oracle = 0;     % initial value of line search in oracle
  DEFAULTS.max_boost_iter = 100;   % maximum number of boosting iterations
  DEFAULTS.max_power_iter = 200;   % maximum number of iterations in each call of power iteration
  DEFAULTS.max_try_power_iter = 0; % in how many boosting iterations do we try 
                                   % using power iteration for strong oracle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFAULT PARAMETERS ENDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nargin < 3
    opts = DEFAULTS;
  else  
    opts = getOptions(opts, DEFAULTS);
  end

  if norm(Z, 'fro') < 1e-5
    warning('cone_recover -> norm of Z was too small, randomly setting U and setting Phi = 0!\n');
    x = randn(n, 1);   y = randn(c, 1);
    U = [x/norm(x); y/norm(y)];
    Phi = zeros(1, t);
    return;
  end

  g2 = opts.gamma^2;

  [pm_bfgs_dtnorm, pm_pbm_dtnorm, pm_bfgs_oracle, pm_pbm_oracle,pm_bfgs_tcup] ...
      = get_opt_param ();

  % Try two solvers for line search in computing the dual triple nrom
  dual_tnorm_lbfgs = inf; dual_tnorm_pbm = inf;

  max_rho = 100;
  pm_bfgs_dtnorm.A = [1; -1];
  pm_bfgs_dtnorm.b = [0; -max_rho];
  [rho_lbfgs, dual_tnorm_lbfgs] = lbfgs(@(rho)(dual_norm_rho(rho,Z)), ... 
                                        opts.ls_init_dtnorm, pm_bfgs_dtnorm);

  [rho_pbm, dual_tnorm_pbm, iter, feval, flag] = ...
      pbm(@(rho)(dual_norm_rho(rho,Z)), opts.ls_init_dtnorm, ...
          1e-10, inf, pm_pbm_dtnorm);  

  % Among all the solvers, pick the one giving the best result
  if dual_tnorm_lbfgs < dual_tnorm_pbm
    dual_tnorm = -dual_tnorm_lbfgs;   rho = rho_lbfgs;
  else
    dual_tnorm = -dual_tnorm_pbm;  rho = rho_pbm;
  end

  if rho < 1e-5, rho = 0; end

  Z_frob = norm(Z, 'fro');
  if opts.verbose >= VerboseConst.BASIC_ALG
    fprintf(1, 'dual triple norm = %f, ||Z||_F = %f, rho = %g, gamma = %g\n', ...
            dual_tnorm, Z_frob, rho, opts.gamma);
  end

  % Recover Lambda, the associated null space, and 
  [U,~,V] = svd(D_rho_inv(rho)*Z, 'econ');
  Lambda = D_rho_inv(rho)*U*V';
  lam = 1/(1+g2*rho);
  nu = lam*rho;
  M = diag([lam*ones(n,1); nu*ones(c,1)]) - Lambda * Lambda';

  G = null_space(M);
  Gx = G(1:n,:);
  Gy = G(n+1:end,:);
  [P Sigma] = eig(g2*(Gx'*Gx) - Gy'*Gy);
  Sigma = diag(Sigma);
  k = size(G, 2);

  % Starting main boosting loop
  K = zeros(k,k);
  V = [];
  mu = [];
  GP = G*P;
  AtZBt = GP' * Z * Lambda' * GP; % pre-compute A' Z B'
  temp = Lambda'*GP;
  BBt = temp' * temp;
  idx_pos = find(Sigma >= 0); idx_neg = find(Sigma < 0);
  Sigma_pos = Sigma(idx_pos); Sigma_neg = Sigma(idx_neg);
  time_all = 0;
  time_oracle = 0;
  try_power = opts.max_try_power_iter > 0;

  for boost_iter = 1:opts.max_boost_iter   % boosting
    
    start_time = cputime;
    grad = K*BBt - AtZBt;
    grad = (-0.5)*(grad+grad');

    % get a new basis.  First try power iteration. 
    if try_power
      [beta, diff] = solve_sep_prob_power(grad, Sigma_pos, Sigma_neg, ...
                                          idx_pos, idx_neg, opts.max_power_iter);
    else
      diff = 1;
    end
    if diff > 4e-3  % if power iteratino fails
      if opts.verbose >= VerboseConst.DETAILED_ALG
        if try_power
          cprintf('r', 'switching to line search as power iteration failed\n');
        else
          fprintf(1, 'using strict line search\n')
        end
      end
      beta = solve_sep_prob_ls(grad, Sigma, opts.ls_init_oracle, pm_pbm_oracle);
    end
    time_oracle = time_oracle + (cputime - start_time);

    % Now totally corrective update on mu.  Precompute some matrices for QP. 
    V = [V, beta];
    S = GP*V;
    T = Lambda'*S;
    S2 = S'*S;
    T2 = T'*T;
    F = S2.*T2;
    w = sum((Z*T).*S)';
    
    mu0 = [mu; 1];
    lb = zeros(length(mu0), 1);
    pm_bfgs_tcup.A = 1;
    pm_bfgs_tcup.b = lb;  
    [mu_opt, f, flag, iter] = lbfgs(@total_correct, mu0, pm_bfgs_tcup);

    idx = find(mu_opt > 1e-5);  
    V = V(:,idx);
    mu = mu_opt(idx);
    
    if ~isempty(mu)
      K = (repmat(mu', k, 1) .* V) * V';
    else
      K = zeros(k, k);
    end
    time_all = time_all + (cputime - start_time);
    
    % print discrepancy  
    U = GP * V * sqrt(1+g2);
    Phi = repmat(mu/(1+g2), 1, t) .* (U'*Lambda);
    diff_fro = norm(U*Phi - Z, 'fro');  
    rel_fro = diff_fro * 100 / Z_frob;
    diff_tr = sum(sqrt(sum(Phi.^2, 2))) - dual_tnorm;
    rel_tr = abs(diff_tr) * 100 / dual_tnorm;
    if opts.verbose >= VerboseConst.DETAILED_ALG
      fprintf(1, '%d: %d, Frob = %.3f, L21 = %.3f, t_all = %.2f, t_orc = %.2f\n', ...
              boost_iter, length(mu), diff_fro, diff_tr, time_all, time_oracle);
      fprintf(1, '\t  Relative diff: Frob = %.3f %%, L21 = %.3f %%, min_norm = %e\n', ...
              rel_fro, rel_tr, min(mu_opt));  
    end
    if rel_fro < 5 || boost_iter > opts.max_try_power_iter
      try_power = 0;
    end
    
    if rel_fro < opts.diff_thrd*100 || opts.latent_dim < length(mu), break;  end  
  end

  [mu_sort, mu_idx] = sort(mu, 'descend');
  U = U(:, mu_idx);
  Phi = Phi(mu_idx, :);

  function res = D_rho_inv(rho)
    if rho < 1e-5
      res = diag([ones(n, 1); zeros(c,1)]);
    else
      res = diag([(1/sqrt(1+g2*rho))*ones(n, 1); (1/sqrt(g2+1/rho))*ones(c,1)]);
    end
  end

  function [f,g] = dual_norm_rho(rho, Z)
    [f, G] = trace_norm([1.0/sqrt(rho*g2+1) * Z(1:n,:); ...
                        sqrt(rho)/sqrt(g2*rho+1) * Z(n+1:end,:)]);

    g = sum(sum(G .* [(g2 /(-2*(1+g2*rho)^(3/2))) * Z(1:n,:); ...
                      1/(2*(g2*rho+1)^(3/2)*sqrt(rho)) * Z(n+1:end,:)]));
    f = -f;    g = -g;
  end

  function [f, g] = total_correct(mu)    
    g = F*mu;
    f = 0.5*mu'*g - mu'*w;
    g = g - w;
  end
end
