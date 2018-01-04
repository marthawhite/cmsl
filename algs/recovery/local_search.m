function [B, W, Phi] = local_search(X_tr, Y_tr, B_ori, W_ori, Phi_ori, opts)

  DEFAULTS = [];
  DEFAULTS.gamma = 1;
  DEFAULTS.lbfgs_param = [];
  DEFAULTS.recover_batch_x = true;
  DEFAULTS.recover_online_xy = false;
  DEFAULTS.recover_online_x = false;
  DEFAULTS.use_local_search = 1;

  DEFAULTS.reg_wgt = 1e-1;
  DEFAULTS.L2_wgt = 1;
  DEFAULTS.L1 = @smooth_L11; % reconstruction loss on X
  DEFAULTS.L2 = @smooth_L11; % classification or reconstruction loss on Y
  DEFAULTS.verbose = VerboseConst.NONE;      
  DEFAULTS.latent_dim = 5;

  if nargin < 6
    opts = DEFAULTS;
  else
    opts = getOptions(opts, DEFAULTS);
  end
  
  
  B = B_ori;  W = W_ori;  Phi = Phi_ori;
  [k, t] = size(Phi);
  
  prePhi = Phi;
  if opts.verbose >= VerboseConst.BASIC_ALG, cprintf('blue', 'Starting local search\n'); end
  
  % TODO: Why do we have reg_wgt = 0?
	Wopts = opts; Wopts.radius_ball = opts.gamma; Wopts.reg_wgt = 0; Wopts.init_L = 1/t;
	Bopts = opts; Bopts.radius_ball = 1; Bopts.reg_wgt = 0; Bopts.init_L = 1/t;
  
  for i = 1 : 40
    if opts.verbose > 1
      f1 = opts.L1(X_tr, B*Phi);
      f2 = opts.L2(Y_tr, W*Phi);
      obj = f1 + opts.L2_wgt*f2;      
      fprintf('iter: %d, obj = %g\n', i, obj);
    end

    Phivec = lbfgs(@obj_train_Phi, Phi(:), opts.lbfgs_param);
    Phi = reshape(Phivec, [k, t]);
    if norm(prePhi - Phi, 'fro') < 1e-4
      break;
    else
      prePhi = Phi;
    end

    B = solve_Nesterov_generic(@obj_B, [], @prox_op_L2ball, B, Bopts);
    W = solve_Nesterov_generic(@obj_W, [], @prox_op_L2ball, W,  Wopts);
    %[B, loss] = minH_L2_Nesterov(X_tr, Phi, B, opts.L1, 1, (opts.verbose > 1));
    %[W, loss] = minH_L2_Nesterov(Y_tr, Phi, W, opts.L2, opts.gamma, (opts.verbose > 1));
  end
  
  
  function [f, g]= obj_train_Phi(Phivec)
    Phimat = reshape(Phivec, [k, t]);
    [f1 g1] = opts.L1(X_tr, B*Phimat);
    [f2 g2] = opts.L2(Y_tr, W*Phimat);
    f = f1 + opts.L2_wgt*f2;
    g = B'*g1 + opts.L2_wgt*W'*g2;
    g = g(:);
  end
  
  function [f, g] = obj_W(W)
    [f, g] = opts.L1(Y_tr, W, Phi, 1);
  end

  function [f, g] = obj_B(B)
    [f, g] = opts.L2(X_tr, B, Phi, 1);
  end

end
