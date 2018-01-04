function [X_tr_recons, Y_tr_recons, X_te_recons, Y_te_recons, B, W, Phi, pobj] = ...
      adjust_recovery(X_tr, Y_tr, X_te, Y_te, B,W,Phi,opts)

  if nargin < 7
    error('getRecoveryRecons requires at least 7 argumets');
  end

  DEFAULTS = [];

  DEFAULTS.recover_batch_xy = false;
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

  if nargin < 8
    opts = DEFAULTS;
  else
    opts = getOptions(opts, DEFAULTS);
  end

  k = opts.latent_dim;
  tu = size(X_te,2);

  
  lbfgs_param = [];
  lbfgs_param.maxiter = 1000;      % max number of iterations
  lbfgs_param.funTol = 1e-4;  % max number of calling the function

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if k < 0 || size(B,2) < k
    k = size(B, 2);
  else
    B = B(:, 1:k);     W = W(:, 1:k);     Phi = Phi(1:k, :);
  end
 
  if opts.use_local_search
    [B, W, Phi] = local_search(X_tr, Y_tr, B, W, Phi, opts);
    X_tr_recons = B*Phi;  Y_tr_recons = W*Phi;
  end
  pobj = opts.L1(X_tr, X_tr_recons) + opts.L2(Y_tr, Y_tr_recons) ...
         + opts.reg_wgt*L21_loss(Phi);

  % Next do out-of-sample prediction: Given a new matrix X Y, solve
  % \min_phi  L1(X, B Phi) + L2_wgt*L2(Y, W Phi) + reg_wgt ||Phi||_2
  % And recover Xhat = B Phi and Yhat = W Phi

  if isempty(X_te)
    X_te_recons = []; Y_te_recons = [];
  
  elseif opts.recover_batch_xy || opts.recover_batch_x
    [X_te_recons, Y_te_recons, tag] = recover_latent_batch();

  elseif opts.recover_online
    [X_te_recons, Y_te_recons, tag] = recover_latent_online();            
  end
  

  function [X_te_recons, Y_te_recons, tag] = recover_latent_batch()

    opts.reg_wgt = opts.reg_wgt / 2;
    Phi_init = randn(k*tu, 1);
    
    if opts.recover_batch_xy
      Phi = lbfgs(@obj_outofsample_xy, Phi_init, lbfgs_param);
      tag = 'Testing on Y using batch duo view';
    elseif opts.recover_batch_x
      Phi = lbfgs(@obj_outofsample_x, Phi_init, lbfgs_param);
      tag = 'Testing on Y using batch X view';
    end
    
    Phi = reshape(Phi, [k, tu]);
    X_te_recons = B*Phi;
    Y_te_recons = W*Phi;
  end

  function [X_te_recons, Y_te_recons] = recover_latent_online()
    phi_init = rand(k, 1);
    X_te_recons = zeros(nx, tu);        
    Y_te_recons = zeros(ny, tu);
    for i = 1 : tu
      x = X_te(:,i);
      y = Y_te(:,i);  
      if opts.recover_online_xy
        [phi, obj, iter, feval, msg] = lbfgs(@(phi)obj_online_xy(phi, x, y), phi_init, lbfgs_param);
      elseif opts.recover_online_x
        [phi, obj, iter, feval, msg] = lbfgs(@(phi)obj_online_x(phi, x), phi_init, lbfgs_param);
      end
      X_te_recons(:, i) = B*phi;    Y_te_recons(:, i) = W*phi;    
    end
    
    Phi = reshape(Phi, k, tu);
    X_te_recons = B*Phi;        Y_te_recons = W*Phi;
    if opts.recover_online_xy
      tag = 'Testing on Y using online duo view';
    elseif opts.recover_online_x
      tag = 'Testing on Y using online X view';
    end
  end

  function [f, g] = obj_online_xy(phi, x, y)
    z = B*phi;
    [f1 g1] = opts.L1(x, z);
    g1 = B'*g1;
    z = W*phi;
    [f2 g2] = opts.L2(y, z);
    g2 = W'*g2;
    f3 = norm(phi);
    if f3 > 1e-5
      g3 = phi / f3;
    else
      g3 = zeros(size(phi));
    end
    f = f1 + opts.L2_wgt*f2 + opts.reg_wgt*f3;
    g = g1 + opts.L2_wgt*g2 + opts.reg_wgt*g3;
  end

  function [f, g] = obj_online_x(phi, x)
    z = B*phi;
    [f1 g1] = opts.L1(x, z);
    g1 = B'*g1;
    f3 = norm(phi);
    if f3 > 1e-5
      g3 = phi / f3;
    else
      g3 = zeros(size(phi));
    end
    f = f1 + opts.reg_wgt*f3;
    g = g1 + opts.reg_wgt*g3;
  end

  function [f, g] = obj_outofsample_xy(Phi)
    Phi = reshape(Phi, [k, tu]);
    R = B*Phi;
    [f1 G1] = opts.L1(X_te, R);
    G1 = B'*G1;
    R = W*Phi;
    [f2 G2] = opts.L2(Y_te, R);
    G2 = W'*G2;

    norms = sqrt(sum(Phi.^2, 2));
    f3 = sum(norms);
    norms(norms < 1e-5) = inf;
    G3 = Phi ./ repmat(norms, 1, size(Phi, 2));

    f = f1 + opts.L2_wgt*f2 + opts.reg_wgt*f3;
    G = G1 + opts.L2_wgt*G2 + opts.reg_wgt*G3;
    g = G(:);
  end

  function [f, g] = obj_outofsample_x(Phi)
    Phi = reshape(Phi, [k, tu]);
    R = B*Phi;
    [f1 G1] = opts.L1(X_te, R);
    G1 = B'*G1;

    norms = sqrt(sum(Phi.^2, 2));
    f3 = sum(norms);
    norms(norms < 1e-5) = inf;
    G3 = Phi ./ repmat(norms, 1, size(Phi, 2));

    f = f1 + opts.reg_wgt*f3;
    G = G1 + opts.reg_wgt*G3;
    g = G(:);
  end

  function [f, g]= obj_train_Phi(Phivec)
    Phimat = reshape(Phivec, [k, tl]);
    [f1 g1] = opts.L1(X_tr, B*Phimat);
    [f2 g2] = opts.L2(Y_tr, W*Phimat);
    f = f1 + opts.L2_wgt*f2;
    g = vec(B'*g1 + opts.L2_wgt*W'*g2);
  end

end


