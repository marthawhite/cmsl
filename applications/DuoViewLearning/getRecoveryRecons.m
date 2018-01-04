function [X_tr_recons, Y_tr_recons, X_te_recons, Y_te_recons, B, W, Phi, pobj] = getRecoveryRecons(X_tr, Y_tr, X1_tr, Y1_tr, X_te, Y_te, X1_te, Y1_te, B,W,Phi,opts)

if nargin < 11
    error('getRecoveryRecons requires at least 11 argumets');
end

DEFAULTS = [];

DEFAULTS.recover_batch_xy = false;
DEFAULTS.recover_batch_x = true;
DEFAULTS.recover_online_xy = false;
DEFAULTS.recover_online_x = false;
DEFAULTS.use_local_search = 1;

DEFAULTS.alpha = 1e-1;
DEFAULTS.L2_wgt = 1;
DEFAULTS.gamma = 1;
DEFAULTS.L1 = @smooth_L11;
DEFAULTS.L2 = @smooth_L11;
DEFAULTS.Recover = 1;
DEFAULTS.verbose = 0;       % 0: nothing
                        % 1: message of the outer solver (line search)
                        % 2: message of the inner solver (eg. neseterov)

DEFAULTS.out_var = 'eta';         % 'eta' or 'rho'
DEFAULTS.out_solver = 'pbm';   % 'pbm', 'lbfgsb', or 'fmin_BFGS'
DEFAULTS.in_solver = 'adal';  % 'pbm', 'lbfgsb', 'adal', or 'nesterov'
DEFAULTS.ssl_solver = 'adal';  % 'pbm', 'lbfgsb', 'adal', or 'nesterov'
DEFAULTS.latent_dim = 5;

if nargin < 12
    opts = DEFAULTS;
else
    opts = getOptions(opts, DEFAULTS);
end

k = opts.latent_dim;
tu = size(X_te,2);


lbfgsb_param = [];
lbfgsb_param.maxIter = 200;      % max number of iterations
lbfgsb_param.maxFnCall = 1000;  % max number of calling the function
lbfgsb_param.relCha = 1e-5;     % tolerance of constraint satisfaction
lbfgsb_param.tolPG = 1e-5;      % final objective function accuracy parameter
lbfgsb_param.m = 50;


X_te_recons = []; Y_te_recons = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if k < 0 || size(B,2) < k
            k = size(B, 2);
        else
            B = B(:, 1:k);     W = W(:, 1:k);     Phi = Phi(1:k, :);
        end
        if opts.use_local_search
    %         local_search1();    % locally adjust B, W, Phi
            [B, W, Phi] = local_search(X_tr, Y_tr, B, W, Phi, lbfgsb_param, opts);
            X_tr_recons = B*Phi;  Y_tr_recons = W*Phi;
            show_result(X_tr_recons, Y_tr_recons, X1_tr, Y1_tr, 'Training with duo view after local search');
        end
        pobj = opts.L1(X_tr, X_tr_recons) + opts.L2(Y_tr, Y_tr_recons) ...
                        + opts.alpha*L21_loss(Phi);
        % Print the noise level


        % Next do out-of-sample prediction: Given a new matrix X Y, solve
        % \min_phi  L1(X, B Phi) + L2_wgt*L2(Y, W Phi) + alpha ||Phi||_2
        % And recover Xhat = B Phi and Yhat = W Phi

        if opts.Recover > 0        
            if opts.recover_batch_xy || opts.recover_batch_x
                [X_te_recons, Y_te_recons, tag] = recover_latent_batch();

            elseif opts.recover_online
                [X_te_recons, Y_te_recons, tag] = recover_latent_online();            
            end

            show_result(0, Y_te_recons, 0, Y1_te, tag);
        end

    function [] = show_result (X_recons, Y_recons, X, Y, msg)
        noise = norm(X_recons - X, 'fro')^2 + norm(Y_recons - Y, 'fro')^2;
        signal = norm(X, 'fro')^2 + norm(Y, 'fro')^2;

		fprintf(1,'%s: snr = %g, error = %g\n', msg, signal / noise, noise);
    end

    function [X_te_recons, Y_te_recons, tag] = recover_latent_batch()

        opts.alpha = opts.alpha / 2;
        Phi_init = randn(k*tu, 1);
        ub = Inf(k*tu, 1);    lb = -ub;
        
        if opts.recover_batch_xy
            Phi = lbfgsb(Phi_init, lb, ub, @obj_outofsample_xy, [], [], lbfgsb_param);
            tag = 'Testing on Y using batch duo view';
        elseif opts.recover_batch_x
            Phi = lbfgsb(Phi_init, lb, ub, @obj_outofsample_x, [], [], lbfgsb_param);
            tag = 'Testing on Y using batch X view';
        end
        
        Phi = reshape(Phi, [k, tu]);
        X_te_recons = B*Phi;        Y_te_recons = W*Phi;
    end

    function [X_te_recons, Y_te_recons] = recover_latent_online()
        phi_init = rand(k, 1);   ub = Inf(k, 1);   lb = -ub;
        X_te_recons = zeros(nx, tu);        Y_te_recons = zeros(ny, tu);
        for i = 1 : tu
          x = X_te(:,i);    y = Y_te(:,i);  
          if opts.recover_online_xy
            [phi, obj, iter, feval, msg] = lbfgsb(phi_init, lb, ub, ...
                              @(phi)obj_online_xy(phi, x, y), [], [], lbfgsb_param);
          elseif opts.recover_online_x
            [phi, obj, iter, feval, msg] = lbfgsb(phi_init, lb, ub, ...
                              @(phi)obj_online_x(phi, x), [], [], lbfgsb_param);
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
        f = f1 + opts.L2_wgt*f2 + opts.alpha*f3;
        g = g1 + opts.L2_wgt*g2 + opts.alpha*g3;
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
        f = f1 + opts.alpha*f3;
        g = g1 + opts.alpha*g3;
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

        f = f1 + opts.L2_wgt*f2 + opts.alpha*f3;
        G = G1 + opts.L2_wgt*G2 + opts.alpha*G3;
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

        f = f1 + opts.alpha*f3;
        G = G1 + opts.alpha*G3;
        g = G(:);
    end

    function [] = local_search1()
        ub = Inf(k*tl, 1);    lb = -ub;
        prePhi = Phi;
        cprintf('blue', 'Starting local search\n');
        for i = 1 : 100
            f1 = opts.L1(X_tr, B*Phi);
            f2 = opts.L2(Y_tr, W*Phi);
            obj = f1 + opts.L2_wgt*f2;
            fprintf('iter: %d, obj = %g\n', i, obj);
            
            Phivec = lbfgsb(Phi(:), lb, ub, @obj_train_Phi, [], [], lbfgsb_param);
            Phi = reshape(Phivec, [k, tl]);
            if norm(prePhi - Phi, 'fro') < 1e-4
                break;
            else
                prePhi = Phi;
            end
        
            [B, loss] = minH_L2_Nesterov(X_tr, Phi, B, opts.L1, 1, opts.verbose);
            [W, loss] = minH_L2_Nesterov(Y_tr, Phi, W, opts.L2, opts.gamma, opts.verbose);
        end
    end


    function [f, g]= obj_train_Phi(Phivec)
        Phimat = reshape(Phivec, [k, tl]);
        [f1 g1] = opts.L1(X_tr, B*Phimat);
        [f2 g2] = opts.L2(Y_tr, W*Phimat);
        f = f1 + opts.L2_wgt*f2;
        g = vec(B'*g1 + opts.L2_wgt*W'*g2);
    end

end


