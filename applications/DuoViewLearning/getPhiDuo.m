function [Phi, X_recons, Y_recons] = getPhiDuo(X, Y, B, W, opts)
% opts.L
% opts.alpha
% If only X and B and W are provided, then finds Phi using only X and
% then computes Y_recons using this reconstruction

	lbfgsb_param = [];
	lbfgsb_param.maxIter = 200;      % max number of iterations
	lbfgsb_param.maxFnCall = 1000;  % max number of calling the function
	lbfgsb_param.relCha = 1e-5;     % tolerance of constraint satisfaction
	lbfgsb_param.tolPG = 1e-5;      % final objective function accuracy parameter
	lbfgsb_param.m = 50;

    k = size(B, 2);
    tu = size(X,2);
    
    %     opts.alpha = opts.alpha * 80;
    Phi_init = randn(k*tu, 1);
    ub = Inf(k*tu, 1);    lb = -ub;
        
    % If no Y provided, then only reconstruct with X     
    if isempty(Y)
        Phi = lbfgsb(Phi_init, lb, ub, @obj_generateY, [], [], lbfgsb_param);        
    else    
        Phi = lbfgsb(Phi_init, lb, ub, @obj_outofsample, [], [], lbfgsb_param);
    end
    Phi = reshape(Phi, k, tu);    
    X_recons = B*Phi;
    Y_recons = W*Phi;
    
%     phi_init = rand(k, 1);
%     ub = Inf(k, 1);   lb = -ub;
%     X_te_recons = zeros(size(B, 1), tu);
%     Y_te_recons = zeros(size(W, 1), tu);
%     for i = 1 : tu
%       x = X_te(:,i);    y = Y_te(:,i);    
%       [phi, obj, iter, feval, msg] = lbfgsb(phi_init, lb, ub, ...
%                           @(phi)obj_one_sample(phi, x, y), [], [], lbfgsb_param);
%       X_te_recons(:, i) = B*phi;    Y_te_recons(:, i) = W*phi;    


  function [f, g] = obj_outofsample(Phi)
    Phi = reshape(Phi, [numel(Phi)/tu, tu]);
    R = B*Phi;
    [f1 G1] = opts.L(X, R);
    G1 = B'*G1;
    R = W*Phi;
    [f2 G2] = opts.L(Y, R);
    G2 = W'*G2;
    
    norms = sqrt(sum(Phi.^2, 2));
    f3 = sum(norms);
    norms(norms < 1e-5) = inf;
    G3 = Phi ./ repmat(norms, 1, size(Phi, 2));
        
    f = f1 + opts.L2_wgt*f2 + opts.alpha*f3;
    G = G1 + opts.L2_wgt*G2 + opts.alpha*G3;
    g = G(:);
end

  function [f, g] = obj_generateY(Phi)
    Phi = reshape(Phi, [numel(Phi)/tu, tu]);
    R = B*Phi;
    [f1 G1] = opts.L(X, R);
    G1 = B'*G1;
    
    norms = sqrt(sum(Phi.^2, 2));
    f3 = sum(norms);
    norms(norms < 1e-5) = inf;
    G3 = Phi ./ repmat(norms, 1, size(Phi, 2));
        
    f = f1 + opts.alpha*f3;
    G = G1 + opts.alpha*G3;
    g = G(:);
  end

end