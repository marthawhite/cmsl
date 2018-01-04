function [Phi,flag] = getRepresentation(X, Y, B, W, opts)
% Obtain representation for fixed H
%
%    minimize_Phi   L1(X,B*Phi) + beta*reg_loss(Phi_i:)
%
%   Example: L1 is the l2 loss and regularizer_loss is l1 loss
%
%   X has each column as a sample
%
%   Requires that all losses also return a gradient, though it is
%   optional to supply the losses
%   Default: L1 = l2 loss, reg_loss = 2-1 norm
%


t = size(X,2);
flag = 0;

lbfgs_params.funTol = 1e-6;
lbfgs_params.m = 100;
lbfgs_params.maxiter = 1000;
optimizer = @(fun,x)(lbfgs(fun,x,lbfgs_params));

% Initialize Phi (according to Lee&Ng's work)
Phi = rand(size(B,2), t) - 0.5;

[Phivec,f,flag,iter1] = optimizer(@(Phi)(loss1(X,Y,B,W,Phi)), Phi(:));
Phi = reshape(Phivec, size(B,2), t);

  function [f,g] = loss1(X,Y,B,W,Phivec)
    Phimat = reshape(Phivec,size(B,2),t);
    [f1,g1] = opts.L1(X,B,Phimat,2);
    [f3,g3] = opts.reg_loss(Phimat); 
    if ~isempty(Y)
      [f2,g2] = opts.L2(Y,W,Phimat,2);      
      f = f1+opts.L2_wgt*f2+opts.reg_wgt*f3;
      g = g1+opts.L2_wgt*g2+opts.reg_wgt*g3;
    else
      f = f1+opts.reg_wgt*f3;
      g = g1+opts.reg_wgt*g3;      
    end
    g = g(:);
  end   

end  
  
