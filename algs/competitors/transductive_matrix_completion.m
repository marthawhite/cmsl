function [X_recover Y_recover obj runtime B W Phi] = transductive_matrix_completion(X, Yl, opts)
% Runs TMC; see below function for description
% opts.reg_wgt = regularization parameter
% Note: Currently only works for semi-supervised learning
% not multi-view or single-view subspace learning

  t = size(X,2);
  tl = size(Yl,2);
  tu = t-tl;
  x_dim = size(X,1);
  y_dim = size(Yl,1);
  
  DEFAULTS.eta = 0.25;
  DEFAULTS.L1 = @smooth_L11;
  DEFAULTS.L2 = @smooth_L11;
  DEFAULTS.maxiter = 300;
  DEFAULTS.method = 1;
  DEFAULTS.mu = 1e-5;
  DEFAULTS.Oy = zeros(y_dim,t);
  DEFAULTS.Ox = ones(x_dim,t);
  DEFAULTS.reg_wgt = 0.1; 
  DEFAULTS.verbose = VerboseConst.NONE;

  if nargin < 3
    opts = DEFAULTS;
  else  
    opts = getOptions(opts, DEFAULTS);
  end

  if opts.verbose >= VerboseConst.BASIC_ALG, fprintf(1,'\n\ntransductive_matrix_completion -> Starting...\n\n'); end

  % Run TMC
  start_time = cputime;
  [Z, fval, b] = TMC([[Yl zeros(y_dim,tu)]; X], y_dim, opts);

  % Recover X_0 and Y_0
  B = repmat(b,1,t);
  Y_recover = sign(Z(1:y_dim,:)+B);
  X_recover = Z((y_dim+1):end,:);

  runtime = cputime - start_time;

  % This is silly, but computing objective that alternator and staged
  % used, kinda
  obj = opts.L1(X,X_recover) + opts.L2_wgt*opts.L2(Yl,Y_recover(:,1:tl)) + opts.reg_wgt*trace_norm(Z);

  % If nargout is greater than 4, user requesting models, so do recovery
  if nargout > 4
    [B W Phi] = recover_TMC(Z, y_dim);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMC algorithm    
function [Z, fval, b] = TMC(Z, y_dim, params)
% Transductive Matrix Completion

% Input:     Z: initial matrix ([Y;X]), matrix,
%            y_dim: # of classes,           scalar,
%       reg_wgt: regularization,         scalar,
%       params: parameters, structure;
%               Ox: mask of features, index,        default: [];
%               Oy: mask of labels,   index,        default: [];
%               mu: regularization,   scalar,       default: 1e-5;
%               eta: annealing,       scalar,       default: 1/4;
%               maxiter: max # of iteration,        default: 100;
%               method: w/o bias,     boolean,      default: 1.
%
% Output:    Z: completed matrix (hopefully),  matrix,
%         fval: objective values,              vector,
%            b: bias,                          vector,
% Reference: Transduction with Matrix Completion: Three Birds with One Stone
%                Andrew Goldberg et. al.
%
% Implemented by Yao-Liang Yu
% Enviorentment: Matlab 2010a
% Date: 2011-01-30

  Y = Z(1:y_dim,:);
  X = Z(y_dim+1:end,:);
  noy = sum(opts.Oy(:));
  nox = sum(opts.Ox(:));


  [U,sigma,V] = svds(Z,1);
  Z = sigma*U*V';
  b = zeros(y_dim,1);
  mu0 = sigma;
  L = ceil(log(mu0/opts.mu)/log(1/opts.eta));
  taub = 3.8*noy/opts.reg_wgt/t;
  tauz = min(3.8*noy/opts.reg_wgt, nox);
  tol = 1e-5;

  fval = zeros(L*opts.maxiter+1,1);
  fval(1) = inf;
  num = 1;
  for ii = 1:L
    mu0 = max(mu0*opts.eta, opts.mu);
    for jj = 1:opts.maxiter  
      num = num + 1;
      if opts.method == 0
        Z(end,:) = 1; % projection
      end
      [gZ, f] = grad(Y, X, Z, b, opts.Ox, opts.Oy, opts.reg_wgt, noy, nox);
      A = Z - tauz * gZ;
      [U, Reg_Wgt, V] = svd(A,'econ');
      dd = max(diag(Reg_Wgt)-tauz*mu0, 0);
      Z = U * diag(dd) * V';
      if opts.method == 1
        gb = sum(gZ(1:y_dim,:),2);
        b = b - taub * gb;
      end
      fval(num) = f + mu0*sum(dd);
      if abs(fval(num) - fval(num-1)) <= tol
        break;
      end
    end
  end
  fval = fval(1:num);
  if ii == opts.maxiter
    warning('TMC -> Optimization hit max iterations %u, terminated with tolerance %g\n',...
            opts.maxiter, abs(fval(num) - fval(num-1)));    
  end
end

function [gZ, f] = grad(Y, X, Z, b, Ox, Oy, reg_wgt, noy, nox)
% compute the complete grad, then apply the masks
%     not efficient but works well for Matlab
	if ~isempty(b)
    Z(1:y_dim,:) = Z(1:y_dim,:) + repmat(b,1,size(Y,2));
 end
 tmp = exp(-Y.*Z(1:y_dim,:));
 f = sum(sum(Oy .* log(1+tmp)))*reg_wgt/noy;
 gZ(1:y_dim,:) = Oy .* (-Y .* tmp ./ (tmp+1)) * reg_wgt/noy;
 tmp = Z((1+y_dim):end,:) - X;
 f = f + 1/2/nox*norm(Ox.*tmp,'fro');
 gZ= [gZ; Ox .* tmp / nox]; 
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMC recovery algorithm

function [B, W, Phi] = recover_TMC(Z, y_dim)
% recover the basis, prediction weight, and features for
% transductive matrix completion

% input: Z: completion returned by TMC.m
%        y_dim: size(Y,1), # of classes

  [U,S,V] = svd(Z, 'econ');
  % recall that in TMC, Z=[Y;X]
  W = U(1:y_dim,:);
  B = U(1+y_dim:end,:);
  Phi = S*V';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
