function [f,g] = smooth_L11(X, B, Phi, var)
% Computes a smooth version of L11 loss
% 
% When nargin == 2, then B is treated as Xhat and gradient
% returned for Xhat
%
% When nargin == 4, var says for which variable to return gradient
% var = 0 means no gradient
% var = 1 means gradient in B
% var = 2 means gradient in Phi

  
  global sigma_smooth_L11;
  
  sigma = sigma_smooth_L11;

  if nargin == 2  % now B is the output of the classifier 
                    % used for global optimization
    Y = B - X;
    %    idx = find (abs(Y) < sigma);
    idx = abs(Y) < sigma;
    Z = abs(Y) - (sigma/2);
    Z(idx) = Y(idx).^2 * (0.5/sigma);
    
    f = sum(sum(Z)); 
    
    if nargout < 2
      g = 0;
    else  
      g = sign(Y);
      g(idx) = Y(idx)/sigma;
    end      
    
  elseif nargin == 4  % now B*Phi is the output of the classifier
    Y = B * Phi - X;
    idx = abs(Y) < sigma; % find (abs(Y) < sigma);
    Z = abs(Y) - sigma/2;
    Z(idx) = Y(idx).^2 * (0.5/sigma);
    
    f = sum(sum(Z));

    if var == 0 || nargout < 2   % Just want the function value
      g = 0;
    elseif var == 1         % gradient wrt B
      g = sign(Y);
      g(idx) = Y(idx)/sigma;
      g = g * Phi';
    elseif var == 2         % gradient wrt Phi
      g = sign(Y);
      g(idx) = Y(idx)/sigma;
      g = B' * g;
    else
      error('var in smooth_L11 can only be 0, 1, or 2.');
    end
  else
    error('smooth_L11 takes 2 or 4 arguments only.');
  end

  denom = size(X,2);
  f = f / denom;
  g = g / denom;

end
