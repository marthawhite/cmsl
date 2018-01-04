function [f,g] = euclidean_loss(X, B, Phi, var)

% When nargin == 4, var says for which variable to return gradient
% var = 0 means no gradient
% var = 1 means gradient in B
% var = 2 means gradient in Phi

  denom = size(X,2);

  if nargin == 2  % now B is the output of the classifier 
                  % used for global optimization
    if nargout == 2
      g = (1.0 / denom) * (B - X);        % gradient wrt B    
      f = (denom / 2) * sum(sum(g.^2));
    else
      f = (0.5 / denom) * sum(sum((B - X).^2));
      g=0;
    end
    
  elseif nargin == 4  % now B * Phi is the output of the classifier
    
    F = B * Phi;
    f = (0.5 / denom) * sum(sum((X - F).^2));
    
    if var == 0 || nargout < 2   % Just want the function value
      g = 0;
    elseif var == 1         % gradient wrt B
      g = (F - X) * Phi';
    elseif var == 2         % gradient wrt Phi
      g = B' * (F - X);
    else
      error('var in euclidean_loss can only be 0, 1, or 2.');
    end
    
    g = g / denom;    
  else
    error('euclidean_loss takes 2 or 4 arguments only.');
  end
end
