function [f,g] = logistic_loss(Y, U, W, var)

% When nargin == 4, var says for which variable to return gradient
% var = 0 means no gradient
% var = 1 means gradient in U
% var = 2 means gradient in W
% Note for fmin_LBFGS to be used, g must be
% a column vector, whereas this function return row g

denom = numel(Y);

if nargin == 2  % now U is the output of the classifier   
                % used for global optimization
    
    f = sum(sum(log(1 + exp(-Y .* U)))) / denom;
    if nargout > 1
        g = Y ./ (1 + exp(Y .* U)) / (-denom);   % a subgradient in U
    end
elseif nargin == 4  % now U*W is the output of the classifier
        
    F = U * W;
    f = sum(sum(log(1 + exp(-Y .* F)))) / denom;
    
    if var == 0         % Just want the function value
        g = 0;
    elseif var == 1     % Compute the subgradient in U
        
        M = Y ./ (1 + exp(Y .* F)) / (-denom);
        g = M * W';
    elseif var == 2     % Compute the subgradient in W 
        
        M = Y ./ (1 + exp(Y .* F)) / (-denom);
        g = U' * M;
    else
        error('var in logistic_loss can only be 0, 1, or 2.');
    end
else
    error('logistic_loss takes 2 or 4 arguments only.');
end


end
