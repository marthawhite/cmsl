function [f,g] = hinge_loss(Y, U, W, var)

% When nargin == 4, var says for which variable to return gradient
% var = 0 means no gradient
% var = 1 means gradient for U
% var = 2 means gradient for W

denom = size(Y, 1) * size(Y, 2);

if nargin == 2  % now U is the output of the classifier
                % used for global optimization
        
    f = sum(sum(max(0, 1 - Y .* U))) / denom;
    g = (-1.0 / denom) * Y .* (Y .* U < 1);    % a subgradient in U
    
elseif nargin == 4    % U * W is the output of the classifier
    
    F = U * W;    
    denom = size(Y, 1) * size(Y, 2);
    f = sum(sum(max(0, 1 - Y .* F))) / denom;
    
    if var == 0         % Just want the function value
        g = 0;
    elseif var == 1     % Compute the subgradient in U
        
        M = (Y .* F <= 1) .* (-Y) / denom;
        g = M * W';
    elseif var == 2     % Compute the subgradient in W 
        
        M = (Y .* F <= 1) .* (-Y) / denom;
        g = U' * M;
    else
        error('var in hinge_loss can only be 0, 1, or 2.');
    end
    
else
    error('hinge_loss takes 2 or 4 arguments only.');
end

end
