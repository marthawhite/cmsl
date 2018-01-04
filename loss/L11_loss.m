function [f,g] = L11_loss(W)
%
%   f = sum_ij | W(i, j) |
    
    f = sum(vec(abs(W)));
    
    if nargout > 1
        g = sign(W);
    end    
end
