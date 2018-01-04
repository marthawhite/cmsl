function [f,g] = L21_loss(W)
%
%   f = sum_i || W(i, :) ||_2

    temp = sqrt(sum(W .^2, 2));
    f = sum(temp);
    
    if nargout > 1
        [num_row, num_col] = size(W);

        g = zeros(size(W));
        for cur_row = 1 : num_row        
            if temp(cur_row) > 1e-6
                g(cur_row, :) = W(cur_row, :) / temp(cur_row);
            else
                g(cur_row, :) = zeros(1, num_col);
            end            
        end
    end    
end
