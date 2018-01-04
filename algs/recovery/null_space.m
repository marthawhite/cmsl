% Find a basis of the null space of M (stored in the columns of G)
% gt_one, if set true, forces to return at least one basis 
function G = null_space(M, gt_one)
  [Q R] = qr(M);
  row_sum = sum(R.^2, 2);
  G = Q(:, row_sum < 1e-20);
  
  if nargin < 2 || gt_one    
    if size(G, 2) < 1
        %warning('manually add a dimension in null space');
      [foo idx] = min(row_sum);
      G = Q(:, idx);
    end
  end
end
