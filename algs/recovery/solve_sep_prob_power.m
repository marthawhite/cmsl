% Use power iteration for strong oracle.

function [basis diff] = solve_sep_prob_power(A, Sigma_pos, Sigma_neg, ...
                                              idx_pos, idx_neg, max_iter)
  
  % now use the heuristic
  
  basis = ones(size(A,1), 1);
  basis_new = basis;
  large_count = 0;
  extra = 100;
  
  for i = 1 : max_iter
    hkp1 = A*basis + extra*basis;
    g = hkp1(idx_pos);
    h = hkp1(idx_neg);
    p = g'*g; q = h'*h;
    r = Sigma_pos' * (g.^2);
    s = -(Sigma_neg' * (h.^2));
    alsq = 1/(p+q*r/s);
    besq = alsq * r / s;
    basis_new(idx_pos) = g * sqrt(alsq);
    basis_new(idx_neg) = h * sqrt(besq);
    diff = norm(basis_new - basis);
    basis = basis_new;
%     fprintf('%f\n', diff);
    if diff < 1e-3, break; end
%     if diff > 1.5
%       large_count = large_count + 1; 
%     else
%       large_count = 0;
%     end
%     if large_count > 10 && i > 100
%       break;
%     end
  end
  
  fprintf('iter_power = %d, diff = %f\n', i, diff);
  
end
